//
// Created by Jordan Dialpuri on 21/04/2025.
//

#ifndef BACKBONE_TRACING_H
#define BACKBONE_TRACING_H
#include <fstream>

#include "predicted-maps.h"
#include "nucleofind-utils.h"
#include <iostream>
#include <unordered_set>
#include <gemmi/neighbor.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/to_mmcif.hpp>
#include <gemmi/to_pdb.hpp>

namespace NucleoFind {

    struct Edge {
        Edge(int source, int target, double distance): source(source), target(target), distance(distance) {}

        int source;
        int target;
        double distance;
    };

    struct Node {
        explicit Node(int index): point_index(index) {}

        void add_edge(std::shared_ptr<Edge> edge) {
            edges.push_back(edge);
        }

        int degree() const {
            return edges.size();
        }

        int point_index;
        std::vector<std::shared_ptr<Edge>> edges;
    };

    struct BackboneGraph {

        BackboneGraph(gemmi::Residue& input, gemmi::Grid<>& xgrid, double max_distance = 8.0): input(input), xgrid(xgrid), max_distance(max_distance) {
            generate_graph();
            // print_stats();
            create_linked_structure();
        };

    private:

        void generate_graph() {
            nodes.reserve(input.atoms.size());

            gemmi::Model model = create_gemmi_model(input);
            gemmi::Structure s = split_atoms_into_structure(input, xgrid.unit_cell, xgrid.spacegroup);

            gemmi::NeighborSearch ns = gemmi::NeighborSearch(model, xgrid.unit_cell, max_distance).populate();

            for (int a = 0; a < input.atoms.size(); a++) {
                nodes.push_back(std::make_shared<Node>(a));
            }

            for (int a = 0; a < input.atoms.size(); a++) {
                auto nearby = ns.find_atoms(input.atoms[a].pos, '*', 0, max_distance);

                for (const auto& near: nearby) {

                    if (near->atom_idx == a) continue;

                    auto symm_copy_near = symmetry_copy_near(near->pos, input.atoms[a].pos, xgrid.unit_cell, xgrid.spacegroup);
                    // input.atoms[near->atom_idx].pos = symm_copy_near;

                    std::cout << "Distance between " << a << "-" << near->atom_idx << " is " << (input.atoms[a].pos - symm_copy_near).length() << std::endl;

                    auto edge = std::make_shared<Edge>(a, near->atom_idx, (input.atoms[a].pos-near->pos).length());
                    edges.push_back(edge);

                    nodes[a]->add_edge(edge);
                    nodes[near->atom_idx]->add_edge(edge);
                }


            }

            auto fs= split_atoms_into_structure(input, xgrid.unit_cell, xgrid.spacegroup);
            std::ofstream ofs;
            ofs.open("final.cif");
            gemmi::cif::write_cif_to_stream(ofs, make_mmcif_document(s));
            ofs.close();
        }


        const std::vector<std::shared_ptr<Node>>& get_nodes() const {
            return nodes;
        }

        const std::vector<std::shared_ptr<Edge>>& get_edges() const {
            return edges;
        }

        std::shared_ptr<Node> get_node(int index) const {
            if (index >= 0 && index < nodes.size()) {
                return nodes[index];
            }
            return nullptr;
        }

        std::optional<gemmi::Position> get_position(int index) const {
            if (index >= 0 && index < nodes.size()) {
                return positions[index];
            }
            return std::nullopt;
        }

        std::vector<int> find_node_by_degree(int degree) const {
            std::vector<int> result;
            for (int i = 0; i < nodes.size(); i++) {
                if (nodes[i]->degree() == degree) {
                    result.push_back(i);
                }
            }
            return result;
        }

        std::vector<int> find_branch_points() const {
            std::vector<int> result;
            for (int i = 0; i < nodes.size(); i++) {
                if (nodes[i]->degree() > 2) {
                    result.push_back(i);
                }
            }
            return result;
        }

        std::vector<int> find_end_points() const {
            return find_node_by_degree(1);
        }

        std::vector<std::vector<int>> find_connected_components() const {
            std::vector<std::vector<int>> components;
            std::unordered_set<int> visited;

            for (size_t i = 0; i < nodes.size(); ++i) {
                if (visited.count(i) > 0) continue;

                // Start a new component
                std::vector<int> component;
                std::queue<int> queue;

                queue.push(i);
                visited.insert(i);

                while (!queue.empty()) {
                    int current = queue.front();
                    queue.pop();
                    component.push_back(current);

                    // Visit all neighbors
                    for (const auto& edge : nodes[current]->edges) {
                        int neighbor = (edge->source == current) ? edge->target : edge->source;

                        if (visited.count(neighbor) == 0) {
                            visited.insert(neighbor);
                            queue.push(neighbor);
                        }
                    }
                }

                components.push_back(component);
            }

            return components;
        }

        void create_linked_structure() {
            auto s = create_gemmi_structure(input);

            for (auto i = 0; i < edges.size(); i++) {
                gemmi::AtomAddress a1 = {"A", input.seqid, input.name, std::to_string(edges[i]->source)};
                gemmi::AtomAddress a2 = {"A", input.seqid, input.name, std::to_string(edges[i]->target)};
                gemmi::Connection connection = {
                    std::to_string(i),
                    "",
                    gemmi::Connection::Covale,
                    gemmi::Asu::Any,
                    a1,
                    a2,
                    edges[i]->distance
                };
                s.connections.push_back(connection);
                std::cout << "Connecting " << edges[i]->source << " to " << edges[i]->target << " " << edges[i]->distance <<  std::endl;
            }

            std::ofstream of;
            of.open("triplets.pdb");
            gemmi::write_pdb(s, of);
            of.close();
        }


        void print_stats() const {
            std::cout << "Graph Statistics:" << std::endl;
            std::cout << "  Nodes: " << nodes.size() << std::endl;
            std::cout << "  Edges: " << edges.size() << std::endl;

            // Degree distribution
            std::unordered_map<int, int> degreeCount;
            for (const auto& node : nodes) {
                degreeCount[node->degree()]++;
            }

            std::cout << "  Degree distribution:" << std::endl;
            for (const auto& pair : degreeCount) {
                std::cout << "    Degree " << pair.first << ": " << pair.second << " nodes" << std::endl;
            }

            // Connected components
            // std::vector<std::vector<int>> components = find_connected_components();
            // std::cout << "  Connected components: " << components.size() << std::endl;
            // for (size_t i = 0; i < components.size(); ++i) {
            //     std::cout << "    Component " << i << ": " << components[i].size() << " nodes" << std::endl;
            // }

            // Branch points
            std::vector<int> branchPoints = find_branch_points();
            std::cout << "  Branch points: " << branchPoints.size() << std::endl;

            // End points
            std::vector<int> endPoints = find_end_points();
            std::cout << "  End points: " << endPoints.size() << std::endl;
        }


    private:
        gemmi::Residue input;
        gemmi::Grid<> xgrid;
        std::vector<gemmi::Position> positions;
        std::vector<std::shared_ptr<Node>> nodes;
        std::vector<std::shared_ptr<Edge>> edges;
        double max_distance;
    };
}

#endif //BACKBONE_TRACING_H
