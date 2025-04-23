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
#include <clipper/minimol/minimol_utils.h>
#include <gemmi/neighbor.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/to_mmcif.hpp>
#include <gemmi/to_pdb.hpp>

#include "fragment-library.h"

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

        void remove_edge(const std::shared_ptr<Edge>& edge_to_remove) {
            auto new_end = std::remove_if(edges.begin(), edges.end(),
                [&](const std::shared_ptr<Edge>& edge) {
                    return edge == edge_to_remove;
                }
            );
            edges.erase(new_end, edges.end());
        }

        int degree() const {
            return edges.size();
        }

        int point_index;
        std::vector<std::shared_ptr<Edge>> edges;
    };

    struct BackboneTracer {
        BackboneTracer(clipper::MiniMol& mol,
                        clipper::Xmap<float>& xgrid,
                        PredictedMaps& predicted_maps): mol(mol), xgrid(xgrid), predicted_maps(predicted_maps) {
            library = TriNucleotideLibrary("/Users/dialpuri/Downloads/collated.cif");
            initialise_tracer();
        };

    private:

        // Initialisation functions
        void initialise_tracer() {
            input = mol[0][0];
            generate_sugar_nonbond();
            generate_graph();
            generate_adjacency_list();
            // print_stats();
            identify_and_resolve_branches();
            create_linked_structure();

        }

        void generate_sugar_nonbond() {
            clipper::Xmap<float> sugar = *predicted_maps.get_sugar_map();
            clipper::MiniMol sugar_mol = MapToPoints::create_mol_at_gridpoints(sugar, 0.1);
            predicted_sugar_positions = clipper::MAtomNonBond(sugar_mol, 1);
        }

        void generate_adjacency_list() {
            for (const auto& edge : edges) {
                adjacency[edge->source].push_back(edge);
            }
        }


        // Graph generation functions
        void generate_graph();

        void determine_edge(int source_atom, int target_atom, clipper::Coord_orth &current_atom_orth,
                   clipper::Coord_orth &target_atom_orth);

        void find_nearby_nodes(const clipper::MAtomNonBond& nb, int a);


        // Branch resolution functions
        void identify_and_resolve_branches();


        // Model building functions
        double fit_and_score_fragment(int n1, int n2, int n3);

        std::vector<double> fit_and_score_fragment_individually(int n1, int n2, int n3);

        double extract_library_fragment_and_score(int l, std::vector<clipper::Coord_orth> &reference_coords);

        double score_monomers(std::vector<clipper::MMonomer>& monomers);

        std::vector<double> score_monomers_individually(std::vector<clipper::MMonomer>& monomers);

        double score_monomer(clipper::MMonomer& monomer);

        static double score_to_grid(const clipper::Coord_orth& coord, const clipper::Xmap<float>* grid);

        // Graph Structure Utility Functions
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

        std::optional<clipper::Coord_orth> get_position(int index) const {
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

        std::shared_ptr<Node> find_node_by_point_index(int point_index) const {
            for (const auto& node : nodes) {
                if (node->point_index == point_index) {
                    return node;
                }
            }
            return nullptr;
        }

        std::vector<int> find_branch_points() const {
            std::vector<int> result;
            for (int i = 0; i < nodes.size(); i++) {
                if (nodes[i]->degree() > 4) {
                    result.push_back(i);
                }
            }
            return result;
        }

        std::vector<int> find_end_points() const {
            return find_node_by_degree(2);
        }


        // Debugging functions

        void create_linked_structure();

        void print_stats() const;



    private:
        // Parameters
        double max_distance = 8.0;

        // Clipper members
        clipper::MMonomer input;
        clipper::MiniMol mol;
        clipper::Xmap<float> xgrid;
        PredictedMaps& predicted_maps;
        clipper::MAtomNonBond predicted_sugar_positions; // Used to determine if an edge is valid

        // Library members
        TriNucleotideLibrary library;

        // Graph members
        std::vector<clipper::Coord_orth> positions;
        std::vector<std::shared_ptr<Node>> nodes;
        std::vector<std::shared_ptr<Edge>> edges;
        std::unordered_map<int, std::vector<std::shared_ptr<Edge>>> adjacency;
    };


}

#endif //BACKBONE_TRACING_H
