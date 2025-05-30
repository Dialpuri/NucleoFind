//
// Created by Jordan Dialpuri on 21/04/2025.
//

#ifndef BACKBONE_TRACING_H
#define BACKBONE_TRACING_H
#include <fstream>

#include "predicted-maps.h"
#include "nucleofind-utils.h"
#include <set>
#include <vector>

#include <clipper/minimol/minimol_utils.h>
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

        void remove_edge(int n) {
            auto new_end = std::remove_if(edges.begin(), edges.end(),
                [&](const std::shared_ptr<Edge>& edge) {
                    return edge->target == n || edge->source == n;
                }
            );
            edges.erase(new_end, edges.end());
        }

        std::vector<std::shared_ptr<Edge>> find_outgoing_edges() const {
            std::vector<std::shared_ptr<Edge>> outgoing_edges;
            for (auto& edge: edges) {
                if (edge->source == point_index) {
                    outgoing_edges.push_back(edge);
                }
            }
            return outgoing_edges;
        }

        int degree() const {
            return edges.size();
        }

        int point_index;
        std::vector<std::shared_ptr<Edge>> edges;
    };

    struct ScoredMonomer {
        ScoredMonomer() = default;
        ScoredMonomer(double score, clipper::MMonomer& monomer): score(score), monomer(monomer) {}
        double score;
        clipper::MMonomer monomer;
    };

    struct FragmentResult {
        FragmentResult(const FragmentResult& other) {
            scores = other.scores;
            monomers = other.monomers;
            terminus = other.terminus;
        }
        FragmentResult() = default;
        FragmentResult(std::vector<double>& scores, std::vector<clipper::MMonomer>& monomers): scores(scores), monomers(monomers) {};

        void append(FragmentResult& other) {
            scores.insert(scores.end(), other.scores.begin(), other.scores.end()-1);
            monomers.insert(monomers.end(), other.monomers.begin(), other.monomers.end()-1);
        }

        double get_total_score() const {
            return std::accumulate(scores.begin(), scores.end(), 0.0);
        }

        [[nodiscard]] bool empty() const {
            return monomers.empty();
        }

        FragmentResult sort_result();

        void add_terminus() {
            monomers.emplace_back(terminus.monomer);
            scores.emplace_back(terminus.score);
        };

        ScoredMonomer terminus = {};
        std::vector<clipper::MMonomer> monomers;
        std::vector<double> scores;
    };


    struct BackboneTracer {
        BackboneTracer(clipper::MiniMol& mol,
                        clipper::Xmap<float>& xgrid,
                        PredictedMaps& predicted_maps,
                        const std::string& database_path): mol(mol), xgrid(xgrid), predicted_maps(predicted_maps) {
            library = TriNucleotideLibrary(database_path);
            initialise_tracer();
        };

        clipper::MiniMol build() {
            identify_and_resolve_branches();
            // print_stats();
            return build_chains();
        }

    private:

        // Initialisation functions
        void initialise_tracer() {
            input = mol[0][0];
            generate_sugar_nonbond();
            generate_input_nonbond();
            generate_graph();
            generate_adjacency_list();
            // create_linked_structure();
        }


        void generate_sugar_nonbond() {
            clipper::Xmap<float> sugar = *predicted_maps.get_sugar_map();
            clipper::MiniMol sugar_mol = MapToPoints::create_mol_at_gridpoints(sugar, 0.1);
            predicted_sugar_positions = clipper::MAtomNonBond(sugar_mol, 1);
        }

        void generate_input_nonbond() {
            initial_nb = clipper::MAtomNonBond(mol, 2);
        }

        void generate_adjacency_list() {
            for (const auto& edge : edges) {
                adjacency[edge->source].push_back(edge);
            }
        }


        // Graph functions
        void generate_graph();

        void determine_edge(int source_atom, int target_atom, clipper::Coord_orth &current_atom_orth,
                   clipper::Coord_orth &target_atom_orth);

        void find_nearby_nodes(const clipper::MAtomNonBond& nb, int a);

        void traverse_chain(std::shared_ptr<Node>& node, std::vector<int> &chain);

        static std::vector<std::vector<int>>  find_unique_chains(const std::vector<std::vector<int>>& chains);

        void find_cyclic_chains(std::vector<std::vector<int>> &chains, std::set<int> &visited);

        void symmetrise_chains(std::vector<std::vector<int>> &chains);


        // Model building main functions
        void identify_and_resolve_branches();

        void identify_and_resolve_clashes(std::vector<clipper::MMonomer> &monomers,
                                          std::vector<int> &chain, std::vector<double> &scores);

        clipper::MiniMol build_chains();

        FragmentResult build_chain(std::vector<int> &chain);


        // Model building utility functions
        FragmentResult fit_best_fragment(int n1, int n2, int n3);

        double fit_and_score_fragment(int n1, int n2, int n3);

        std::vector<double> fit_and_score_fragment_individually(int n1, int n2, int n3);

        double extract_library_fragment_and_score(int l, std::vector<clipper::Coord_orth> &reference_coords);

        FragmentResult extract_library_fragment_score_and_return(int l, std::vector<clipper::Coord_orth> &reference_coords);

        double score_monomers(std::vector<clipper::MMonomer>& monomers);

        std::vector<double> score_monomers_individually(std::vector<clipper::MMonomer>& monomers);

        double score_monomer(clipper::MMonomer& monomer, bool use_predicted_maps, bool use_experimental_map);

        static double score_to_grid(const clipper::Coord_orth& coord, const clipper::Xmap<float>* grid);

        clipper::Coord_orth get_symmetry_copy(clipper::Coord_orth& target, clipper::Coord_orth& reference);

        std::vector<std::vector<int>> find_local_chains(std::vector<int> &chain);

        // Graph Structure Utility Functions
        std::shared_ptr<Node> find_node_by_point_index(int point_index) const {
            for (const auto& node : nodes) {
                if (node->point_index == point_index) {
                    return node;
                }
            }
            return nullptr;
        }

        std::vector<std::shared_ptr<Node>> find_nodes_by_degree(int degree) const {
            std::vector<std::shared_ptr<Node>> results;
            for (const auto& node : nodes) {
                if (node->degree() == degree) {
                    results.push_back(node);
                }
            }
            return results;
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
        clipper::MAtomNonBond initial_nb; // Used to determine overlaps with input

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
