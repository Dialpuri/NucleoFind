//
// Created by Jordan Dialpuri on 21/04/2025.
//

#include "backbone-tracing.h"

#include "src/cpp/nautilus-util.h"

void NucleoFind::BackboneTracer::determine_edge(int source_atom, int target_atom,
                                                clipper::Coord_orth &current_atom_orth,
                                                clipper::Coord_orth &target_atom_orth) {
    float distance = clipper::Coord_orth::length(current_atom_orth, target_atom_orth);
    if (distance < 4 || distance > 8) return;

    clipper::Vec3<> midpoint = (current_atom_orth + target_atom_orth) * 0.5;
    auto midpoint_orth = clipper::Coord_orth(midpoint);
    auto sugar_predictions_nearby = predicted_sugar_positions.atoms_near(midpoint_orth, 1);
    if (sugar_predictions_nearby.empty()) return;

    auto edge = std::make_shared<Edge>(source_atom, target_atom, distance);
    edges.push_back(edge);

    nodes[source_atom]->add_edge(edge);
    nodes[target_atom]->add_edge(edge);
}

void NucleoFind::BackboneTracer::find_nearby_nodes(const clipper::MAtomNonBond &nb, int a) {
    clipper::Coord_orth current_atom_coord_orth = input[a].coord_orth();
    clipper::Coord_frac current_atom_coord_frac = current_atom_coord_orth.coord_frac(xgrid.cell());
    std::vector<clipper::MAtomIndexSymmetry> nearby = nb(current_atom_coord_orth, max_distance);

    for (const auto &near: nearby) {
        if (near.atom() == a) continue;
        clipper::Coord_orth nearby_atom_coord_orth = input[near.atom()].coord_orth();
        clipper::Coord_frac nearby_atom_coord_frac = nearby_atom_coord_orth.coord_frac(xgrid.cell());
        clipper::Coord_frac nearby_atom_symmetry_copy = nearby_atom_coord_frac.symmetry_copy_near(
            xgrid.spacegroup(), xgrid.cell(), current_atom_coord_frac);
        clipper::Coord_orth nearby_atom_symmetry_copy_orth = nearby_atom_symmetry_copy.coord_orth(xgrid.cell());

        determine_edge(a, near.atom(), current_atom_coord_orth, nearby_atom_symmetry_copy_orth);
    }
}

void NucleoFind::BackboneTracer::generate_graph() {
    nodes.reserve(input.size());
    const auto nb = clipper::MAtomNonBond(mol, max_distance);

    for (int a = 0; a < input.size(); a++) {
        nodes.push_back(std::make_shared<Node>(a));
    }

    for (int a = 0; a < input.size(); a++) {
        find_nearby_nodes(nb, a);
    }
}


double NucleoFind::BackboneTracer::score_to_grid(const clipper::Coord_orth &coord, const clipper::Xmap<float> *grid) {
    return grid->interp<clipper::Interp_linear>(coord.coord_frac(grid->cell()));
}

double NucleoFind::BackboneTracer::score_monomer(clipper::MMonomer &monomer) {
    double score = 0.0;

    int ip = monomer.lookup(" P  ", clipper::MM::ANY);
    int io5 = monomer.lookup(" O5'", clipper::MM::ANY);
    int ic5 = monomer.lookup(" C5'", clipper::MM::ANY);
    int ic4 = monomer.lookup(" C4'", clipper::MM::ANY);
    int io4 = monomer.lookup(" O4'", clipper::MM::ANY);
    int ic3 = monomer.lookup(" C3'", clipper::MM::ANY);
    int io3 = monomer.lookup(" O3'", clipper::MM::ANY);
    int ic2 = monomer.lookup(" C2'", clipper::MM::ANY);
    int ic1 = monomer.lookup(" C1'", clipper::MM::ANY);
    int in = monomer.lookup(" N9 ", clipper::MM::ANY);
    if (in < 0) in = monomer.lookup(" N1 ", clipper::MM::ANY);

    if (ip != -1) score += score_to_grid(monomer[ip].coord_orth(), predicted_maps.get_phosphate_map());
    if (io5 != -1) score += score_to_grid(monomer[io5].coord_orth(), predicted_maps.get_sugar_map());
    if (ic5 != -1) score += score_to_grid(monomer[ic5].coord_orth(), predicted_maps.get_sugar_map());
    if (ic4 != -1) score += score_to_grid(monomer[ic4].coord_orth(), predicted_maps.get_sugar_map());
    if (io4 != -1) score += score_to_grid(monomer[io4].coord_orth(), predicted_maps.get_sugar_map());
    if (ic3 != -1) score += score_to_grid(monomer[ic3].coord_orth(), predicted_maps.get_sugar_map());
    if (io3 != -1) score += score_to_grid(monomer[io3].coord_orth(), predicted_maps.get_sugar_map());
    if (ic2 != -1) score += score_to_grid(monomer[ic2].coord_orth(), predicted_maps.get_sugar_map());
    if (ic1 != -1) score += score_to_grid(monomer[ic1].coord_orth(), predicted_maps.get_sugar_map());
    if (in != -1) score += score_to_grid(monomer[in].coord_orth(), predicted_maps.get_sugar_map());

    return score;
}


double NucleoFind::BackboneTracer::score_monomers(std::vector<clipper::MMonomer> &monomers) {
    double score = 0;
    for (auto &monomer: monomers) {
        score += score_monomer(monomer);
    }
    return score;
}

std::vector<double> NucleoFind::BackboneTracer::score_monomers_individually(std::vector<clipper::MMonomer> &monomers) {
    std::vector<double> scores;
    for (auto &monomer: monomers) {
        scores.push_back(score_monomer(monomer));
    }
    return scores;
}

double NucleoFind::BackboneTracer::extract_library_fragment_and_score(
    int l, std::vector<clipper::Coord_orth> &reference_coords) {
    TriNucleotide library_fragment = library[l];
    std::vector<clipper::Coord_orth> library_fragment_positions = library_fragment.get_phosphates();
    clipper::RTop_orth rtop = {library_fragment_positions, reference_coords};
    std::vector<clipper::MMonomer> library_fragment_monomers = library_fragment.transform(rtop);
    return score_monomers(library_fragment_monomers);
}

double NucleoFind::BackboneTracer::fit_and_score_fragment(int n1, int n2, int n3) {
    double score = INT_MIN;
    const clipper::Coord_orth p1 = input[n1].coord_orth();
    const clipper::Coord_orth p2 = input[n2].coord_orth();
    const clipper::Coord_orth p3 = input[n3].coord_orth();

    std::vector<clipper::Coord_orth> trial_fragment_positions_forward = {p1, p2, p3};
    std::vector<clipper::Coord_orth> trial_fragment_positions_backward = {p3, p2, p1};

    for (int l = 0; l < library.size(); l++) {
        double forward_score = extract_library_fragment_and_score(l, trial_fragment_positions_forward);
        double backward_score = extract_library_fragment_and_score(l, trial_fragment_positions_backward);
        double best_score = std::max(forward_score, backward_score);
        score = std::max(score, best_score);
    }
    return score;
}

std::vector<double> NucleoFind::BackboneTracer::fit_and_score_fragment_individually(int n1, int n2, int n3) {
    double score = INT_MIN;
    int best_l = -1;
    bool fwd = false;

    const clipper::Coord_orth p1 = input[n1].coord_orth();
    const clipper::Coord_orth p2 = input[n2].coord_orth();
    const clipper::Coord_orth p3 = input[n3].coord_orth();

    std::vector<clipper::Coord_orth> trial_fragment_positions_forward = {p1, p2, p3};
    std::vector<clipper::Coord_orth> trial_fragment_positions_backward = {p3, p2, p1};

    for (int l = 0; l < library.size(); l++) {
        double forward_score = extract_library_fragment_and_score(l, trial_fragment_positions_forward);
        double backward_score = extract_library_fragment_and_score(l, trial_fragment_positions_backward);
        double best_score = std::max(forward_score, backward_score);

        if (best_score > score) {
            if (forward_score > backward_score) fwd = true;
            else fwd = false;

            score = best_score;
            best_l = l;
        }
    }

    TriNucleotide library_fragment = library[best_l];
    std::vector<clipper::Coord_orth> library_fragment_positions = library_fragment.get_phosphates();
    clipper::RTop_orth rtop;
    if (fwd)
        rtop = {library_fragment_positions, trial_fragment_positions_forward};
    else
        rtop = {library_fragment_positions, trial_fragment_positions_backward};

    std::vector<clipper::MMonomer> library_fragment_monomers = library_fragment.transform(rtop);
    std::string path = "superimposed-fragment-" + std::to_string(n1) + "-" + std::to_string(n2) + "-" + std::to_string(n3) + ".pdb";
    NautilusUtil::save_minimol(create_clipper_minimol(library_fragment_monomers, xgrid), path);
    std::vector<double> scores = score_monomers_individually(library_fragment_monomers);
    return scores;
}


void NucleoFind::BackboneTracer::identify_and_resolve_branches() {

    // identify which nodes are the midpoints of branches
    std::unordered_set<int> branch_nodes = {};
    for (const auto &node: nodes) {
        int node1 = node->point_index;
        auto it1 = adjacency.find(node1);
        if (it1 == adjacency.end()) continue;


        for (const auto &edge1: it1->second) {
            auto node2 = edge1->target;
            if (node1 == node2) continue;

            auto it2 = adjacency.find(node2);
            if (it2 == adjacency.end()) continue;

            bool branching = it2->second.size() > 2;
            if (!branching) continue;
            branch_nodes.insert(node2);
        }
    }

    // go through all surrounding points of all branch nodes, and find the wrong edge(s).
    for (auto &branch_node: branch_nodes) {
        auto node = find_node_by_point_index(branch_node);

        std::unordered_set<int> nearby_nodes_s = {};
        for (auto& edge: node->edges) {
            nearby_nodes_s.insert(edge->source);
            nearby_nodes_s.insert(edge->target);
        }

        std::vector<int> nearby_nodes = {nearby_nodes_s.begin(), nearby_nodes_s.end()};

        std::unordered_set<int> best_triplet = {};
        double best_score = INT_MIN;
        for (int i = 0; i < nearby_nodes.size(); i++) {
            for (int j = i + 1; j < nearby_nodes.size(); j++) {
                if (nearby_nodes[i] == branch_node || nearby_nodes[j] == branch_node) {continue;}
                double score = fit_and_score_fragment(nearby_nodes[i], branch_node, nearby_nodes[j]);
                if (score > best_score) {
                    best_score = score;
                    best_triplet = {nearby_nodes[i], branch_node, nearby_nodes[j]};
                }
            }
        }
        std::unordered_set<int> difference;
        for (const auto& elem : nearby_nodes_s) {
            if (best_triplet.find(elem) == best_triplet.end()) {
                difference.insert(elem);
            }
        }

        auto new_end = std::remove_if(edges.begin(), edges.end(), [&](const std::shared_ptr<Edge>& edge) {
            bool contains_wrong_node = difference.find(edge->target) != difference.end() || difference.find(edge->source) != difference.end();
            bool contains_midpoint = edge->target == branch_node || edge->source == branch_node;
            return contains_wrong_node && contains_midpoint;
       });
       edges.erase(new_end, edges.end());

    }


    // std::cout << "Branch point detected after " << node1 << " and " << node2 << std::endl;
    // std::cout << "choices are: " << std::endl;
    //
    // double best_score = INT_MIN;
    // int best_node = -1;
    // for (auto& edge2 : it2->second) {
    //     auto node3 = edge2->target;
    //     if (node3 == node1) continue;
    //     double score = fit_and_score_fragment(node1, node2, node3);
    //     std::cout << "\t" << node3 << " with score=" << score << std::endl;
    //     if (score > best_score) {
    //         best_score = score;
    //         best_node = node3;
    //     }
    // }
    //
    // std::cout << "Decision:" << node1 << " " << node2 << " -> " << best_node << " with score = " << best_score << std::endl;
    //
    // std::unordered_set<int> incorrect_branched_nodes = {};
    // for (auto& edge2 : it2->second) {
    //     auto node3 = edge2->target;
    //     if (node3 == node1 || node3 == best_node) continue;
    //     incorrect_branched_nodes.insert(node3);
    //     std::shared_ptr<Node> current_node = find_node_by_point_index(node3);
    //     current_node->remove_edge(edge2);
    //
    //     std::cout << "to remove: " << node3 << std::endl;
    // }
    //
    //
    // auto new_end = std::remove_if(edges.begin(), edges.end(), [&](const std::shared_ptr<Edge>& edge) {
    //     // when source == node2 and target == node3
    //     // when source == node3 and target == node2
    //         bool source_is_wrong = incorrect_branched_nodes.find(edge->source) != incorrect_branched_nodes.end();
    //         bool target_is_wrong = incorrect_branched_nodes.find(edge->target) != incorrect_branched_nodes.end();
    //         if (source_is_wrong || target_is_wrong) {
    //             std::cout << "Removing " << edge->source << " - " << edge->target << std::endl;
    //         }
    //         return source_is_wrong || target_is_wrong;
    // });
    // edges.erase(new_end, edges.end());

    // }
}


void NucleoFind::BackboneTracer::create_linked_structure() {
    clipper::MiniMol mol = {xgrid.spacegroup(), xgrid.cell()};
    clipper::MPolymer mp;
    mp.set_id("D");
    for (auto i = 0; i < edges.size(); i++) {
        // std::cout << edges[i]->source << "," << edges[i]->target << " " << edges[i]->distance << std::endl;

        auto source_orth = input[edges[i]->source].coord_orth();
        auto target_orth = input[edges[i]->target].coord_orth();

        clipper::Coord_frac source_frac = source_orth.coord_frac(xgrid.cell());
        clipper::Coord_frac target_frac = target_orth.coord_frac(xgrid.cell());
        clipper::Coord_frac nearby_target = target_frac.symmetry_copy_near(
            xgrid.spacegroup(), xgrid.cell(), source_frac);
        clipper::Coord_orth nearby_target_orth = nearby_target.coord_orth(xgrid.cell());
        clipper::Vec3<> p1p2 = (source_orth - nearby_target_orth) * 0.1;
        clipper::MMonomer monomer;
        monomer.set_type("P");
        monomer.set_id(i);

        for (int j = 1; j <= 10; j++) {
            clipper::Vec3<> p1p2_stepped = nearby_target_orth + (double(j) * p1p2);
            clipper::Coord_orth p1p2_stepped_orth = clipper::Coord_orth(p1p2_stepped);

            clipper::MAtom atom_p1 = NautilusUtil::create_atom(p1p2_stepped_orth, i + j, "O");
            monomer.insert(atom_p1);
        }
        mp.insert(monomer);
    }
    mol.model().insert(mp);
    NautilusUtil::save_minimol(mol, "triplet-drawn.pdb");

    std::ofstream output("graph.txt");
    for (const auto &edge: edges) {
        output << edge->source << "," << edge->target << "," << edge->distance << std::endl;
    }
    output.close();
}

void NucleoFind::BackboneTracer::print_stats() const {
    std::cout << "Graph Statistics:" << std::endl;
    std::cout << "  Nodes: " << nodes.size() << std::endl;
    std::cout << "  Edges: " << edges.size() << std::endl;

    // Degree distribution
    std::unordered_map<int, int> degreeCount;
    for (const auto &node: nodes) {
        degreeCount[node->degree()]++;
    }

    std::cout << "  Degree distribution:" << std::endl;
    for (const auto &pair: degreeCount) {
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

