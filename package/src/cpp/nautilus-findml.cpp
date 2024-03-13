//
// Created by Jordan Dialpuri on 06/03/2024.
//

#include "nautilus-findml.h"
#include "nautilus-refine.h"
#include "nautilus-tools.h"

FindML::FindML(const clipper::MiniMol& mol,
                    const clipper::Xmap<float>& xphospred,
                    const clipper::Xmap<float>& xsugarpred,
                    const clipper::Xmap<float>& xbasepred,
                    const clipper::Xmap<float>& xwrk) {
    this->mol = mol;
    this->xphospred = xphospred;
    this->xsugarpred = xsugarpred;
    this->xbasepred = xbasepred;
    this->xwrk = xwrk;
}

/*
 *
 * PREDICTIONS -> POINTS START
 *
 */

clipper::MiniMol FindML::generate_phosphate_molecule_from_gridpoints(double value_threshold) {
    clipper::MiniMol minimol(xwrk.spacegroup(), xwrk.cell());
    clipper::MModel m_model;
    clipper::MPolymer m_poly;
    clipper::MMonomer m_mon;
    m_mon.set_seqnum(1);
    m_mon.set_type("X");

    clipper::Grid_sampling grid_sampling = xphospred.grid_sampling();
    clipper::Coord_grid g0 = clipper::Coord_grid(0, 0, 0);
    clipper::Coord_grid g1 = clipper::Coord_grid(grid_sampling.nu(),
                                                 grid_sampling.nv(),
                                                 grid_sampling.nw());

    clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;
    i0 = clipper::Xmap_base::Map_reference_coord(xphospred, g0);

    int i = 0;
    for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u()) {
        for (iv = iu; iv.coord().v() <= g1.v(); iv.next_v()) {
            for (iw = iv; iw.coord().w() <= g1.w(); iw.next_w()) {
                if (xphospred[iw] >= value_threshold) {
                    clipper::MAtom m_atom = clipper::MAtom();
                    m_atom.set_id(std::to_string(i));
                    m_atom.set_name("X");
                    m_atom.set_coord_orth(iw.coord_orth());
                    m_mon.insert(m_atom);
                    i++;
                }
            }
        }
    }

    m_poly.insert(m_mon);
    m_model.insert(m_poly);
    minimol.model() = m_model;
    return minimol;
}

clipper::Coord_grid FindML::ascend_grid_gradient(const clipper::Coord_grid& grid_point,
                                                 const clipper::Xmap<float> &xmap) {
    clipper::Coord_grid best_gridpoint = grid_point;

    float max_value = xmap.get_data(grid_point);

    for (int i = -1; i <= 1; i += 1) {
        for (int j = -1; j <= 1; j += 1) {
            for (int k = -1; k <= 1; k += 1) {
                clipper::Coord_grid new_point = grid_point;
                new_point.u() += i;
                new_point.v() += j;
                new_point.w() += k;

                float value = xmap.get_data(new_point);

                if (value > max_value) {
                    max_value = value;
                    best_gridpoint = new_point;
                }
            }
        }
    }

    return best_gridpoint;


}

clipper::MiniMol FindML::find_phosphate_peaks(const clipper::MiniMol& phosphate_mol) const {
    clipper::MiniMol ascended_mol = clipper::MiniMol(xphospred.spacegroup(), xphospred.cell());

    clipper::MMonomer ascended_monomer;
    ascended_monomer.set_type("PA");
    ascended_monomer.set_id(1);

    std::vector <clipper::Coord_grid> atom_peaks;

    for (int poly = 0; poly < phosphate_mol.model().size(); poly++) {
        for (int mon = 0; mon < phosphate_mol.model()[poly].size(); mon++) {
            for (int atom = 0; atom < phosphate_mol.model()[poly][mon].size(); atom++) {
                constexpr int max_iter = 1000;

                clipper::Coord_grid m_atom_grid = phosphate_mol.model()[poly][mon][atom].coord_orth().coord_frac(
                        xphospred.cell()).coord_grid(xphospred.grid_sampling());

                clipper::Coord_grid current_best = m_atom_grid;

                for (int i = 0; i < max_iter; i++) {
                    clipper::Coord_grid best_neighbour = ascend_grid_gradient(current_best, xphospred);
                    if (best_neighbour.u() == current_best.u() &&
                        best_neighbour.v() == current_best.v() &&
                        best_neighbour.w() == current_best.w()) {
                        break;
                    }

                    current_best = best_neighbour;

                    if (i == max_iter - 1) {
                        std::cout << "Max gradient ascent reached, this indicates issues with the predicted map" << std::endl;
                    }
                }

                if (m_atom_grid != current_best) {
                    atom_peaks.emplace_back(current_best);
                }
            }
        }

        clipper::MPolymer mp;
        mp.set_id(1);
        mp.insert(ascended_monomer);
        ascended_mol.model().insert(mp);

        clipper::MiniMol peak_mol = clipper::MiniMol(xphospred.spacegroup(), xphospred.cell());

        clipper::MMonomer peak_monomer;
        peak_monomer.set_type("X");
        peak_monomer.set_id(0);

        for (int i = 0; i < atom_peaks.size(); i++) {
            peak_monomer.insert(
                    NautilusUtil::create_atom(atom_peaks[i].coord_frac(xphospred.grid_sampling()).coord_orth(xphospred.cell()), i,
                                              "X"));
        }

        clipper::MPolymer peak_mp;
        peak_mp.set_id(1);
        peak_mp.insert(peak_monomer);
        peak_mol.insert(peak_mp);

        return peak_mol;
    }

    return {};
}

clipper::MiniMol FindML::assimilate_phosphate_peaks(clipper::MiniMol& phosphate_peaks, float radius, const std::string& name) const {
    clipper::MModel m_model = phosphate_peaks.model();
    std::vector<clipper::MAtom> results;
    std::set<int> checked_atoms;

    clipper::MAtomNonBond non_bond = clipper::MAtomNonBond(phosphate_peaks, radius);

    for (int poly = 0; poly < m_model.size(); poly++) {
        for (int mon = 0; mon < m_model[poly].size(); mon++) {
            for (int atom = 0; atom < m_model[poly][mon].size(); atom++) {

                std::vector<clipper::MAtomIndexSymmetry> atom_list = non_bond.atoms_near(
                        m_model[poly][mon][atom].coord_orth(), radius);
                if (atom_list.empty()) {
                    continue;
                }

                if (checked_atoms.find(atom) != checked_atoms.end()) {
                    continue;
                }

                clipper::Coord_orth centroid_sum = {0, 0, 0};

                int count = 0;
                for (clipper::MAtomIndexSymmetry &atom_index: atom_list) {

                    if (checked_atoms.find(atom_index.atom()) ==
                        checked_atoms.end()) {
                        checked_atoms.insert(atom_index.atom());
                    }

                    clipper::Coord_orth diff =
                            m_model[atom_index.polymer()][atom_index.monomer()][atom_index.atom()].coord_orth() -
                            m_model[poly][mon][atom].coord_orth();

                    if (sqrt(diff.lengthsq()) == 0 || sqrt(diff.lengthsq()) > radius) {
                        continue;
                    }

                    centroid_sum += m_model[atom_index.polymer()][atom_index.monomer()][atom_index.atom()].coord_orth();
                    count++;
                }

                if (count > 0) {
                    clipper::Coord_orth centroid = {centroid_sum.x() / count, centroid_sum.y() / count,
                                                    centroid_sum.z() / count};

                    m_model[poly][mon][atom].set_coord_orth(centroid);
                }
                results.push_back(m_model[poly][mon][atom]);
            }
        }
    }


    clipper::MiniMol filtered_minimol(xwrk.spacegroup(), xwrk.cell());
    clipper::MModel filtered_m_model;
    clipper::MPolymer filtered_m_poly;
    clipper::MMonomer filtered_m_mon;
    filtered_m_mon.set_seqnum(1);
    filtered_m_mon.set_type(name);
    filtered_m_mon.set_id(name);

    int i = 0;
    for (auto &result: results) {
        result.set_id(i);
        filtered_m_mon.insert(result);
        i++;
    }
    filtered_m_poly.insert(filtered_m_mon);

    filtered_m_model.insert(filtered_m_poly);
    filtered_minimol.model() = filtered_m_model;
    return filtered_minimol;
}

clipper::MiniMol FindML::refine_phosphate_peaks(const clipper::MiniMol& phosphate_peaks) const {
    clipper::MiniMol mol_ = phosphate_peaks;

    for (int poly = 0; poly < mol_.model().size(); poly++) {
        for (int mon = 0; mon < mol_.model()[poly].size(); mon++) {
            for (int atom = 0; atom < mol_.model()[poly][mon].size(); atom++) {
                Target_fn_refine_phosphate refinement_target(xwrk, 0.1);
                clipper::Coord_orth refined_atom = refinement_target.refine(mol_.model()[poly][mon][atom].coord_orth());
                mol_.model()[poly][mon][atom].set_coord_orth(refined_atom);
            }
        }
    }

    return mol_;
}

FindML::TripletCoordinates FindML::symmetrise_phosphate_peaks(TripletCoordinates& triplet_coordinates, clipper::MiniMol& phosphate_peaks) const {

    std::set<int> visited = {};
    clipper::Coord_frac reference_point = triplet_coordinates[0][0].second.coord_frac(phosphate_peaks.cell());

    for (int triplet_idx = 0; triplet_idx < triplet_coordinates.size(); triplet_idx++) {
        TripletCoordinate triplet_coordinate = triplet_coordinates[triplet_idx];
        for (int fragment_idx = 0; fragment_idx < triplet_coordinate.size(); fragment_idx++) {
            auto data = triplet_coordinates[triplet_idx][fragment_idx];
            int current_id = data.first;

            if (visited.find(current_id) != visited.end()) {continue;}

            clipper::Coord_orth current_point = data.second;
            clipper::Coord_frac current_point_frac = current_point.coord_frac(phosphate_peaks.cell());
            current_point_frac = current_point_frac.symmetry_copy_near(phosphate_peaks.spacegroup(), phosphate_peaks.cell(), reference_point);
            triplet_coordinates[triplet_idx][fragment_idx].second = current_point_frac.coord_orth(phosphate_peaks.cell());
            // reference_point = current_point_frac;
            visited.insert(current_id);
        }
    }

    return triplet_coordinates;



    // for (int p = 0; p < phosphate_peaks.size(); p++) {
    //     for (int m = 0; m < phosphate_peaks[p].size(); m++) {
    //         for (int a = 0; a < phosphate_peaks[p][m].size(); a++) {
    //             clipper::Coord_orth point = phosphate_peaks[p][m][a].coord_orth();
    //             float distance = point.lengthsq();
    //             float best_distance = distance;
    //             clipper::Coord_orth best_coord = point;
    //
    //             for (int symop = 0; symop < phosphate_peaks.spacegroup().num_symops(); symop++) {
    //                 clipper::Coord_orth new_point = point;
    //                 clipper::RTop_orth new_point_rtop = phosphate_peaks.spacegroup().symop(symop).rtop_orth(phosphate_peaks.cell());
    //                 clipper::Coord_orth transformed_point = new_point.transform(new_point_rtop);
    //                 // if (transformed_point.x() < 0 || transformed_point.y() < 0 || transformed_point.z() < 0) {
    //                 //     continue;
    //                 // }
    //                 float new_distance = transformed_point.lengthsq();
    //                 if (new_distance < best_distance) {
    //                     best_coord = transformed_point;
    //                     best_distance = new_distance;
    //                 }
    //             }
    //
    //             if (best_distance != distance) {
    //                 phosphate_peaks[p][m][a].set_coord_orth(best_coord);
    //             }
    //         }
    //     }
    // }
    //
    // return phosphate_peaks;

    // for (int a = 0; a < phosphate_peaks[0][0].size(); a++) {
    //     clipper::Coord_orth point = phosphate_peaks[0][0][a].coord_orth();
    //     std::cout << point.format() << std::endl;
    //     float distance = point.lengthsq();
    //
    //     for (int symop = 0; symop < phosphate_peaks.spacegroup().num_symops(); symop++) {
    //         clipper::Coord_orth new_point = point;
    //         clipper::RTop_orth new_point_rtop = phosphate_peaks.spacegroup().symop(symop).rtop_orth(phosphate_peaks.cell());
    //         clipper::Coord_orth transformed_point = new_point.transform(new_point_rtop);
    //         float new_distance = transformed_point.lengthsq();
    //
    //         // std::cout << "Distance = " << distance << " " << new_distance << std::endl;
    //     }
    // }



}


clipper::MiniMol FindML::calculate_phosphate_peaks(double value_threshold) {
    clipper::MiniMol phosphate_mol = generate_phosphate_molecule_from_gridpoints(value_threshold);
    clipper::MiniMol phosphate_peak_groups = find_phosphate_peaks(phosphate_mol);
    clipper::MiniMol phosphate_peaks = assimilate_phosphate_peaks(phosphate_peak_groups, 1.5, "P");
    clipper::MiniMol refined_phosphate_peaks = refine_phosphate_peaks(phosphate_peaks);
    return refined_phosphate_peaks;
}

/*
 *
 * PREDICTIONS -> POINTS END
 *
 */

/*
 *
 * PREDICTED POINTS -> TRIPLETS START
 *
 */

FindML::TripletCoordinates FindML::find_triplet_coordinates(const clipper::MiniMol& phosphate_peaks) {
    float radius = 8;
    clipper::MAtomNonBond m_atom_non_bond = clipper::MAtomNonBond(phosphate_peaks, radius);
    clipper::MiniMol mol_ = phosphate_peaks;

    int p = 0;
    int m = 0;
    std::cout << mol_[p][m].size() << " phosphates found" << std::endl;

    TripletCoordinates return_list;

    float target_angle = 150;
    float target_range = 30;

    for (int atom = 0; atom < phosphate_peaks[p][m].size(); atom++) {

        clipper::MAtom m_atom = mol_.model()[p][m][atom];
        mol_.model()[p][m][atom].set_id(atom);
        clipper::Coord_orth m_atom_orth = m_atom.coord_orth();
        clipper::Coord_frac m_atom_frac = m_atom_orth.coord_frac(mol_.cell());
        std::vector<clipper::MAtomIndexSymmetry> atom_list = m_atom_non_bond(m_atom_orth, radius);

        for (auto &first_atom: atom_list) {

            clipper::Coord_frac first_atom_frac = mol_[p][m][first_atom.atom()].coord_orth().coord_frac(mol_.cell());
            first_atom_frac = first_atom_frac.symmetry_copy_near(mol_.spacegroup(), mol_.cell(), m_atom_frac);
            clipper::Coord_orth first_atom_orth = first_atom_frac.coord_orth(mol_.cell());

            for (auto &third_atom: atom_list) {

                if (first_atom.atom() == third_atom.atom()) continue;
                if (first_atom.atom() == atom) continue;
                if (third_atom.atom() == atom) continue;

                clipper::Coord_frac third_atom_frac = mol_[p][m][third_atom.atom()].coord_orth().coord_frac(mol_.cell());
                third_atom_frac = third_atom_frac.symmetry_copy_near(mol_.spacegroup(), mol_.cell(), m_atom_frac);
                clipper::Coord_orth third_atom_orth = third_atom_frac.coord_orth(mol_.cell());

                float first_second_distance = (m_atom_orth - first_atom_orth).lengthsq();
                float second_third_distance = (third_atom_orth - m_atom_orth).lengthsq();

                if (first_second_distance < 1 || second_third_distance < 1) continue;

                clipper::Vec3<> p1_v = first_atom_orth;
                clipper::Vec3<> p2_v = m_atom_orth;
                clipper::Vec3<> p3_v = third_atom_orth;

                clipper::Vec3<> AB = p1_v - p2_v;
                clipper::Vec3<> CB = p3_v - p2_v;

                clipper::Vec3<> AB_unit = AB.unit();
                clipper::Vec3<> CB_unit = CB.unit();

                double AB_x_BC = clipper::Vec3<>::dot(AB_unit, CB_unit);

                double angle = acos(AB_x_BC);
                double angle_d = clipper::Util::rad2d(angle);

                if (target_angle - target_range < angle_d && angle_d < target_angle + target_range) {
                    return_list.push_back({
                                                  {first_atom.atom(), first_atom_orth},
                                                  {atom,              m_atom_orth},
                                                  {third_atom.atom(), third_atom_orth}
                                          });
                }


            }
        }
    }

    return return_list;
}

/*
 *
 * PREDICTED POINTS -> TRIPLETS END
 *
 */

clipper::Coord_orth FindML::calculate_com(NucleicAcidDB::ChainFull &chain) {

    clipper::Coord_orth cm(0.0, 0.0, 0.0);
    double sm = 0.0;
    for (int i = 0; i < chain.size(); i++) {
        clipper::MMonomer monomer = chain[i].get_mmonomer();

        for (int atom = 0; atom < monomer.size(); atom++) {
            cm += monomer[atom].coord_orth();
            sm += 1.0;
        }
    }
    cm = (1.0 / sm) * cm;

    return cm;
}

float FindML::score_density(NucleicAcidDB::NucleicAcidFull &chain, clipper::Xmap<float> &xmap, bool terminal) {
    float score = 0.0f;

    if (terminal) {
        if (!chain.P.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.P.coord_frac(xmap.cell()));
        return score;
    }
    if (!chain.O5p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.O5p1.coord_frac(xmap.cell()));
    if (!chain.C5p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.C5p1.coord_frac(xmap.cell()));
    if (!chain.C4p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.C4p1.coord_frac(xmap.cell()));
    if (!chain.O4p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.O4p1.coord_frac(xmap.cell()));
    if (!chain.C3p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.C3p1.coord_frac(xmap.cell()));
    if (!chain.O3p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.O3p1.coord_frac(xmap.cell()));
    if (!chain.C2p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.C2p1.coord_frac(xmap.cell()));
    if (!chain.C1p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.C1p1.coord_frac(xmap.cell()));
    if (!chain.N1_1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.N1_1.coord_frac(xmap.cell()));
    if (!chain.N9_1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.N9_1.coord_frac(xmap.cell()));

    return score;
}

float FindML::score_sugar(NucleicAcidDB::NucleicAcidFull &chain) {
    float score = 0.0f;
    if (!chain.O5p1.is_null()) score += xsugarpred.interp<clipper::Interp_cubic>(chain.O5p1.coord_frac(xsugarpred.cell()));
    if (!chain.C5p1.is_null()) score += xsugarpred.interp<clipper::Interp_cubic>(chain.C5p1.coord_frac(xsugarpred.cell()));
    if (!chain.C4p1.is_null()) score += xsugarpred.interp<clipper::Interp_cubic>(chain.C4p1.coord_frac(xsugarpred.cell()));
    if (!chain.O4p1.is_null()) score += xsugarpred.interp<clipper::Interp_cubic>(chain.O4p1.coord_frac(xsugarpred.cell()));
    if (!chain.C3p1.is_null()) score += xsugarpred.interp<clipper::Interp_cubic>(chain.C3p1.coord_frac(xsugarpred.cell()));
    if (!chain.O3p1.is_null()) score += xsugarpred.interp<clipper::Interp_cubic>(chain.O3p1.coord_frac(xsugarpred.cell()));
    if (!chain.C2p1.is_null()) score += xsugarpred.interp<clipper::Interp_cubic>(chain.C2p1.coord_frac(xsugarpred.cell()));
    if (!chain.C1p1.is_null()) score += xsugarpred.interp<clipper::Interp_cubic>(chain.C1p1.coord_frac(xsugarpred.cell()));
    return score;
}

float FindML::score_base(NucleicAcidDB::NucleicAcidFull &chain) {
    float score = 0.0f;
    if (!chain.N1_1.is_null()) score += xbasepred.interp<clipper::Interp_cubic>(chain.N1_1.coord_frac(xbasepred.cell()));
    if (!chain.N2_1.is_null()) score += xbasepred.interp<clipper::Interp_cubic>(chain.N2_1.coord_frac(xbasepred.cell()));
    if (!chain.N3_1.is_null()) score += xbasepred.interp<clipper::Interp_cubic>(chain.N3_1.coord_frac(xbasepred.cell()));
    if (!chain.N4_1.is_null()) score += xbasepred.interp<clipper::Interp_cubic>(chain.N4_1.coord_frac(xbasepred.cell()));
    if (!chain.N6_1.is_null()) score += xbasepred.interp<clipper::Interp_cubic>(chain.N6_1.coord_frac(xbasepred.cell()));
    if (!chain.N7_1.is_null()) score += xbasepred.interp<clipper::Interp_cubic>(chain.N7_1.coord_frac(xbasepred.cell()));
    if (!chain.N9_1.is_null()) score += xbasepred.interp<clipper::Interp_cubic>(chain.N9_1.coord_frac(xbasepred.cell()));
    return score;
}


float FindML::score_fragment(NucleicAcidDB::ChainFull &fragment, clipper::Xmap<float> &xmap) {

    float total_score = 0.0f;
    for (int i = 0; i < fragment.size(); i++) {
        float score = score_density(fragment[i], xmap, i==fragment.size()-1);
        // score += score_sugar(fragment[i]);
        // score += score_base(fragment[i]);
        fragment[i].score = score;
        total_score += score;
    }
    fragment.chain_score = total_score;
    return total_score;
}

NucleicAcidDB::ChainFull
FindML::refine_fragment(NucleicAcidDB::ChainFull &original_fragment, float translation_range, float translation_step) {

    float best_score = -1e8;
    NucleicAcidDB::ChainFull best_chain;

    std::vector<std::vector<float>> translation_list;

    int steps = (2 * translation_range) / translation_step;

    for (int i = 0; i < steps; i++) {
        for (int j = 0; j < steps; j++) {
            for (int k = 0; k < steps; k++) {
                translation_list.push_back({
                                                   -translation_range + (i * translation_step),
                                                   -translation_range + (j * translation_step),
                                                   -translation_range + (k * translation_step),
                                           });
            }
        }
    }


    for (auto &translation: translation_list) {

        clipper::Coord_orth com = calculate_com(original_fragment);
        clipper::Coord_orth translated_com = {com.x() + translation[0], com.y() + translation[2],
                                              com.z() + translation[1]};
        clipper::RTop_orth com_rtop = {clipper::Mat33<>::identity(),
                                       -NautilusUtil::coord_to_vec(translated_com)};
        NucleicAcidDB::ChainFull com_fragment = original_fragment;
        com_fragment.transform(com_rtop);

        Target_fn_refine_fragment fragment_refiner = Target_fn_refine_fragment(xwrk, translated_com,
                                                                               com_fragment, 2);
        clipper::RTop_orth refined_rtop = fragment_refiner.refine();

        NucleicAcidDB::ChainFull test_fragment = com_fragment;
        test_fragment.transform(refined_rtop);

        float refined_score = score_fragment(test_fragment, xwrk);

        if (refined_score > best_score) {
            best_chain = test_fragment;
            best_score = refined_score;
        }
//                goto break_;
    }

    best_chain[0].set_triplet_id(0);
    best_chain[1].set_triplet_id(1);
    best_chain[2].set_triplet_id(2);

    return best_chain;
}

clipper::MiniMol FindML::remove_bases(clipper::MiniMol &mol) {

    std::vector<std::string> safe_list_pyr = {
        " C1'",
        " C2'",
        " C3'",
        " C4'",
        " C5'",
        " O3'",
        " O4'",
        " O5'",
        " P  ",
        " N1 ",
};

    std::vector<std::string> safe_list_pur = {
        " C1'",
        " C2'",
        " C3'",
        " C4'",
        " C5'",
        " O3'",
        " O4'",
        " O5'",
        " P  ",
        " N9 ",
};

    clipper::MiniMol safe_mol = {mol.spacegroup(), mol.cell()};


    for (int poly = 0; poly < mol.model().size(); poly++) {
        for (int mon = 0; mon < mol.model()[poly].size(); mon++) {
            clipper::MPolymer mp;
            clipper::MMonomer return_monomer;
            return_monomer.set_type(mol.model()[poly][mon].type());
            return_monomer.set_id(mol.model()[poly][mon].id());

            for (int atom = 0; atom < mol.model()[poly][mon].size(); atom++) {
                clipper::MAtom m_atom = mol.model()[poly][mon][atom];
                if (mol.model()[poly][mon].lookup(" N9 ", clipper::MM::ANY) >= 0) {
                    if (std::find(safe_list_pur.begin(), safe_list_pur.end(), m_atom.name()) != safe_list_pur.end()) {
                        return_monomer.insert(m_atom);
                    }
                } else if (mol.model()[poly][mon].lookup(" N1 ", clipper::MM::ANY) >= 0) {
                    if (std::find(safe_list_pyr.begin(), safe_list_pyr.end(), m_atom.name()) != safe_list_pyr.end()) {
                        return_monomer.insert(m_atom);
                    }
                }

            }

            mp.insert(return_monomer);
            safe_mol.model().insert(mp);

        }
    }

    return safe_mol;
}


clipper::MiniMol FindML::filter_and_form_chain(PossibleFragments &fragments) const {

    clipper::MiniMol mol = clipper::MiniMol(xwrk.spacegroup(), xwrk.cell());
    clipper::MPolymer mp;

    std::map<std::pair<int, int>, std::vector<NucleicAcidDB::NucleicAcidFull>> assimilated_fragments;

    int i = 0;
    for (auto & [fst, snd]: fragments) {
        if (assimilated_fragments.find(fst) != assimilated_fragments.end()) {
            assimilated_fragments[fst].insert(assimilated_fragments[fst].end(), snd.begin(), snd.end());
            continue;
        }

        auto fst_swap = std::make_pair(fst.second, fst.first);
        if (assimilated_fragments.find(fst_swap) != assimilated_fragments.end()) {
            assimilated_fragments[fst_swap].insert(assimilated_fragments[fst_swap].end(), snd.begin(), snd.end());
            continue;
        }

        if (assimilated_fragments.find(fst) == assimilated_fragments.end()) {
            assimilated_fragments[fst] = snd;
            continue;
        }
    }

    //
    // clipper::MiniMol mol = clipper::MiniMol(xwrk.spacegroup(), xwrk.cell());
    // clipper::MPolymer mp;
    //
    // int i = 0;
    // for (auto &frag: fragments) {
    //     for (auto &monomers: frag.second) {
    //         auto f1 = monomers.get_mmonomer();
    //         std::string id1 = "F" + std::to_string(frag.first.first) + std::to_string(frag.first.second);
    //         f1.set_type(id1);
    //         // f1.set_id(i + 1);
    //         f1.set_id(monomers.m_triplet_id);
    //         auto rounded_score = static_cast<float>(static_cast<int>(monomers.score * 10.0) /
    //                                                 10.0);
    //         for (int a = 0; a < f1.size(); a++) {
    //             f1[a].set_u_iso(rounded_score );
    //             f1[a].set_occupancy(1.0);
    //         }
    //
    //         i += 1;
    //         mp.insert(f1);
    //     }
    // }
    //
    // mol.model().insert(mp);
    // NautilusUtil::save_minimol(mol, "formed.pdb");
    //
    // return mol;

    // clipper::MiniMol mol = clipper::MiniMol(xwrk.spacegroup(), xwrk.cell());
    // clipper::MPolymer mp;



    for (auto &frag: assimilated_fragments) {
        std::sort(frag.second.begin(), frag.second.end());

        auto f1 = frag.second[frag.second.size() - 1].get_mmonomer();
        f1.set_type(frag.second[frag.second.size() - 1].get_type());
        // f1.set_id(i);
        std::string id1 = "F" + std::to_string(frag.first.first) + std::to_string(frag.first.second);

        f1.set_id(id1);
        auto rounded_score = static_cast<float>(static_cast<int>(frag.second[frag.second.size() - 1].score * 10.0) /
                                                10.0);
        for (int a = 0; a < f1.size(); a++) {
            f1[a].set_u_iso(rounded_score);
            f1[a].set_occupancy(1.0);
        }

        i += 1;
        mp.insert(f1);
    }

    mol.model().insert(mp);
    return mol;
}

void FindML::find_chains(int current_index, std::map<int, std::vector<int>> &connections, ChainData &chain) {
    chain.lookup_list.insert(current_index);
    chain.ordered_chain.emplace_back(current_index);

    for (auto const& nearby_P: connections.at(current_index)) {
        if (chain.lookup_list.find(nearby_P) != chain.lookup_list.end())
            continue;

        find_chains(nearby_P, connections, chain);
    }
}


clipper::MiniMol FindML::organise_to_chains(clipper::MiniMol &mol) {

    clipper::MAtomNonBond nb = clipper::MAtomNonBond(mol, 2);

    std::map<int, std::vector<int>> connections;
    std::map<int, std::vector<int>> connections_O3;

    for (int c = 0; c < mol.size(); c++) {
        connections[c] = {};

        int p = mol[c][0].lookup(" P  ", clipper::MM::UNIQUE);
        if (p == -1) continue;

        clipper::MAtom P = mol[c][0].find(" P  ", clipper::MM::UNIQUE);
        std::vector<clipper::MAtomIndexSymmetry> atoms_near = nb(P.coord_orth(), 2.5);

        for (auto const &atom: atoms_near) {
            if (mol[atom.polymer()][atom.monomer()][atom.atom()].name().trim() == "O3'")
                connections[c].emplace_back(atom.polymer());
        }

        connections_O3[c] = {};

        int o3 = mol[c][0].lookup(" O3'", clipper::MM::UNIQUE);
        if (o3 == -1) continue;

        clipper::MAtom O3 = mol[c][0].find(" O3'", clipper::MM::UNIQUE);
        std::vector<clipper::MAtomIndexSymmetry> atoms_near_o3 = nb(O3.coord_orth(), 2.5);

        for (auto const &atom: atoms_near_o3) {
            if (mol[atom.polymer()][atom.monomer()][atom.atom()].name().trim() == "P")
                connections_O3[c].emplace_back(atom.polymer());
        }
    }

    std::unordered_set<int> visited;
    std::vector<std::vector<int>> found_chains;
    clipper::MiniMol output_mol = {mol.spacegroup(), mol.cell()};

    int mon_id = 0;

    std::vector<int> seed_points = {};

    for (auto const &connection: connections) {
        if (connection.second.empty()) {
            seed_points.emplace_back(connection.first);
        }

    }

    for (auto& seed: seed_points) {
        ChainData chain;
        find_chains(seed, connections_O3, chain);
        clipper::MPolymer mp;

        for (auto const& id: chain.ordered_chain) {
            mp.insert(mol[id][0]);
        }
        output_mol.insert(mp);
    }
    NucleicAcidTools::chain_label(output_mol, clipper::MMDBManager::PDB);
//    NautilusUtil::save_minimol(output_mol, "debug/1u9s/chained_find.pdb");

    return output_mol;
}


clipper::MiniMol FindML::find() {
    clipper::MiniMol phosphate_peaks = calculate_phosphate_peaks(0.01);
    TripletCoordinates phosphate_triplets = find_triplet_coordinates(phosphate_peaks);
    // TripletCoordinates symmetrised_triplets = symmetrise_phosphate_peaks(phosphate_triplets, phosphate_peaks);
//    NautilusUtil::save_minimol(phosphate_peaks, "phosphate_peaks.pdb");
    std::map<std::pair<int, int>, std::vector<NucleicAcidDB::NucleicAcidFull>> placed_fragments;

    for (int i = 0; i < phosphate_triplets.size(); i++) {

        clipper::Coord_orth p1 = phosphate_triplets[i][0].second;
        clipper::Coord_orth p2 = phosphate_triplets[i][1].second;
        clipper::Coord_orth p3 = phosphate_triplets[i][2].second;

        std::vector<clipper::Coord_orth> triplet_pos = {p1, p2, p3};
        // if (phosphate_triplets[i][0].first != 6  && phosphate_triplets[i][0].first != 8) { continue;}

        float max_score = -1e8f;
        NucleicAcidDB::ChainFull best_fragment;

        for (int j = 0; j < nadb.size(); j++) {
            NucleicAcidDB::ChainFull fragment = nadb.extract(j, 3);

            if (!fragment.is_continuous()) continue;

            for (int x = 0; x < phosphate_triplets[i].size(); x++) {
                fragment[x].set_triplet_id(phosphate_triplets[i][x].first);
                // fragment[x].set_triplet_id(std::stoi(std::to_string(phosphate_triplets[i][0].first) + std::to_string(phosphate_triplets[i][1].first) +  std::to_string(phosphate_triplets[i][2].first)));
            }

            std::vector<clipper::Coord_orth> fragment_pos = {fragment[0].P, fragment[1].P, fragment[2].P};

            clipper::RTop_orth align = clipper::RTop_orth(fragment_pos, triplet_pos);

            fragment.transform(align);
            fragment.alignment = align;

            float total_score = score_fragment(fragment, xwrk);
            // std::cout << "Total Score = " << total_score << " - Max Score = " << max_score << std::endl;
            if (total_score > max_score) {
                max_score = total_score;
                best_fragment = fragment;
                // placed_fragments[phosphate_triplets[i][0].first].emplace_back(fragment[0]);
                // placed_fragments[phosphate_triplets[i][1].first].emplace_back(fragment[1]);
                // placed_fragments[phosphate_triplets[i][2].first].emplace_back(fragment[2]);
            }


        }

        NucleicAcidDB::ChainFull refined_fragment = refine_fragment(best_fragment, 1, 1);

        placed_fragments[std::make_pair(phosphate_triplets[i][0].first, phosphate_triplets[i][1].first)].emplace_back(best_fragment[0]);
        placed_fragments[std::make_pair(phosphate_triplets[i][1].first, phosphate_triplets[i][2].first)].emplace_back(best_fragment[1]);
        // placed_fragments[phosphate_triplets[i][2].first].emplace_back(best_fragment[2]);
        // break;

    }

    clipper::MiniMol filtered_chain = filter_and_form_chain(placed_fragments);
    clipper::MiniMol base_removed_mol = remove_bases(filtered_chain);

    clipper::MiniMol mol_ = organise_to_chains(base_removed_mol);
    for (int p = 0; p < mol_.size(); p++) {
        mol.insert(mol_[p]);
    }

    NautilusUtil::save_minimol(mol, "mlfind.pdb");
    return mol;
}
