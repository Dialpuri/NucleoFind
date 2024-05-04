//
// Created by Jordan Dialpuri on 06/03/2024.
//

#include "nautilus-findml.h"
#include "nautilus-refine.h"
#include "nautilus-tools.h"

FindML::FindML(const clipper::MiniMol &mol, const clipper::Xmap<float> &xwrk, PredictedMaps predictions)
        : predictions(predictions) {
    this->mol = mol;
    this->xwrk = xwrk;
    if (!predictions.get_phosphate_map().has_value()) {
        throw std::runtime_error("Cannot predict without an inputted phosphate predicted map");
    }
    xphospred = predictions.get_phosphate_map().value();
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
    if (!predictions.get_sugar_map().has_value()) return 0;
    clipper::Xmap<float> xsugarpred = predictions.get_sugar_map().value();
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
    if (!predictions.get_base_map().has_value()) return 0;
    clipper::Xmap<float> xbasepred = predictions.get_base_map().value();
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

NucleicAcidDB::ChainFull FindML::refine_fragment_coordinates(NucleicAcidDB::ChainFull &original_fragment) {


    clipper::Coord_orth com = calculate_com(original_fragment);
    clipper::RTop_orth com_rtop = {clipper::Mat33<>::identity(),-com};
    NucleicAcidDB::ChainFull com_fragment = original_fragment;
    com_fragment.transform(com_rtop);

    RefineFragmentCoordinates refiner = RefineFragmentCoordinates(xwrk, xphospred, com_fragment, com);

    clipper::RTop_orth refined_rtop = refiner.refine();

    NucleicAcidDB::ChainFull new_fragment = com_fragment;
    new_fragment.transform(refined_rtop);

//    for (int i = 0; i < new_fragment.size(); i++) {
//        std::cout << "P Atom moved  " << (original_fragment[i].P-new_fragment[i].P).lengthsq()<< std::endl;
//    }

    return new_fragment;
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
        clipper::MPolymer mp;
        mp.set_id(mol[poly].id());
        for (int mon = 0; mon < mol.model()[poly].size(); mon++) {
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
        }
        safe_mol.model().insert(mp);

    }

    return safe_mol;
}
clipper::MiniMol FindML::filter_and_form_bidirectional_chain(PossibleFragments& fragments) const {
    clipper::MiniMol mol_ = clipper::MiniMol(xwrk.spacegroup(), xwrk.cell());
    clipper::MPolymer mp;

    int i = 0;
    for (auto & [fst, snd]: fragments) {
        std::sort(snd.begin(), snd.end());

        auto f1 = snd[snd.size() - 1].get_mmonomer();
        f1.set_type(snd[snd.size() - 1].get_type());
        std::string id1 = "F" + std::to_string(fst.first) + std::to_string(fst.second);

        // f1.set_id(i);
        auto rounded_score = static_cast<float>(static_cast<int>(snd[snd.size() - 1].score * 10.0) /
                                                10.0);
        for (int a = 0; a < f1.size(); a++) {
            f1[a].set_u_iso(rounded_score);
            f1[a].set_occupancy(1.0);
        }

        i += 1;
        mp.insert(f1);
    }

    mol_.model().insert(mp);
    return mol_;
}

clipper::MiniMol FindML::form_organised_chains(PossibleFragments& fragments, std::vector<std::vector<int>>& fragment_indices) const {
    clipper::MiniMol mol_ = clipper::MiniMol(xwrk.spacegroup(), xwrk.cell());

    for (const auto& chain: fragment_indices) {
        clipper::MPolymer mp;
        mp.set_id("");
        for (int i = 0; i < chain.size()-1; i++) {
            std::pair<int, int> pair = std::make_pair(chain[i], chain[i+1]);

            if (fragments.find(pair) == fragments.end()) {
                std::cout << "not found a pair for " << pair.first << " " << pair.second << std::endl;
                continue;
            }

            std::vector<NucleicAcidDB::NucleicAcidFull> fragment_list = fragments[pair];
            std::sort(fragment_list.begin(), fragment_list.end());
            auto f1 = fragment_list[fragment_list.size() - 1].get_mmonomer();
            f1.set_type(fragment_list[fragment_list.size() - 1].get_type());
            f1.set_id(i);
            mp.insert(f1);
        }

        mol_.insert(mp);
    }

    NucleicAcidTools::chain_label(mol_, clipper::MMDBManager::PDB );
    return mol_;
}


clipper::MiniMol FindML::form_chain(PossibleFragments& fragments) const {
    clipper::MiniMol mol = clipper::MiniMol(xwrk.spacegroup(), xwrk.cell());
    clipper::MPolymer mp;

    int i = 0;
    for (auto &frag: fragments) {
        for (auto &monomers: frag.second) {
            auto f1 = monomers.get_mmonomer();
            std::string id1 = "F" + std::to_string(frag.first.first) + std::to_string(frag.first.second);
            f1.set_type(id1);
            // f1.set_id(i + 1);
            f1.set_id(monomers.m_triplet_id);
            auto rounded_score = static_cast<float>(static_cast<int>(monomers.score * 10.0) /
                                                    10.0);
            for (int a = 0; a < f1.size(); a++) {
                f1[a].set_u_iso(rounded_score );
                f1[a].set_occupancy(1.0);
            }

            i += 1;
            mp.insert(f1);
        }
    }
    mol.insert(mp);
    return mol;
}


clipper::MiniMol FindML::filter_and_form_single_chain(PossibleFragments&fragments) const {
    clipper::MiniMol mol = clipper::MiniMol(xwrk.spacegroup(), xwrk.cell());
    clipper::MPolymer ascending_mp;
    clipper::MPolymer descending_mp;

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
    clipper::MPolymer mp;

    for (auto &frag: assimilated_fragments) {
        std::sort(frag.second.begin(), frag.second.end());

        auto f1 = frag.second[frag.second.size() - 1].get_mmonomer();
        f1.set_type(frag.second[frag.second.size() - 1].get_type());
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

FindML::PairedChainIndices FindML::organise_triplets_to_chains(TripletCoordinates&triplets) {

    std::vector<std::vector<int>> chains;
    std::set<std::vector<int>> completed_triplets;

    /**Go through each triplet and look to see whether there is a chain where it can be added to the front or back, if there is
     * do so. If not, create a new chain for that triplet.
     * Should create a large chain of triplets, but the output can be split depending on the start points.
     ***/
    for (TripletCoordinate& triplet: triplets) {
        bool found_triplet = false;
        std::vector<int> triplet_ids = {triplet[0].first, triplet[1].first, triplet[2].first};
        if (completed_triplets.find(triplet_ids) != completed_triplets.end()) {continue;}

        for (auto & chain : chains) {
            int chain_length = chain.size();

            if (triplet_ids[0] == chain[chain_length-2] && triplet_ids[1] == chain[chain_length-1]) {
                chain.push_back(triplet_ids[2]);
                found_triplet = true;
            }
            if (triplet_ids[1] == chain[0] && triplet_ids[2] == chain[1]) {
                chain.insert(chain.begin(), triplet_ids[0]);
                found_triplet = true;
            }
        }
        if (!found_triplet) {
            chains.push_back({triplet[0].first, triplet[1].first, triplet[2].first});
        }

        completed_triplets.insert(triplet_ids);
    }

    /*chains could have become split depending on starting point, bring them together*/
    bool merging = true;
    while (merging) {
        merging = false;
        for (auto it1 = chains.begin(); it1 != chains.end(); ++it1) {
            for (auto it2 = it1+1; it2 != chains.end(); ++it2) {
                if (std::equal((*it1).end()-3, (*it1).end(), (*it2).begin())) {
                    (*it1).insert((*it1).end(), (*it2).begin() + 3, (*it2).end());
                    it2 = chains.erase(it2);
                    merging = true;
                    break;
                }
                if (std::equal((*it1).begin(), (*it1).begin()+3, (*it2).end()-3)) {
                    (*it2).insert((*it2).end(), (*it1).begin()+3, (*it1).end());
                    it1 = chains.erase(it1);
                    merging = true;
                    break;
                }
            }
            if (merging) {
                    break;
            }
        }
    }

    std::map<std::set<int>, int> paired_map = {};
    PairedChainIndices pairs;

    for (int chain_idx = 0; chain_idx < chains.size(); chain_idx++) {
        std::set<int> set = {chains[chain_idx].begin(), chains[chain_idx].end()};
        if (paired_map.find(set) == paired_map.end()) {
            paired_map.insert({set, chain_idx});
        } else {
            pairs.emplace_back(chains[paired_map[set]], chains[chain_idx]);
        }
    }

    return pairs;
}


clipper::MiniMol FindML::organise_to_chains(clipper::MiniMol &mol) {

    clipper::MAtomNonBond nb = clipper::MAtomNonBond(mol, 2);

    std::map<int, std::vector<int>> connections;
    std::map<int, std::vector<int>> connections_O3;

    for (int c = 0; c < mol.size(); c++) {
        connections[c] = {};

        int p = mol[c][0].lookup(" P  ", clipper::MM::UNIQUE);
        if (p == -1) {std::cout << "Couldn't find P" << std::endl; continue;}

        clipper::MAtom P = mol[c][0][p];
        std::vector<clipper::MAtomIndexSymmetry> atoms_near = nb(P.coord_orth(), 2.5);

        for (auto const &atom: atoms_near) {
            if (mol[atom.polymer()][atom.monomer()][atom.atom()].name().trim() == "O3'")
                connections[c].emplace_back(atom.polymer());
        }

        connections_O3[c] = {};

        int o3 = mol[c][0].lookup(" O3'", clipper::MM::UNIQUE);
        if (o3 == -1) continue;

        clipper::MAtom O3 = mol[c][0][o3];
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

    if (seed_points.empty()) {
        seed_points.emplace_back(connections.begin()->first);
    }

    for (auto& seed: seed_points) {
        ChainData chain;
        find_chains(seed, connections_O3, chain);
        clipper::MPolymer mp;

        int i = 0;
        for (auto const& id: chain.ordered_chain) {
            clipper::MMonomer mon = mol[id][0];
            mon.set_id(i++);
            mp.insert(mon);
        }
        output_mol.insert(mp);
    }
    NucleicAcidTools::chain_label(output_mol, clipper::MMDBManager::PDB);
    // NautilusUtil::save_minimol(output_mol, "organised.pdb");

    return output_mol;
}


clipper::MiniMol FindML::remove_clashing_protein(clipper::MiniMol& na_chain) {
    clipper::MiniMol molwrk = mol;
    clipper::MAtomNonBond nb = clipper::MAtomNonBond(molwrk, 4);

    std::set<std::vector<std::string>> to_remove;

    clipper::MiniMol mol_final = {mol.spacegroup(), mol.cell()};

    for (int p = 0; p < na_chain.size(); p++) {
        clipper::MPolymer mp;
       mp.set_id(na_chain[p].id());
        for (int m = 0; m < na_chain[p].size(); m++) {
            for (int a = 0; a < na_chain[p][m].size(); a++) {
                auto nearby_atoms = nb.atoms_near(na_chain[p][m][a].coord_orth(), 1);

                for (const auto &nearby_atom: nearby_atoms) {
                    clipper::MPolymer chain = molwrk[nearby_atom.polymer()];
                    clipper::MMonomer residue = chain[nearby_atom.monomer()];
                    NucleicAcidDB::NucleicAcid na = {residue};
                    if (na.flag() == NucleicAcidDB::NucleicAcid::NONE) {
                        to_remove.insert({
                                                 chain.id(), residue.type(), std::to_string(residue.seqnum())
                                         });
                    }
                }
            }

            mp.insert(na_chain[p][m]);
        }
        mol_final.insert(mp);
    }

    std::set<std::string> allowed_atoms = {"CA", "C", "N", "O"};

    for (int p = 0; p < mol.size(); p++) {
        clipper::MPolymer mp;
        mp.set_id(mol[p].id());
        int count = 0;
        for (int m = 0; m < mol[p].size(); m++) {
            clipper::MPolymer chain = mol[p];
            clipper::MMonomer residue = mol[p][m];
            std::vector<std::string> key = {chain.id(), residue.type(), std::to_string(residue.seqnum())};
            if (to_remove.find(key) != to_remove.end()) {
                // clipper::MMonomer backbone_only;
                // backbone_only.set_type(mol[p][m].type());
                // backbone_only.set_id(mol[p][m].id());

                // for (int a = 0; a < mol[p][m].size(); a++) {
                //     if (allowed_atoms.find(mol[p][m][a].id().trim()) != allowed_atoms.end()) {
                //         backbone_only.insert(mol[p][m][a]);
                //     }
                // }

                // if (backbone_only.size() != 0) {
                    // mp.insert(backbone_only);
                    // count += 1;
                // }

                continue;
            }
            count += 1;
            mp.insert(mol[p][m]);
        }
        if (count > 0) {
            mol_final.insert(mp);

        }
    }


    // NautilusUtil::save_minimol(mol_final, "clash_mol.pdb");
    // Need to flag chains so join knows which is protein and which is NA.
    mol_final = NucleicAcidTools::flag_chains(mol_final);

    return mol_final;
}

clipper::MiniMol FindML::remove_low_confidence(clipper::MiniMol & mol) {
    NucleicAcidTools::chain_label(mol, clipper::MMDBManager::PDB);

    // Set the IDs of the residues, for use as keys later
    NucleicAcidTools::residue_label(mol);

    std::map<std::pair<std::string, std::string>, double> rsrs = NautilusUtil::per_residue_rsrz(mol, xwrk, m_resolution);
    double rsr_threshold = -1.5;
    clipper::MiniMol mol_final = {mol.spacegroup(), mol.cell()};

    for (int p = 0; p < mol.size(); p++) {
        clipper::MPolymer mp;
        mp.set_id(mol[p].id());
        int count = 0;
        for (int m = 0; m < mol[p].size(); m++) {
            std::pair<std::string, std::string> residue_key = std::make_pair(mol[p].id().trim(), mol[p][m].id().trim());
            if (rsrs.find(residue_key) != rsrs.end()) {
                double rsrz = rsrs[residue_key];
                if (rsrz < rsr_threshold) {continue;}
            }
            else {
                std::cout << "Residue found that was not in the RSR calculation " << mol[p].id().trim() << " " << mol[p][m].id().trim() << std::endl;
                continue;
            }

            clipper::MPolymer chain = mol[p];
            clipper::MMonomer residue = mol[p][m];
            std::vector<std::string> key = {chain.id(), residue.type(), std::to_string(residue.seqnum())};
            count += 1;
            mp.insert(mol[p][m]);
        }
        if (count > 0)
            mol_final.insert(mp);
    }

    return mol_final;
}

PlacedFragmentResult FindML::place_fragments(const clipper::MiniMol& phosphate_peaks, const std::vector<int>& positions) {
    float score = 0;
    std::map<std::pair<int, int>, std::vector<NucleicAcidDB::NucleicAcidFull>> placed_fragments;

    for (int i = 0; i < positions.size()-2; i++) {
        clipper::Coord_orth p1 = phosphate_peaks[0][0][positions[i]].coord_orth();
        clipper::Coord_orth p2 = phosphate_peaks[0][0][positions[i+1]].coord_orth();
        clipper::Coord_orth p3 = phosphate_peaks[0][0][positions[i+2]].coord_orth();
        p1 = p1.coord_frac(mol.cell()).symmetry_copy_near(mol.spacegroup(), mol.cell(), p2.coord_frac(mol.cell())).coord_orth(mol.cell());
        p3 = p3.coord_frac(mol.cell()).symmetry_copy_near(mol.spacegroup(), mol.cell(), p2.coord_frac(mol.cell())).coord_orth(mol.cell());

        std::vector<clipper::Coord_orth> triplet_pos = {p1, p2, p3};

        float max_score = -1e8f;
        NucleicAcidDB::ChainFull best_fragment;

        for (int j = 0; j < nadb.size() - 2; j++) {
            NucleicAcidDB::ChainFull fragment = nadb.extract(j, 3);

            if (!fragment.is_continuous()) continue;

            std::vector<clipper::Coord_orth> fragment_pos = {fragment[0].P, fragment[1].P, fragment[2].P};
            clipper::RTop_orth align = clipper::RTop_orth(fragment_pos, triplet_pos);

            fragment.transform(align);
            fragment.alignment = align;

            float total_score = score_fragment(fragment, xwrk);
            if (total_score > max_score) {
                max_score = total_score;
                best_fragment = fragment;
            }
        }

//        NucleicAcidDB::ChainFull refined_fragment = refine_fragment(best_fragment, 1, 0.2);
        NucleicAcidDB::ChainFull refined_fragment = refine_fragment_coordinates(best_fragment);

        score += max_score;
        placed_fragments[std::make_pair(positions[i], positions[i+1])].emplace_back(refined_fragment[0]);
        placed_fragments[std::make_pair(positions[i+1], positions[i+2])].emplace_back(refined_fragment[1]);
    }

    return {score, placed_fragments};
}

clipper::MiniMol FindML::find() {
    clipper::MiniMol phosphate_peaks = calculate_phosphate_peaks(0.1);
    TripletCoordinates phosphate_triplets = find_triplet_coordinates(phosphate_peaks);
     draw_triplets(phosphate_triplets, phosphate_peaks, "triplets-ext.pdb");
//     NautilusUtil::save_minimol(phosphate_peaks, "phosphate_peaks.pdb");
    std::cout << phosphate_triplets.size() << " phosphate triplets found\n";

    PairedChainIndices pairs = organise_triplets_to_chains(phosphate_triplets);

//    for (const auto& pair : pairs) {
//        std::cout << "Pair:" << std::endl;
//        std::cout << "First vector: ";
//        for (const auto& element : pair.first) {
//            std::cout << element << " ";
//        }
//        std::cout << std::endl;
//
//        std::cout << "Second vector: ";
//        for (const auto& element : pair.second) {
//            std::cout << element << " ";
//        }
//        std::cout << std::endl << std::endl;
//    }

    std::map<std::pair<int, int>, std::vector<NucleicAcidDB::NucleicAcidFull>> placed_fragments;
    std::vector<std::vector<int>> placed_fragment_indices;
    for (const auto& [fwd, bck]: pairs) {
        PlacedFragmentResult forward_result = place_fragments(phosphate_peaks, fwd);
        PlacedFragmentResult backward_result = place_fragments(phosphate_peaks, bck);

        if (forward_result.score > backward_result.score) {
            placed_fragments.insert(forward_result.fragments.begin(), forward_result.fragments.end());
            placed_fragment_indices.emplace_back(fwd);
        } else {
            placed_fragments.insert(backward_result.fragments.begin(), backward_result.fragments.end());
            placed_fragment_indices.emplace_back(bck);
        }
    }

    clipper::MiniMol filtered_chain = form_organised_chains(placed_fragments, placed_fragment_indices);
    clipper::MiniMol base_removed_mol = remove_bases(filtered_chain);
    clipper::MiniMol low_confidence_removed_model = remove_low_confidence(base_removed_mol);
    clipper::MiniMol clash_removed_mol = remove_clashing_protein(low_confidence_removed_model);
    return clash_removed_mol;
}
