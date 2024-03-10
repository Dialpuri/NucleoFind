//
// Created by jordan on 06/06/23.
//

#include "nautilus-mlfind.h"

float PredictionBuilder::score_density(NucleicAcidDB::NucleicAcidFull &chain, clipper::Xmap<float> &xmap) {
    float score = 0.0f;

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


clipper::MiniMol PredictionBuilder::density_to_atom(clipper::Xmap<float> &xmap, float threshold_value) {
    clipper::MiniMol minimol(xmap.spacegroup(), xmap.cell());
    clipper::MModel m_model;
    clipper::MPolymer m_poly;
    clipper::MMonomer m_mon;
    m_mon.set_seqnum(1);
    m_mon.set_type("X");
    clipper::Coord_grid g0 = clipper::Coord_grid(0, 0, 0);
    clipper::Coord_grid g1 = clipper::Coord_grid(xmap.cell().a(),
                                                 xmap.cell().b(),
                                                 xmap.cell().c());

    int i = 0;
    clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;
    i0 = clipper::Xmap_base::Map_reference_coord(xmap, g0);
    for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u()) {
        for (iv = iu; iv.coord().v() <= g1.v(); iv.next_v()) {
            for (iw = iv; iw.coord().w() <= g1.w(); iw.next_w()) {
                if (xmap[iw] > threshold_value) {
                    clipper::MAtom m_atom = clipper::MAtom();
                    m_atom.set_id(std::to_string(i));
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


clipper::MiniMol PredictionBuilder::generate_phos_minimol(float threshold_value) {
    clipper::Xmap<float> &phos_map = m_prediction.phosphate_map;

    clipper::MiniMol minimol(p_xmap.spacegroup(), p_xmap.cell());
    clipper::MModel m_model;
    clipper::MPolymer m_poly;
    clipper::MMonomer m_mon;
    m_mon.set_seqnum(1);
    m_mon.set_type("X");

    clipper::Coord_grid g0 = clipper::Coord_grid(0, 0, 0);
    clipper::Coord_grid g1 = clipper::Coord_grid(phos_map.grid_sampling().nu(),
                                                 phos_map.grid_sampling().nv(),
                                                 phos_map.grid_sampling().nw());

    clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;
    i0 = clipper::Xmap_base::Map_reference_coord(phos_map, g0);

    int i = 0;
    for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u()) {
        for (iv = iu; iv.coord().v() <= g1.v(); iv.next_v()) {
            for (iw = iv; iw.coord().w() <= g1.w(); iw.next_w()) {
                if (phos_map[iw] >= threshold_value) {
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


clipper::MiniMol PredictionBuilder::generate_phosphate_positions() {
    clipper::MiniMol mol = generate_phos_minimol(0.01);

    clipper::MiniMol peaks = find_peaks(mol, m_prediction.phosphate_map);

    clipper::MiniMol middle_point = find_middlepoint(peaks, 1.5, "P");

    clipper::MiniMol refined_pos = refine_phosphates(middle_point);

    return refined_pos;

}


clipper::MiniMol PredictionBuilder::find_middlepoint(clipper::MiniMol &mol, float radius, const std::string &name) {
    clipper::MModel m_model = mol.model();
    std::vector<clipper::MAtom> results;
    std::vector<int> checked_atoms;

    clipper::MAtomNonBond non_bond = clipper::MAtomNonBond(mol, radius);

    for (int poly = 0; poly < m_model.size(); poly++) {
        for (int mon = 0; mon < m_model[poly].size(); mon++) {
            for (int atom = 0; atom < m_model[poly][mon].size(); atom++) {

                std::vector<clipper::MAtomIndexSymmetry> atom_list = non_bond.atoms_near(
                        m_model[poly][mon][atom].coord_orth(), radius);
                if (atom_list.empty()) {
                    continue;
                }

                if (std::find(checked_atoms.begin(), checked_atoms.end(), atom) != checked_atoms.end()) {
                    continue;
                }

                clipper::Coord_orth centroid_sum = {0, 0, 0};

                int count = 0;
                for (clipper::MAtomIndexSymmetry &atom_index: atom_list) {

                    if (std::find(checked_atoms.begin(), checked_atoms.end(), atom_index.atom()) ==
                        checked_atoms.end()) {
                        checked_atoms.push_back(atom_index.atom());
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


    clipper::MiniMol filtered_minimol(p_xmap.spacegroup(), p_xmap.cell());
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


clipper::Coord_grid
PredictionBuilder::ascend_gradient(clipper::Coord_grid &grid_point, clipper::Xmap<float> &xmap) const {
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


clipper::MiniMol PredictionBuilder::find_peaks(clipper::MiniMol &mol, clipper::Xmap<float> &xmap) {
    clipper::MiniMol mol_ = mol;
    clipper::MiniMol ascended_mol = clipper::MiniMol(xmap.spacegroup(), xmap.cell());

    clipper::MMonomer ascended_monomer;
    ascended_monomer.set_type("PA");
    ascended_monomer.set_id(1);

    std::vector <clipper::Coord_grid> atom_peaks;

    int max_iter = 1000;
    for (int poly = 0; poly < mol_.model().size(); poly++) {
        for (int mon = 0; mon < mol_.model()[poly].size(); mon++) {
            for (int atom = 0; atom < mol_.model()[poly][mon].size(); atom++) {

                clipper::Coord_grid m_atom_grid = mol_.model()[poly][mon][atom].coord_orth().coord_frac(
                        xmap.cell()).coord_grid(xmap.grid_sampling());

                // float phosphate_value = m_prediction.phosphate_map.get_data(m_atom_grid);
//                if (phosphate_value == 0) {
//                    continue;
//                }

                clipper::Coord_grid current_best = m_atom_grid;
                for (int i = 0; i < max_iter; i++) {

                    clipper::Coord_grid best_neighbour = ascend_gradient(current_best, xmap);
                    if (best_neighbour.u() == current_best.u() &&
                        best_neighbour.v() == current_best.v() &&
                        best_neighbour.w() == current_best.w()) {
//                        std::cout << m_atom_grid.format() << "->" << best_neighbour.format() << " " << i<< std::endl;
                        break;
                    }

                    current_best = best_neighbour;

                    if (i == max_iter - 1) {
                        std::cout << "Max iter reached" << std::endl;
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

//    return ascended_mol;

        clipper::MiniMol peak_mol = clipper::MiniMol(xmap.spacegroup(), xmap.cell());

        clipper::MMonomer peak_monomer;
        peak_monomer.set_type("X");
        peak_monomer.set_id(0);

        for (int i = 0; i < atom_peaks.size(); i++) {
            peak_monomer.insert(
                    NautilusUtil::create_atom(atom_peaks[i].coord_frac(xmap.grid_sampling()).coord_orth(xmap.cell()), i,
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


clipper::MiniMol PredictionBuilder::refine_phosphates(clipper::MiniMol &mol) const {
    clipper::MiniMol mol_ = mol;

    for (int poly = 0; poly < mol_.model().size(); poly++) {
        for (int mon = 0; mon < mol_.model()[poly].size(); mon++) {
            for (int atom = 0; atom < mol_.model()[poly][mon].size(); atom++) {
                Target_fn_refine_phosphate refinement_target(p_xmap, 0.1);
                clipper::Coord_orth refined_atom = refinement_target.refine(mol_.model()[poly][mon][atom].coord_orth());
                mol_.model()[poly][mon][atom].set_coord_orth(refined_atom);
            }
        }
    }

    return mol_;
}



void MLFind::load_library_from_file(const std::string &path) {
    na_db.add_pdb(path);
    std::cout << "Loaded library from file " << path << std::endl;
}

// clipper::MiniMol MLFind::find_with_fingerprints() {
//     clipper::MiniMol mol_new = mol;
//
//   // make a list of rotations
//   std::vector<clipper::RTop_orth> rots;
//   // make a list of rotation ops to try
//   float glim = 360.0;  // gamma
//   float blim = 180.0;  // beta
//   float alim = 360.0;  // alpha
//   // do a uniformly sampled search of orientation space
//   float anglim = clipper::Util::min( alim, glim );
//   for ( float bdeg=step/2; bdeg < 180.0; bdeg += step ) {
//     float beta = clipper::Util::d2rad(bdeg);
//     float spl = anglim/clipper::Util::intf(cos(0.5*beta)*anglim/step+1);
//     float smi = anglim/clipper::Util::intf(sin(0.5*beta)*anglim/step+1);
//     for ( float thpl=spl/2; thpl < 720.0; thpl += spl )
//       for ( float thmi=smi/2; thmi < 360.0; thmi += smi ) {
//         float adeg = clipper::Util::mod(0.5*(thpl+thmi),360.0);
//         float gdeg = clipper::Util::mod(0.5*(thpl-thmi),360.0);
//         if ( adeg <= alim && bdeg <= blim && gdeg <= glim ) {
//           float alpha = clipper::Util::d2rad(adeg);
//           float gamma = clipper::Util::d2rad(gdeg);
//           clipper::Euler_ccp4 euler( alpha, beta, gamma );
//           rots.push_back(clipper::RTop_orth(clipper::Rotation(euler).matrix()));
//         }
//       }
//   }
//
//   // get cutoff (for optimisation)
//   clipper::Map_stats stats( xmap );
//   double sigcut = stats.mean() + 1.0*stats.std_dev();
//
//   // feature search
//   if ( found_s.size() == 0 || found_p.size() == 0 ) {
//     SSfind ssfind;
//     ssfind.prep_xmap( xmap, std::max( target_sugar().radius(),
//                                     target_phosphate().radius() ) + 1.0 );
//     ssfind.prep_search( xmap );
//     found_s = ssfind.search( target_sugar().target(), rots, sigcut, 0.0 );
//     found_p = ssfind.search( target_phosphate().target(), rots, sigcut, 0.0 );
//
//     std::sort( found_s.begin(), found_s.end() );
//     std::reverse( found_s.begin(), found_s.end() );
//     std::sort( found_p.begin(), found_p.end() );
//     std::reverse( found_p.begin(), found_p.end() );
//     //std::cout << found_s.size() << "\t" << found_p.size() << std::endl;
//     //for ( int i = 0; i < clipper::Util::min( int(found_p.size()), 100 ); i++ )
//     //  std::cout << i << ":\tSgr: " << found_s[i].score << "\t" << found_s[i].rot << " " << found_s[i].trn << "\tPho: " << found_p[i].score << "\t" << found_p[i].rot << " " << found_p[i].trn << std::endl;
//   }
//
//   const clipper::Grid_sampling& grid = xmap.grid_sampling();
//   clipper::MAtomNonBond nb( mol_new, 4.0 );
//
//   // filter lists on translation
//   std::vector<SearchResult> filter_s, filter_p;
//   for ( int i = 0; i < found_s.size(); i++ ) {
//     int it = found_s[i].trn;
//     clipper::Coord_orth trn( xmap.coord_orth( grid.deindex(it).coord_map() ) );
//     std::vector<clipper::MAtomIndexSymmetry> atoms = nb( trn, 4.0 );
//     if ( atoms.size() == 0 ) filter_s.push_back( found_s[i] );
//   }
//   for ( int i = 0; i < found_p.size(); i++ ) {
//     int it = found_p[i].trn;
//     clipper::Coord_orth trn( xmap.coord_orth( grid.deindex(it).coord_map() ) );
//     std::vector<clipper::MAtomIndexSymmetry> atoms = nb( trn, 4.0 );
//     if ( atoms.size() == 0 ) filter_p.push_back( found_p[i] );
//   }
//   //std::cout << "Filter: " << mol.atom_list().size() << std::endl;
//   //std::cout << found_s.size() << " " << filter_s.size() << std::endl;
//   //std::cout << found_p.size() << " " << filter_p.size() << std::endl;
//   //std::cout << mol_new.size() << std::endl;
//
//   // build mono-units on sugars from db fragments
//   for ( int i = 0; i < std::min(int(filter_s.size()),nsugar); i++ ) {
//     std::vector<clipper::Coord_orth> v1(3), v2(3);
//     int ir = filter_s[i].rot;
//     int it = filter_s[i].trn;
//     clipper::Coord_orth trn( xmap.coord_orth( grid.deindex(it).coord_map() ) );
//     clipper::RTop_orth rtop( rots[ir].rot(), trn );
//     NucleicAcidDB::NucleicAcid na = narepr;
//     v1[0] = rtop * target_s.standard()[0]; // C3'
//     v1[1] = rtop * target_s.standard()[1]; // C1'
//     v1[2] = rtop * target_s.standard()[2]; // C4'
//     v2[0] = na.coord_c3();
//     v2[1] = na.coord_c1();
//     v2[2] = na.coord_c4();
//     clipper::RTop_orth rtdb( v2, v1 );
//     na.transform( rtdb );
//     clipper::MPolymer mp;
//     na.set_type( '?' );
//     mp.insert( na.mmonomer() );
//     mol_new.insert( mp );
//   }
//   //std::cout << mol_new.size() << std::endl;
//
//   // build bi-units on phosphates from db fragments
//   for ( int i = 0; i < std::min(int(filter_p.size()),nphosp); i++ ) {
//     std::vector<clipper::Coord_orth> v1(3), v2(3);
//     int ir = filter_p[i].rot;
//     int it = filter_p[i].trn;
//     clipper::Coord_orth trn( xmap.coord_orth( grid.deindex(it).coord_map() ) );
//     clipper::RTop_orth rtop( rots[ir].rot(), trn );
//     v1[0] = rtop * target_p.standard()[0]; // O3'
//     v1[1] = rtop * target_p.standard()[1]; // P
//     v1[2] = rtop * target_p.standard()[2]; // O5'
//     float smax = -1.0e20;
//     clipper::MPolymer mpmax;
//     for ( int j = 0; j < nadb.size()-1; j++ ) {
//       NucleicAcidDB::Chain frag = nadb.extract( j, 2 );
//       if ( frag.is_continuous() ) {
//         v2[0] = frag[0].coord_o3();
//         v2[1] = frag[1].coord_p();
//         v2[2] = frag[1].coord_o5();
//         clipper::RTop_orth rtdb( v2, v1 );
//         frag.transform( rtdb );
//         float score = ( score_sugar( xmap, frag[0] ) +
//                         score_sugar( xmap, frag[1] ) );
//         if ( score > smax ) {
//           clipper::MPolymer mpx;
//           frag[0].set_type( '?' );
//           frag[1].set_type( '?' );
//           mpx.insert( frag[0].mmonomer() );
//           mpx.insert( frag[1].mmonomer() );
//           smax = score;
//           mpmax = mpx;
//         }
//       }
//     }
//     mol_new.insert( mpmax );
//   }
//   //std::cout << mol_new.size() << std::endl;
//
//   return mol_new;
// }

clipper::MiniMol MLFind::find() {

    std::cout << "Beginning ML Find" << std::endl;
    clipper::MiniMol phosphate_positions = generate_phosphate_positions();

    TripletCoordinates triplet_coordinates = find_phosphate_triplets(phosphate_positions);

    clipper::MiniMol output_mol = {p_xmap.spacegroup(), p_xmap.cell()};
    clipper::MPolymer mp;
    mp.set_id("X");
    int count = 0;

    std::map<int, std::vector<NucleicAcidDB::NucleicAcidFull>> placed_fragments;

    char chars[] = {'-', '\\', '|', '/'};
    for (int i = 0; i < triplet_coordinates.size(); i++) {
        //    std::cout << int(i / triplets.size() * 100.0) << " %\n";
        std::cout << (100 * i / triplet_coordinates.size()) << "% " << chars[i % sizeof(chars)] << "\r";
        std::cout.flush();

        clipper::Coord_orth p1 = triplet_coordinates[i][0].second;
        clipper::Coord_orth p2 = triplet_coordinates[i][1].second;
        clipper::Coord_orth p3 = triplet_coordinates[i][2].second;

        std::vector<clipper::Coord_orth> triplet_pos = {p1, p2, p3};

        float max_score = -1e8f;
        NucleicAcidDB::ChainFull best_fragment;

        for (int j = 0; j < na_db.size() - 2; j++) {
            NucleicAcidDB::ChainFull fragment = na_db.extract(j, 3);

            if (!fragment.is_continuous()) continue;

            for (int x = 0; x < triplet_coordinates[i].size(); x++) {
                fragment[x].set_triplet_id(i * j + x);
            }

            std::vector<clipper::Coord_orth> fragment_pos = {fragment[0].P, fragment[1].P, fragment[2].P};

            clipper::RTop_orth align = clipper::RTop_orth(fragment_pos, triplet_pos);

            fragment.transform(align);
            fragment.alignment = align;

            float total_score = score_fragment(fragment, p_xmap);

            if (total_score > max_score) {
                max_score = total_score;
                best_fragment = fragment;
            }
        }

        NucleicAcidDB::ChainFull refined_fragment = refine_fragment(best_fragment, 1, 1);

        placed_fragments[triplet_coordinates[i][0].first].emplace_back(refined_fragment[0]);
        placed_fragments[triplet_coordinates[i][1].first].emplace_back(refined_fragment[1]);
        placed_fragments[triplet_coordinates[i][2].first].emplace_back(refined_fragment[2]);

    }
    output_mol.insert(mp);

    clipper::MiniMol filtered_chain = filter_and_form_chain(placed_fragments);
    clipper::MiniMol base_removed_mol = remove_bases(filtered_chain);

    clipper::MiniMol mol = organise_to_chains(base_removed_mol);

    NautilusUtil::save_minimol(mol, "mlfind.pdb");
    return mol;
}



MLFind::TripletCoordinates MLFind::find_phosphate_triplets(clipper::MiniMol &mol) {

    float radius = 8;
    clipper::MAtomNonBond m_atom_non_bond = clipper::MAtomNonBond(mol, radius);

    int p = 0;
    int m = 0;
    std::cout << mol[p][m].size() << " phosphates found" << std::endl;

    TripletCoordinates return_list;

    float target_angle = 150;
    float target_range = 30;

    for (int atom = 0; atom < mol[p][m].size(); atom++) {

        clipper::MAtom m_atom = mol.model()[p][m][atom];
        mol.model()[p][m][atom].set_id(atom);
        clipper::Coord_orth m_atom_orth = m_atom.coord_orth();
        clipper::Coord_frac m_atom_frac = m_atom_orth.coord_frac(mol.cell());
        std::vector<clipper::MAtomIndexSymmetry> atom_list = m_atom_non_bond(m_atom_orth, radius);

        for (auto &first_atom: atom_list) {

            clipper::Coord_frac first_atom_frac = mol[p][m][first_atom.atom()].coord_orth().coord_frac(mol.cell());
            first_atom_frac = first_atom_frac.symmetry_copy_near(mol.spacegroup(), mol.cell(), m_atom_frac);
            clipper::Coord_orth first_atom_orth = first_atom_frac.coord_orth(mol.cell());

            for (auto &third_atom: atom_list) {

                if (first_atom.atom() == third_atom.atom()) continue;
                if (first_atom.atom() == atom) continue;
                if (third_atom.atom() == atom) continue;

                clipper::Coord_frac third_atom_frac = mol[p][m][third_atom.atom()].coord_orth().coord_frac(mol.cell());
                third_atom_frac = third_atom_frac.symmetry_copy_near(mol.spacegroup(), mol.cell(), m_atom_frac);
                clipper::Coord_orth third_atom_orth = third_atom_frac.coord_orth(mol.cell());

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


float MLFind::score_fragment(NucleicAcidDB::ChainFull &fragment, clipper::Xmap<float> &xmap) {

    float total_score = 0.0f;
    for (int i = 0; i < fragment.size(); i++) {
        float score = score_density(fragment[i], xmap);
        fragment[i].score = score;
        total_score += score;
    }
    fragment.chain_score = total_score;
    return total_score;
}

clipper::MiniMol MLFind::filter_and_form_chain(MLFind::PossibleFragments &fragments) {

    clipper::MiniMol mol = clipper::MiniMol(p_xmap.spacegroup(), p_xmap.cell());
    clipper::MPolymer mp;

    int i = 0;
    for (auto &frag: fragments) {
        std::sort(frag.second.begin(), frag.second.end());

        auto f1 = frag.second[frag.second.size() - 1].get_mmonomer();
        f1.set_type(frag.second[frag.second.size() - 1].get_type());
        f1.set_id(i);

        auto rounded_score = static_cast<float>(static_cast<int>(frag.second[frag.second.size() - 1].score * 10.0) /
                                                10.0);
        for (int a = 0; a < f1.size(); a++) {
            f1[a].set_u_iso(rounded_score / 5);
            f1[a].set_occupancy(1.0);
        }

        i += 1;
        mp.insert(f1);
    }

    mol.model().insert(mp);
    return mol;
}


clipper::MiniMol MLFind::remove_bases(clipper::MiniMol &mol) {

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

NucleicAcidDB::ChainFull
MLFind::refine_fragment(NucleicAcidDB::ChainFull &original_fragment, float translation_range, float translation_step) {

    float best_score = -1e8;
    NucleicAcidDB::ChainFull best_chain;

    std::vector<std::vector<float>> translation_list;

    float steps = (2 * translation_range) / translation_step;

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

        Target_fn_refine_fragment fragment_refiner = Target_fn_refine_fragment(p_xmap, translated_com,
                                                                               com_fragment,
                                                                               &score_fragment, 2);
        clipper::RTop_orth refined_rtop = fragment_refiner.refine();

        NucleicAcidDB::ChainFull test_fragment = com_fragment;
        test_fragment.transform(refined_rtop);

        float refined_score = score_fragment(test_fragment, p_xmap);

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

clipper::Coord_orth MLFind::calculate_com(NucleicAcidDB::ChainFull &chain) {

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


clipper::MiniMol MLFind::organise_to_chains(clipper::MiniMol &mol) {

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

//        std::cout << connection.first << " has size " << connection.second.size() << std::endl;
        if (connection.second.empty()) {
            seed_points.emplace_back(connection.first);
//            mp.insert(mol[connection.first][0]);
        }

//        auto chains = find_n_chains(connection.first, connections, visited);
//
//        found_chains.insert(found_chains.end(), chains.begin(), chains.end());
    }

//    output_mol.insert(mp);
//    NautilusUtil::save_minimol(output_mol, "debug/1u9s/starts_of_chains.pdb");

//
//    clipper::MiniMol output_mol = {mol.spacegroup(), mol.cell()};

    for (auto& seed: seed_points) {
        Chain chain;
        find_chains(seed, connections_O3, chain);
        clipper::MPolymer mp;

        for (auto const& id: chain.ordered_chain) {
            mp.insert(mol[id][0]);
        }
        output_mol.insert(mp);
//        std::cout << "Chain info: lookup list size is " << chain.lookup_list.size() << " unordered_chain_size is " << chain.ordered_chain.size() << std::endl;
    }
    NucleicAcidTools::chain_label(output_mol, clipper::MMDBManager::PDB);

//    NautilusUtil::save_minimol(output_mol, "debug/1u9s/reordered_chain.pdb");



//    clipper::MiniMol output_mol = {mol.spacegroup(), mol.cell()};
//
//    int mon_id = 0;
//    for (int i = 0; i < found_chains.size(); i++) {
//
//        clipper::MPolymer mp;
//        if (found_chains[i].empty()) {
//            clipper::MMonomer mon = mol[i][0];
//            mon.set_id(mon_id++);
//            mp.insert(mon);
//            continue;
//        }
//
//        for (int pol_id: found_chains[i]) {
//            clipper::MMonomer mon = mol[pol_id][0];
//            mon.set_id(mon_id++);
//            mp.insert(mon);
//        }
//
//        output_mol.insert(mp);
//    }


    NucleicAcidTools::chain_label(output_mol, clipper::MMDBManager::PDB);
//    NautilusUtil::save_minimol(output_mol, "debug/1u9s/chained_find.pdb");

    return output_mol;
}

void MLFind::find_chains(int current_index, std::map<int, std::vector<int>> &connections, Chain &chain) {
    chain.lookup_list.insert(current_index);
    chain.ordered_chain.emplace_back(current_index);

    for (auto const& nearby_P: connections.at(current_index)) {
        if (chain.lookup_list.find(nearby_P) != chain.lookup_list.end())
            continue;

        find_chains(nearby_P, connections, chain);
    }

}
