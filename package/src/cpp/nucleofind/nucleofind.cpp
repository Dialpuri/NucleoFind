//
// Created by Jordan Dialpuri on 06/04/2025.
//

#include "nucleofind.h"

#include "backbone-tracing.h"
#include "src/cpp/nautilus-tools.h"
#include "src/cpp/nautilus-util.h"

clipper::MiniMol NucleoFind::Find::find(clipper::MiniMol &mol_wrk) {
    clipper::MiniMol phosphate_peaks = MapToPoints::locate_peaks(m_xwrk, *m_phosphate, 0.1, true);
    // NautilusUtil::save_minimol(phosphate_peaks, "phosphate_peaks.pdb");
    BackboneTracer b = {phosphate_peaks, m_xwrk, m_predicted_maps};
    clipper::MiniMol find_result = b.build();
    return std::move(aggregate(find_result, mol_wrk));
}

clipper::MiniMol NucleoFind::Find::aggregate(clipper::MiniMol &find_result, clipper::MiniMol &mol) {

    clipper::MAtomNonBond nb = {mol, 4};

    std::set<std::tuple<int,int>> to_remove = {};
    for (int c = 0; c < find_result.size(); c++) {
        for (int r = 0; r < find_result[c].size(); r++) {
            for (int a = 0; a < find_result[c][r].size(); a++) {
                auto atoms_near = nb(find_result[c][r][a].coord_orth(), 3);
                for (auto& atom: atoms_near) {
                    to_remove.insert({atom.polymer(), atom.monomer()});
                }
            }
        }
    }

    for (int c = 0; c < mol.size(); c++) {
        clipper::MPolymer mp;;
        mp.set_id(mol[c].id());
        for (int r = 0; r < mol[c].size(); r++) {
            if (to_remove.count({c, r}) != 0) continue;
            mp.insert(mol[c][r]);
        }
        find_result.insert(mp);
    }

    find_result = NucleicAcidTools::flag_chains(find_result);

    return find_result;
}
