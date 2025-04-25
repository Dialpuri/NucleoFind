//
// Created by Jordan Dialpuri on 06/04/2025.
//

#include "nucleofind.h"

#include "backbone-tracing.h"
#include "src/cpp/nautilus-tools.h"
#include "src/cpp/nautilus-util.h"

clipper::MiniMol NucleoFind::Find::find(clipper::MiniMol &mol_wrk) {
    clipper::MiniMol phosphate_peaks = MapToPoints::locate_peaks(m_xwrk, *m_phosphate, 0.1, true);
    NautilusUtil::save_minimol(phosphate_peaks, "phosphate_peaks.pdb");
    BackboneTracer b = {phosphate_peaks, m_xwrk, m_predicted_maps};
    clipper::MiniMol find_result = b.build();
    return std::move(aggregate(find_result, mol_wrk));
}

clipper::MiniMol NucleoFind::Find::aggregate(clipper::MiniMol &find_result, clipper::MiniMol &mol) {

    if (mol.size() == 0) {
        NucleicAcidTools::chain_label(find_result, clipper::MMDBManager::Default);
        NucleicAcidTools::residue_label(find_result);
        return NucleicAcidTools::flag_chains(find_result);
    }
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

    std::set<std::string> backbone_atoms = {"CA", "C", "O", "N"};
    for (int c = 0; c < mol.size(); c++) {
        clipper::MPolymer mp;;
        for (int r = 0; r < mol[c].size(); r++) {
            if (to_remove.count({c, r}) != 0) {
                clipper::MMonomer mon;
                mon.set_id( mol[c][r].id());
                mon.set_type(mol[c][r].type());
                for (int a = 0; a < mol[c][r].size(); a++) {
                    if (backbone_atoms.count(mol[c][r][a].name().trim()) > 0) {
                        mon.insert(mol[c][r][a]);
                    }
                }
                if (mon.size() > 0) {
                    mp.insert(mon);
                    continue;
                }
            };
            mp.insert(mol[c][r]);
        }
        if (mp.size() == 0) continue;
        find_result.insert(mp);
    }

    NucleicAcidTools::chain_label(find_result, clipper::MMDBManager::Default);
    NucleicAcidTools::residue_label(find_result);
    find_result = NucleicAcidTools::flag_chains(find_result);

    return find_result;
}
