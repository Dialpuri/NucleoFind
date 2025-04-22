//
// Created by Jordan Dialpuri on 06/04/2025.
//

#ifndef NUCLEOFIND_UTILS_H
#define NUCLEOFIND_UTILS_H
#include <iostream>


inline gemmi::Atom create_gemmi_atom(gemmi::Position& position, std::string name = "X") {
    gemmi::Atom atom;
    atom.pos = position;
    atom.name = name;
    atom.b_iso = 20;
    atom.element = gemmi::Element("P");
    atom.occ = 1;
    return atom;
}

inline gemmi::Residue create_gemmi_residue(gemmi::Atom& atom, int seqid = 0, std::string name = "X") {
    gemmi::Residue residue;
    residue.atoms.push_back(atom);
    residue.name = name;
    residue.seqid = gemmi::SeqId(std::to_string(seqid));
    return residue;
}

inline gemmi::Model create_gemmi_model(gemmi::Residue& residue) {
    gemmi::Model model;
    model.num = 0;
    gemmi::Chain chain;
    chain.name = "A";
    chain.residues.push_back(residue);
    model.chains.push_back(chain);
    return model;
}

inline gemmi::Structure create_gemmi_structure(gemmi::Residue& residue) {
    gemmi::Structure structure;
    const gemmi::Model model = create_gemmi_model(residue);
    structure.models.push_back(model);
    return structure;
}

inline gemmi::Structure split_atoms_into_structure(gemmi::Residue& residue, gemmi::UnitCell& cell, const gemmi::SpaceGroup* spg) {
    gemmi::Structure structure;
    structure.cell = cell;
    structure.spacegroup_hm = spg->hm;
    gemmi::Chain chain = gemmi::Chain("A");

    for (int a = 0; a < residue.atoms.size(); a++) {
        gemmi::Residue r = create_gemmi_residue(residue.atoms[a], a, "X");
        chain.residues.push_back(r);
    }
    gemmi::Model model;
    model.num = 0;
    model.chains.push_back(chain);
    structure.models.push_back(model);
    return structure;
}


// Coord_frac c, cmin(*this);
// double d2, d2min(1.0e12);
// for ( int k = 0; k < spgr.num_symops(); k++ ) {
//     c = spgr.symop(k) * (*this);
//     c = c.lattice_copy_near( n );
//     d2 = ( c - n ).lengthsq( cell );
//     if ( d2 < d2min ) {
//         d2min = d2;
//         cmin = c;
//     }
// }
// return cmin;

inline gemmi::Fractional lattice_copy_zero(gemmi::Fractional& target) {
    return gemmi::Fractional(target.x - rint(target.x), target.y - rint(target.y), target.z - rint(target.z));
}

inline gemmi::Fractional lattice_copy_near(gemmi::Fractional& target, gemmi::Fractional& reference) {
    auto diff = target-reference;
    return lattice_copy_zero(diff) + reference;
}

inline gemmi::Position symmetry_copy_near(gemmi::Position& target, gemmi::Position& reference, gemmi::UnitCell& cell, const gemmi::SpaceGroup* spg) {
    gemmi::Fractional target_f = cell.fractionalize(target);
    gemmi::Fractional reference_f = cell.fractionalize(reference);
    gemmi::Fractional cmin = target_f;
    gemmi::Fractional c;
    double d2min = 1e12;
    double d2;
    auto symops = spg->operations().sym_ops;

    for (int k = 0; k < symops.size(); k++) {
        auto translated_v = symops[k].apply_to_xyz({target_f.x, target_f.y, target_f.z});
        gemmi::Fractional translated_f = {translated_v[0], translated_v[1], translated_v[2]};
        c = lattice_copy_near(translated_f, reference_f);
        d2 = (cell.orthogonalize(c) - cell.orthogonalize(reference_f)).length_sq();
        if (d2 < d2min) {
            d2min = d2;
            cmin = c;
        }
    }

    return cell.orthogonalize(cmin);
}

#endif //NUCLEOFIND_UTILS_H
