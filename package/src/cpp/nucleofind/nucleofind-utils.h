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


inline clipper::MAtom create_clipper_atom(clipper::Coord_orth& position, std::string name = "X") {
    clipper::MAtom atom;
    atom.set_coord_orth(position);
    atom.set_element("P");
    atom.set_occupancy(1);
    atom.set_u_iso(clipper::Util::b2u(20));
    atom.set_name(name);
    return atom;
}

inline clipper::MMonomer create_clipper_monomer(std::vector<clipper::MAtom>& atoms, int seqid = 0, std::string name = "X") {
    clipper::MMonomer monomer;
    monomer.set_id(seqid);
    monomer.set_seqnum(seqid);
    monomer.set_type(name);
    for (int a = 0; a < atoms.size(); a++) {
        atoms[a].set_id(a);
        monomer.insert(atoms[a]);
    }
    return monomer;
}

inline clipper::MPolymer create_clipper_polymer(clipper::MMonomer& monomer) {
    clipper::MPolymer mpol;
    mpol.insert(monomer);
    mpol.set_id("A");
    return mpol;
}

inline clipper::MPolymer create_clipper_polymer(std::vector<clipper::MMonomer>& monomers, std::string id = "A") {
    clipper::MPolymer mpol;
    for (auto& monomer: monomers) {
        mpol.insert(monomer);
    }
    mpol.set_id(id);
    return mpol;
}


inline clipper::MiniMol create_clipper_minimol(clipper::MMonomer& monomer, clipper::Xmap<float>& grid) {
    clipper::MiniMol mol = {grid.spacegroup(), grid.cell()};
    clipper::MPolymer mpol = create_clipper_polymer(monomer);
    mol.insert(mpol);
    return mol;
}

inline clipper::MiniMol create_clipper_minimol(std::vector<clipper::MAtom>& atoms, clipper::Xmap<float>& grid) {
    clipper::MiniMol mol = {grid.spacegroup(), grid.cell()};
    clipper::MMonomer monomer = create_clipper_monomer(atoms);
    clipper::MPolymer mpol = create_clipper_polymer(monomer);
    mol.insert(mpol);
    return mol;
}

inline clipper::MiniMol create_clipper_minimol(std::vector<clipper::MMonomer>& monomers, clipper::Xmap<float>& grid) {
    clipper::MiniMol mol = {grid.spacegroup(), grid.cell()};
    clipper::MPolymer mpol = create_clipper_polymer(monomers);
    mol.insert(mpol);
    return mol;
}

inline std::string nth_letter(int index) {
    std::string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    return std::string(1, alphabet[index]);
}
// inline gemmi::Structure split_atoms_into_structure(gemmi::Residue& residue, gemmi::UnitCell& cell, const gemmi::SpaceGroup* spg) {
//     gemmi::Structure structure;
//     structure.cell = cell;
//     structure.spacegroup_hm = spg->hm;
//     gemmi::Chain chain = gemmi::Chain("A");
//
//     for (int a = 0; a < residue.atoms.size(); a++) {
//         gemmi::Residue r = create_gemmi_residue(residue.atoms[a], a, "X");
//         chain.residues.push_back(r);
//     }
//     gemmi::Model model;
//     model.num = 0;
//     model.chains.push_back(chain);
//     structure.models.push_back(model);
//     return structure;
// }


#endif //NUCLEOFIND_UTILS_H
