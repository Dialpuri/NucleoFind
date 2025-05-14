//
// Created by Jordan Dialpuri on 06/04/2025.
//

#ifndef NUCLEOFIND_UTILS_H
#define NUCLEOFIND_UTILS_H
#include <iostream>


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
