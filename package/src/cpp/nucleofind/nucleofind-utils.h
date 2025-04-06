//
// Created by Jordan Dialpuri on 06/04/2025.
//

#ifndef NUCLEOFIND_UTILS_H
#define NUCLEOFIND_UTILS_H


inline gemmi::Atom create_gemmi_atom(gemmi::Position& position, std::string name = "X") {
    gemmi::Atom atom;
    atom.pos = position;
    atom.name = name;
    atom.b_iso = 20;
    atom.element = gemmi::Element("P");
    atom.occ = 1;
    return atom;
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

#endif //NUCLEOFIND_UTILS_H
