//
// Created by Jordan Dialpuri on 06/04/2025.
//

#include "nucleofind.h"

#include "backbone-tracing.h"
#include "src/cpp/nautilus-util.h"

void NucleoFind::Find::find() {
    gemmi::Residue phosphate_peaks = MapToPoints::locate_peaks(m_xwrk, *m_phosphate, 0.1);

    auto s = create_gemmi_structure(phosphate_peaks);
    s.cell = m_xwrk.unit_cell;

    std::ofstream fout;
    fout.open("phosphate_peaks.pdb");
    gemmi::write_pdb(s, fout);
    fout.close();

    BackboneGraph b = {phosphate_peaks, m_xwrk};

    clipper::MAtomNonBond ma;

}
