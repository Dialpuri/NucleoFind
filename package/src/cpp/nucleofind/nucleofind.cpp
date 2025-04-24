//
// Created by Jordan Dialpuri on 06/04/2025.
//

#include "nucleofind.h"

#include "backbone-tracing.h"
#include "src/cpp/nautilus-util.h"

void NucleoFind::Find::find() {
    clipper::MiniMol phosphate_peaks = MapToPoints::locate_peaks(m_xwrk, *m_phosphate, 0.1, false);
    NautilusUtil::save_minimol(phosphate_peaks, "phosphate_peaks.pdb");


    BackboneTracer b = {phosphate_peaks, m_xwrk, m_predicted_maps};
    b.build();

}
