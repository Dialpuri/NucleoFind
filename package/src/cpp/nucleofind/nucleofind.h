//
// Created by Jordan Dialpuri on 06/04/2025.
//

#ifndef NUCLEOFIND_H
#define NUCLEOFIND_H

#include <clipper/core/xmap.h>

#include "predicted-maps.h"

namespace NucleoFind {


    class Find {
    public:
        Find(gemmi::Grid<>& xwrk, PredictedMaps& predicted_maps): m_xwrk(xwrk), m_predicted_maps(predicted_maps) {
            m_phosphate = const_cast<gemmi::Grid<> *>(predicted_maps.get_phosphate_map());
            m_sugar = const_cast<gemmi::Grid<> *>(predicted_maps.get_sugar_map());
            m_base = const_cast<gemmi::Grid<> *>(predicted_maps.get_base_map());
        };

        void find();

    private:
        gemmi::Grid<> m_xwrk;
        gemmi::Grid<>* m_phosphate;
        gemmi::Grid<>* m_sugar;
        gemmi::Grid<>* m_base;

        PredictedMaps m_predicted_maps;
    };


}

#endif //NUCLEOFIND_H
