//
// Created by Jordan Dialpuri on 06/04/2025.
//

#ifndef NUCLEOFIND_H
#define NUCLEOFIND_H

#include <clipper/core/xmap.h>

#include "predicted-maps.h"
#include "fragment-library.h"

namespace NucleoFind {


    class Find {
    public:
        Find(clipper::Xmap<float>& xwrk, PredictedMaps& predicted_maps): m_xwrk(xwrk), m_predicted_maps(predicted_maps) {
            m_phosphate = const_cast<clipper::Xmap<float> *>(predicted_maps.get_phosphate_map());
            m_sugar = const_cast<clipper::Xmap<float> *>(predicted_maps.get_sugar_map());
            m_base = const_cast<clipper::Xmap<float>* >(predicted_maps.get_base_map());
        };

        void find();

    private:
        clipper::Xmap<float> m_xwrk;
        clipper::Xmap<float>* m_phosphate;
        clipper::Xmap<float>* m_sugar;
        clipper::Xmap<float>* m_base;

        PredictedMaps m_predicted_maps;
    };


}

#endif //NUCLEOFIND_H
