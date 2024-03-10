//
// Created by jordan on 06/06/23.
//

#ifndef NAUTILUS_NAUTILUS_REFINE_H
#define NAUTILUS_NAUTILUS_REFINE_H

#include <utility>
#include <vector>
#include <clipper/clipper.h>
#include <functional>
#include "nucleicacid_db.h"
#include "nautilus-predict.h"

/*! Abstract base class for zero-th order function. */
class Target_fn_order_zero {
public:
    Target_fn_order_zero() {}

    virtual ~Target_fn_order_zero() {}

    virtual int num_params() const = 0;

    virtual double operator()(const std::vector<double> &args) const = 0;
};


/*! Simplex optimiser. */
class Optimiser_simplex {
public:
    enum TYPE {
        NORMAL, GRADIENT
    };

    Optimiser_simplex(double tolerance = 0.001, int max_cycles = 50, TYPE type =
    NORMAL, bool debug = false);

    std::vector<double> operator()(const Target_fn_order_zero &target_fn, const
    std::vector<std::vector<double> > &args) const;

    void debug() { debug_mode = true; }

private:
    double tolerance_, max_cycles_;
    TYPE type_;
    bool debug_mode;
    std::vector<double> params_;
};

class Target_fn_refine_phosphate : Target_fn_order_zero {
public:

    Target_fn_refine_phosphate() = default;

    Target_fn_refine_phosphate(const clipper::Xmap<float> &xmap, float step) : m_xmap(&xmap), m_step(step) {};

    inline int num_params() const override { return 3; }

    double operator()(const std::vector<double> &args) const;

    clipper::Coord_orth refine(const clipper::Coord_orth &coord);

private:
    const clipper::Xmap<float> *m_xmap{};
    float m_step = 0.1;
};


class Target_fn_refine_fragment : Target_fn_order_zero {
public:

    // typedef float (*score_function)(NucleicAcidDB::ChainFull&, clipper::Xmap<float>&);

    Target_fn_refine_fragment() = default;

    Target_fn_refine_fragment(clipper::Xmap<float> &xmap, clipper::Vec3<> translation,
                              NucleicAcidDB::ChainFull &chain, float step)
            : m_xmap(&xmap), m_translation(translation), m_chain(chain), m_step(step) {
    };

    static float score_density(NucleicAcidDB::NucleicAcidFull &chain, clipper::Xmap<float> &xmap) {
        float score = 0.0f;

        if (!chain.O5p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.O5p1.coord_frac(xmap.cell()));
        if (!chain.C5p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.C5p1.coord_frac(xmap.cell()));
        if (!chain.C4p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.C4p1.coord_frac(xmap.cell()));
        if (!chain.O4p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.O4p1.coord_frac(xmap.cell()));
        if (!chain.C3p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.C3p1.coord_frac(xmap.cell()));
        if (!chain.O3p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.O3p1.coord_frac(xmap.cell()));
        if (!chain.C2p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.C2p1.coord_frac(xmap.cell()));
        if (!chain.C1p1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.C1p1.coord_frac(xmap.cell()));
        if (!chain.N1_1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.N1_1.coord_frac(xmap.cell()));
        if (!chain.N9_1.is_null()) score += xmap.interp<clipper::Interp_cubic>(chain.N9_1.coord_frac(xmap.cell()));


        return score;
    }


    float score_fragment(NucleicAcidDB::ChainFull &fragment, clipper::Xmap<float> &xmap) const {
        float total_score = 0.0f;
        for (int i = 0; i < fragment.size(); i++) {
            float score = score_density(fragment[i], xmap);
            fragment[i].score = score;
            total_score += score;
        }
        fragment.chain_score = total_score;
        return total_score;
    }

    inline int num_params() const override { return 3; }

    double operator()(const std::vector<double> &args) const;

    clipper::RTop_orth refine();

private:
    clipper::Xmap<float> *m_xmap{};
    clipper::Vec3<> m_translation;
    NucleicAcidDB::ChainFull m_chain;

    float m_step = clipper::Util::d2rad(1);
};

class Target_fn_refine_fragment_trn : Target_fn_order_zero {
public:

    Target_fn_refine_fragment_trn() = default;

    Target_fn_refine_fragment_trn(const clipper::Xmap<float> &xmap, clipper::Vec3<> translation,
                              NucleicAcidDB::ChainFull &chain, float step)
            : m_xmap(&xmap), m_translation(translation), m_chain(chain), m_step(step) {};

    inline int num_params() const override { return 3; }

    double operator()(const std::vector<double> &args) const;

    clipper::RTop_orth refine();

private:
    const clipper::Xmap<float> *m_xmap{};
    clipper::Vec3<> m_translation;
    NucleicAcidDB::ChainFull m_chain;
    float m_step = clipper::Util::d2rad(1);
};

#endif //NAUTILUS_NAUTILUS_REFINE_H
