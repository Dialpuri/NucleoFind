//
// Created by Jordan Dialpuri on 14/05/2025.
//

#ifndef REFINE_H
#define REFINE_H

#include "vector"
#include <clipper/clipper.h>

class SimplexOptimiser {
public:
    enum TYPE {
        NORMAL, GRADIENT
    };

    explicit SimplexOptimiser(const int n_params, const double tolerance = 0.001, const int max_cycles = 50, const TYPE type = NORMAL, const bool debug = false): tolerance_(tolerance), max_cycles_(max_cycles), type_(type), debug_mode(debug), n_params(n_params) {}

    std::vector<double> operator()(std::function<double(std::vector<double>)> &target_fn, const std::vector<std::vector<double> > &args) const;

private:
    double tolerance_, max_cycles_;
    TYPE type_;
    bool debug_mode;
    std::vector<double> params_;
    int n_params;
};

#endif //REFINE_H
