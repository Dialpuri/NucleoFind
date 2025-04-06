//
// Created by Jordan Dialpuri on 06/04/2025.
//

#ifndef PREDICTED_MAPS_H
#define PREDICTED_MAPS_H

#include <gemmi/ccp4.hpp>
#include <gemmi/model.hpp>
#include "nucleofind-utils.h"
#include "nelder-mead.h"

namespace NucleoFind {
    class PredictedMaps {
    public:
        PredictedMaps(const gemmi::Grid<> &phosphate_map, const gemmi::Grid<> &sugar_map,
                      const gemmi::Grid<> &base_map)
                : phosphate(phosphate_map), sugar(sugar_map), base(base_map) {}


        [[nodiscard]] const gemmi::Grid<>* get_phosphate_map() const {
            return &phosphate;
        }

        [[nodiscard]] const gemmi::Grid<>* get_sugar_map() const {
            return &sugar;
        }

        [[nodiscard]] const gemmi::Grid<>* get_base_map() const {
            return &base;
        }


    private:
        gemmi::Grid<> phosphate;
        gemmi::Grid<> sugar;
        gemmi::Grid<> base;
    };


    struct PointHash {
        std::size_t operator()(const gemmi::Grid<>::Point& t) const noexcept {
            std::size_t h1 = std::hash<int>()(t.u);
            std::size_t h2 = std::hash<int>()(t.v);
            std::size_t h3 = std::hash<int>()(t.w);
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };

    struct PointEqual {
        bool operator()(const gemmi::Grid<>::Point& p1, const gemmi::Grid<>::Point& p2) const noexcept {
            return p1.u == p2.u && p1.v == p2.v && p1.w == p2.w;
        }
    };

    class MapToPoints {
    public:
        static gemmi::Residue locate_peaks(gemmi::Grid<>& xwrk, gemmi::Grid<>& xpred, double threshold);

    private:
        static std::vector<gemmi::Grid<>::Point> create_atoms_at_gridpoints(gemmi::Grid<> &grid, double threshold);
        static gemmi::Residue find_peaks(std::vector<gemmi::Grid<>::Point> &points, gemmi::Grid<> &grid);
        static gemmi::Grid<>::Point gradient_ascent(gemmi::Grid<>::Point& point, gemmi::Grid<>& grid, int iteration);
        static gemmi::Residue assimilate_peaks(gemmi::Residue &residue, gemmi::Grid<> &grid);
        static gemmi::Residue refine_peaks(gemmi::Residue &residue, gemmi::Grid<> &grid, gemmi::Grid<>& xwrk);
    };

    class DensityRefiner {
    public:
        DensityRefiner(gemmi::Grid<>& xwrk, gemmi::Grid<>& restraint_map): xwrk(xwrk), restraint_map(restraint_map) {};
        gemmi::Position refine_position(gemmi::Position& position, bool use_restraints=true);
    private:
        gemmi::Grid<>& xwrk;
        gemmi::Grid<>& restraint_map;
    };
}

#endif //PREDICTED_MAPS_H
