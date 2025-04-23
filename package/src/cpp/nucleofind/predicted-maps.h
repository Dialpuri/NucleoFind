//
// Created by Jordan Dialpuri on 06/04/2025.
//

#ifndef PREDICTED_MAPS_H
#define PREDICTED_MAPS_H

#include <clipper/core/xmap.h>
#include <clipper/minimol/minimol.h>
#include <gemmi/ccp4.hpp>
#include <gemmi/model.hpp>
#include "nucleofind-utils.h"
#include "nelder-mead.h"

namespace NucleoFind {
    class PredictedMaps {
    public:
        PredictedMaps(const clipper::Xmap<float>& phosphate_map, const clipper::Xmap<float>& sugar_map,
                      const clipper::Xmap<float>& base_map)
                : phosphate(phosphate_map), sugar(sugar_map), base(base_map) {}


        [[nodiscard]] const clipper::Xmap<float>* get_phosphate_map() const {
            return &phosphate;
        }

        [[nodiscard]] const clipper::Xmap<float>* get_sugar_map() const {
            return &sugar;
        }

        [[nodiscard]] const clipper::Xmap<float>* get_base_map() const {
            return &base;
        }


    private:
        clipper::Xmap<float> phosphate;
        clipper::Xmap<float> sugar;
        clipper::Xmap<float> base;
    };


    struct GridPointHash {
        std::size_t operator()(const clipper::Coord_grid& t) const noexcept {
            std::size_t h1 = std::hash<int>()(t.u());
            std::size_t h2 = std::hash<int>()(t.v());
            std::size_t h3 = std::hash<int>()(t.w());
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };

    struct GridPointEqual {
        bool operator()(const clipper::Coord_grid& p1, const clipper::Coord_grid& p2) const noexcept {
            return p1.u() == p2.u() && p1.v() == p2.v() && p1.w() == p2.w();
        }
    };

    class MapToPoints {
    public:
        static clipper::MiniMol locate_peaks(clipper::Xmap<float>& xwrk, clipper::Xmap<float>& xpred, double threshold, bool refine);

        static clipper::MiniMol create_mol_at_gridpoints(clipper::Xmap<float> &grid, double threshold);

    private:
        static std::vector<clipper::Coord_grid> create_atoms_at_gridpoints(clipper::Xmap<float> &grid, double threshold);
        static clipper::MiniMol find_peaks(std::vector<clipper::Coord_grid> &points, clipper::Xmap<float> &grid);
        static clipper::Coord_grid gradient_ascent(clipper::Coord_grid &point, clipper::Xmap<float> &grid, int iteration);
        static clipper::MiniMol assimilate_peaks(clipper::MiniMol & mol, clipper::Xmap<float> &grid);
        static clipper::MiniMol refine_peaks(clipper::MiniMol &mol, clipper::Xmap<float> &grid, clipper::Xmap<float> &xwrk);
    };

    class DensityRefiner {
    public:
        DensityRefiner(clipper::Xmap<float>& xwrk, clipper::Xmap<float>& restraint_map): xwrk(xwrk), restraint_map(restraint_map) {};
        clipper::Coord_orth refine_position(const clipper::Coord_orth &position, bool use_restraints=true);
    private:
        clipper::Xmap<float>& xwrk;
        clipper::Xmap<float>& restraint_map;
    };
}

#endif //PREDICTED_MAPS_H
