//
// Created by Jordan Dialpuri on 06/04/2025.
//

#include "predicted-maps.h"
#include "refine.h"

#include <iostream>
#include <set>
#include <unordered_set>
#include <clipper/core/coords.h>
#include <clipper/core/map_interp.h>
#include <clipper/minimol/minimol_utils.h>
#include <gemmi/neighbor.hpp>

#include "src/cpp/nautilus-util.h"


clipper::MiniMol NucleoFind::MapToPoints::locate_peaks(clipper::Xmap<float>& xwrk, clipper::Xmap<float>& xpred,
                                                       double threshold, bool refine) {
    std::vector<clipper::Coord_grid> points = create_atoms_at_gridpoints(xpred, threshold);
    clipper::MiniMol peaks = find_peaks(points, xpred);
    clipper::MiniMol assimilated_peaks = assimilate_peaks(peaks, xpred);
    // assimilated_peaks = assimilate_peaks(assimilated_peaks, xpred);
    if (!refine)
        return std::move(assimilated_peaks);
    return std::move(refine_peaks(assimilated_peaks, xpred, xwrk));;
}

clipper::MiniMol NucleoFind::MapToPoints::create_mol_at_gridpoints(clipper::Xmap<float> &grid, double threshold) {
    std::vector<clipper::Coord_grid> points = create_atoms_at_gridpoints(grid, threshold);
    clipper::MMonomer peaks;
    peaks.set_seqnum('0');
    peaks.set_id("X");
    for (auto & point : points) {
        clipper::Coord_orth point_orth = point.coord_frac(grid.grid_sampling()).coord_orth(grid.cell());
        clipper::Atom atom = create_clipper_atom(point_orth);
        peaks.insert(atom);
    }
    return create_clipper_minimol(peaks, grid);
}

std::vector<clipper::Coord_grid> NucleoFind::MapToPoints::create_atoms_at_gridpoints(clipper::Xmap<float> &grid, double threshold) {
    std::vector<clipper::Coord_grid> points;

    clipper::Grid_sampling grid_sampling = grid.grid_sampling();
    const auto g0 = clipper::Coord_grid(0, 0, 0);
    auto g1 = clipper::Coord_grid(grid_sampling.nu(), grid_sampling.nv(), grid_sampling.nw());
    const auto i0 = clipper::Xmap_base::Map_reference_coord(grid, g0);

    for (clipper::Xmap_base::Map_reference_coord iu = i0; iu.coord().u() <= g1.u(); iu.next_u()) {
        for (clipper::Xmap_base::Map_reference_coord iv = iu; iv.coord().v() <= g1.v(); iv.next_v()) {
            for (clipper::Xmap_base::Map_reference_coord iw = iv; iw.coord().w() <= g1.w(); iw.next_w()) {
                if (grid[iw] >= threshold) {
                    points.emplace_back(iw.coord());
                }
            }
        }
    }
    return points;
}


clipper::MiniMol NucleoFind::MapToPoints::find_peaks(std::vector<clipper::Coord_grid> &points,
                                                   clipper::Xmap<float> &grid) {
    std::unordered_set<clipper::Coord_grid, GridPointHash, GridPointEqual> visited_points = {};
    clipper::MMonomer peaks;
    peaks.set_seqnum('0');
    peaks.set_id("X");
    for (auto & point : points) {
        if (clipper::Coord_grid highest_point = gradient_ascent(point, grid, 0); visited_points.find(highest_point) == visited_points.end()) {
            visited_points.insert(highest_point);
            clipper::Coord_orth highest_point_orth = highest_point.coord_frac(grid.grid_sampling()).coord_orth(grid.cell());
            clipper::Atom atom = create_clipper_atom(highest_point_orth);
            peaks.insert(atom);
        }
    }
    return create_clipper_minimol(peaks, grid);
}

clipper::MiniMol NucleoFind::MapToPoints::assimilate_peaks(clipper::MiniMol & mol, clipper::Xmap<float> &grid) {
    float radius = 2.0;
    clipper::MAtomNonBond nb = clipper::MAtomNonBond(mol, radius);

    // We can assume that the MiniMol has one MPolymer, with one MMonomer, containing all the points.
    clipper::MMonomer peaks = mol[0][0];

    std::unordered_set<int> checked_points;

    std::vector<clipper::Coord_orth> results;
    for (int atom = 0; atom < peaks.size(); atom++) {
        auto atom_list = nb.atoms_near(peaks[atom].coord_orth(), radius);
        if (atom_list.empty()) {continue;}
        if (checked_points.find(atom) != checked_points.end()) {continue;}

        int count = 0;
        clipper::Coord_orth centroid_sum = {0, 0, 0};

        for (auto& atom_index: atom_list) {
            if (checked_points.find(atom_index.atom()) == checked_points.end()) {
                checked_points.insert(atom_index.atom());
            }

            if (clipper::Coord_orth::length(peaks[atom].coord_orth(),  peaks[atom_index.atom()].coord_orth()) > radius) {
                continue;
            }

            centroid_sum += peaks[atom].coord_orth();
            count++;
        }

        if (count > 0) {
            auto centroid = clipper::Coord_orth(centroid_sum * (1.0 / count));
            results.emplace_back(centroid);
        }
    }

    std::vector<clipper::MAtom> atoms;
    atoms.reserve(results.size());
    for (int r = 0; r < results.size(); r++) {
        atoms.emplace_back(create_clipper_atom(results[r], std::to_string(r)));
    }
    return create_clipper_minimol(atoms, grid);
}

clipper::MiniMol NucleoFind::MapToPoints::refine_peaks(clipper::MiniMol &mol, clipper::Xmap<float> &grid,
                                                     clipper::Xmap<float> &xwrk) {
    DensityRefiner refiner(xwrk, grid);
    // We can assume that the MiniMol has one MPolymer, with one MMonomer, containing all the points.
    clipper::MMonomer peaks = mol[0][0];
    for (int a = 0; a < peaks.size(); a++) {
        auto x = refiner.refine_position(peaks[a].coord_orth(), true);
        mol[0][0][a].set_coord_orth(x);
    }
    return mol;
}

clipper::Coord_orth NucleoFind::DensityRefiner::refine_position(const clipper::Coord_orth &position, bool use_restraints) {
    std::function<double(std::vector<double>)> lambda = [&](std::vector<double> x) {
        clipper::Coord_orth current_position = {x[0], x[1], x[2]};
        clipper::Coord_frac current_position_f = current_position.coord_frac(this->xwrk.cell());
        if (use_restraints) {
            return this->restraint_map.interp<clipper::Interp_cubic>(current_position_f) > 0 ? -this->xwrk.interp<clipper::Interp_cubic>(current_position_f) : INT_MAX;
        }
        return this->xwrk.interp<clipper::Interp_cubic>(current_position_f);
    };

    std::vector<double> initial_values = {position.x(), position.y(), position.z()};
    std::vector<std::vector<double>> initial_simplex = {
        {initial_values[0], initial_values[1], initial_values[2]},
        {initial_values[0]+0.1, initial_values[1], initial_values[2]},
        {initial_values[0], initial_values[1]+0.1, initial_values[2]},
        {initial_values[0], initial_values[1], initial_values[2]+0.1}
    };
    SimplexOptimiser simplex = SimplexOptimiser(3);
    std::vector<double> final_values = simplex(lambda, initial_simplex);
    clipper::Coord_orth final_position = {final_values[0], final_values[1], final_values[2]};

    return final_position;
}


clipper::Coord_grid NucleoFind::MapToPoints::gradient_ascent(clipper::Coord_grid &point, clipper::Xmap<float> &grid, int iteration) {
    constexpr int max_iter = 100;
    if (iteration > max_iter) {
        return point;
    }

    clipper::Coord_grid best_gridpoint = point;
    float best_value = grid.get_data(point);

    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int k = -1; k <= 1; k++) {
                clipper::Coord_grid next_point = {point.u() + i, point.v() + j, point.w() + k};
                if (const float value = grid.get_data(next_point); value > best_value) {
                    best_value = value;
                    best_gridpoint = next_point;
                }
            }
        }
    }

    if (GridPointEqual()(best_gridpoint, point)) {
        return best_gridpoint;
    }

    return gradient_ascent(best_gridpoint, grid, ++iteration);

}
