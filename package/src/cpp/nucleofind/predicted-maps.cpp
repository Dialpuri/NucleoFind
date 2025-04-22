//
// Created by Jordan Dialpuri on 06/04/2025.
//

#include "predicted-maps.h"

#include <iostream>
#include <set>
#include <unordered_set>
#include <clipper/core/coords.h>
#include <gemmi/neighbor.hpp>


gemmi::Residue NucleoFind::MapToPoints::locate_peaks(gemmi::Grid<> &xwrk, gemmi::Grid<> &xpred,
    double threshold) {
    std::vector<gemmi::Grid<>::Point> points = create_atoms_at_gridpoints(xpred, threshold);
    gemmi::Residue peaks = find_peaks(points, xpred);
    gemmi::Residue assimilated_peaks = assimilate_peaks(peaks, xpred);
    assimilated_peaks = assimilate_peaks(assimilated_peaks, xpred);
    return std::move(refine_peaks(assimilated_peaks, xpred, xwrk));;
}

std::vector<gemmi::Grid<>::Point> NucleoFind::MapToPoints::create_atoms_at_gridpoints(gemmi::Grid<> &grid, double threshold) {
    std::vector<gemmi::Grid<>::Point> points;
    for (int i = 0; i < grid.nu; i++) {
        for (int j = 0; j < grid.nv; j++) {
            for (int k = 0; k < grid.nw; k++) {
                gemmi::Grid<>::Point point = grid.get_point(i, j, k);
                if (*point.value < threshold) continue;
                points.emplace_back(point);
            }
        }
    }

    return points;
}


gemmi::Residue NucleoFind::MapToPoints::find_peaks(std::vector<gemmi::Grid<>::Point> &points,
                                                           gemmi::Grid<> &grid) {
    std::unordered_set<gemmi::Grid<>::Point, PointHash, PointEqual> visited_points = {};
    gemmi::Residue peaks;
    peaks.seqid = gemmi::SeqId("0");
    peaks.name = "X";
    for (auto & point : points) {
        gemmi::Grid<>::Point highest_point = gradient_ascent(point, grid, 0);
        if (visited_points.find(highest_point) == visited_points.end()) {
            visited_points.insert(highest_point);
            gemmi::Position highest_position = grid.point_to_position(highest_point);
            gemmi::Atom atom = create_gemmi_atom(highest_position);
            peaks.atoms.emplace_back(atom);
        }

    }
    return peaks;
}

gemmi::Residue NucleoFind::MapToPoints::assimilate_peaks(gemmi::Residue &residue, gemmi::Grid<> &grid) {
    gemmi::Model model = create_gemmi_model(residue);
    gemmi::NeighborSearch ns = {model, grid.unit_cell, 2.0};
    ns.populate();

    std::set<int> checked_atoms;
    gemmi::Residue peaks;
    peaks.seqid = gemmi::SeqId("0");
    peaks.name = "X";
    for (int a = 0; a < residue.atoms.size(); a++) {
        auto nearby = ns.find_atoms(residue.atoms[a].pos, '\0', 0, 1.5);
        if (checked_atoms.find(a) != checked_atoms.end()) {
            continue;
        }

        gemmi::Position centroid = residue.atoms[a].pos;
        int count = 1;

        for (const auto& near: nearby) {
            if (near->atom_idx == a) continue;
            auto im = symmetry_copy_near(residue.atoms[a].pos, near->pos, grid.unit_cell, grid.spacegroup);
            // auto im = grid.unit_cell.find_nearest_pbc_position(residue.atoms[a].pos, near->pos, near->image_idx);
            centroid += im;
            count++;
            checked_atoms.insert(near->atom_idx);
        }

        centroid /= count;
        peaks.atoms.emplace_back(create_gemmi_atom(centroid, std::to_string(a)));
    }
    return peaks;
}

gemmi::Residue NucleoFind::MapToPoints::refine_peaks(gemmi::Residue &residue, gemmi::Grid<> &grid,
    gemmi::Grid<> &xwrk) {
    DensityRefiner refiner(xwrk, grid);
    for (auto & atom : residue.atoms) {
        atom.pos = refiner.refine_position(atom.pos, true);
    }
    return residue;
}

gemmi::Position NucleoFind::DensityRefiner::refine_position(gemmi::Position &position, bool use_restraints) {

    auto lambda = [&](std::vector<double>& x) -> double {
        gemmi::Position current_position = {x[0], x[1], x[2]};
        if (use_restraints) {
            return this->restraint_map.interpolate_value(current_position) > 0 ? -this->xwrk.interpolate_value(current_position) : INT_MAX;
        }
        return this->xwrk.interpolate_value(current_position);
    };
    std::vector<double> initial_values = {position.x, position.y, position.z};
    std::vector<std::vector<double>> initial_simplex = {
        {initial_values[0], initial_values[1], initial_values[2]},
        {initial_values[0]+0.1, initial_values[1], initial_values[2]},
        {initial_values[0], initial_values[1]+0.1, initial_values[2]},
        {initial_values[0], initial_values[1], initial_values[2]+0.1}
    };
    std::vector<double> final_values = nelder_mead::find_min(lambda, initial_values, false, initial_simplex, 1e-8, 1e-8, 100, 100000);
    gemmi::Position final_position = {final_values[0], final_values[1], final_values[2]};
    return final_position;
}


gemmi::Grid<>::Point NucleoFind::MapToPoints::gradient_ascent(gemmi::Grid<>::Point &point, gemmi::Grid<> &grid, int iteration) {
    constexpr int max_iter = 100;
    if (iteration > max_iter) {
        return point;
    }

    gemmi::Grid<>::Point best_gridpoint = point;

    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int k = -1; k <= 1; k++) {
                gemmi::Grid<>::Point next_point = grid.get_point(point.u + i, point.v + j, point.w + k);
                if (*next_point.value > *best_gridpoint.value) {
                    best_gridpoint = next_point;
                }
            }
        }
    }


    if (PointEqual()(best_gridpoint, point)) {
        return best_gridpoint;
    }

    return gradient_ascent(best_gridpoint, grid, ++iteration);

}
