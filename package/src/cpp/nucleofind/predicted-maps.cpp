//
// Created by Jordan Dialpuri on 06/04/2025.
//

#include "predicted-maps.h"

#include <iostream>
#include <set>
#include <unordered_set>
#include <gemmi/neighbor.hpp>


std::vector<gemmi::Grid<>::Point> NucleoFind::PredictedMapToPoint::create_atoms_at_gridpoints(gemmi::Grid<> &grid, double threshold) {
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


gemmi::Residue NucleoFind::PredictedMapToPoint::find_peaks(std::vector<gemmi::Grid<>::Point> &points,
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

gemmi::Residue NucleoFind::PredictedMapToPoint::assimilate_peaks(gemmi::Residue &residue, gemmi::Grid<> &grid) {
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
            centroid += near->pos;
            count++;
            checked_atoms.insert(near->atom_idx);
        }

        centroid /= count;
        float length = (centroid - residue.atoms[a].pos).length();
        std::cout << a << " " << length << std::endl;
        peaks.atoms.emplace_back(create_gemmi_atom(centroid, std::to_string(a)));

    }

    return peaks;
}


gemmi::Grid<>::Point NucleoFind::PredictedMapToPoint::gradient_ascent(gemmi::Grid<>::Point &point, gemmi::Grid<> &grid, int iteration) {
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
