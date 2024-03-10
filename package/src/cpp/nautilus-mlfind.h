//
// Created by Jordan Dialpuri on 13/09/2023.
//

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>

#include <string>
#include <map>
#include <unordered_set>
#include <algorithm>

#include "nautilus-util.h"
#include "nautilus-predict.h"
#include "nautilus-refine.h"
#include "nautilus-join.h"
#include "nucleicacid_db.h"
#include "nautilus-tools.h"

#ifndef NAUTILUS_NAUTILUS_MLFIND_H
#define NAUTILUS_NAUTILUS_MLFIND_H

class PredictionBuilder {
public:
    PredictionBuilder(Predictions &prediction) : m_prediction(prediction) {};

    /// @brief Generic function which scores all of ChainFull with supplied xmap
    /// @param chain
    /// @param xmap
    /// @return
    static float score_density(NucleicAcidDB::NucleicAcidFull &chain, clipper::Xmap<float> &xmap);

protected:
    /// @brief Transform all density over a specific atom into a MiniMol with 1 MModel, 1 MPolymer and 1 MMonomer
    /// @param xmap
    /// @param threshold_value
    /// @return
    static clipper::MiniMol density_to_atom(clipper::Xmap<float> &xmap, float threshold_value);

    /// @brief Using predicted phosphate map values over threshold value, generate MiniMol of positions
    /// @param threshold_value
    /// @return
    clipper::MiniMol generate_phos_minimol(float threshold_value = 0.9);

    /// @brief Using gradient ascent based followed by distance clustering algorithms to group all phosphate atoms into single points
    /// @return
    clipper::MiniMol generate_phosphate_positions();


    /// @brief Find peak density values in phosphate atoms
    /// @param mol
    /// @return
    clipper::MiniMol find_peaks(clipper::MiniMol& mol, clipper::Xmap<float>& xmap);

    // TODO: Make this function generic
    /// @brief Find the next coord_grid position which ascends the gradient of the predicted phosphate map
    /// @param grid_point
    /// @return
    clipper::Coord_grid ascend_gradient(clipper::Coord_grid& grid_point, clipper::Xmap<float>& xmap) const;

    /// @brief Using distance based algorithm, join atoms within a specified radius
    /// @param mol
    /// @param radius
    /// @param name
    /// @return
    clipper::MiniMol find_middlepoint(clipper::MiniMol &mol, float radius, const std::string &name);

    /// @brief Refine phosphate positions into experimental density using simplex refinement
    /// @param mol
    /// @return
    clipper::MiniMol refine_phosphates(clipper::MiniMol &mol) const;

public:
    clipper::Xmap<float> p_xmap;
    NucleicAcidDB::Chain na_db;

protected:
    clipper::CCP4MTZfile m_mtz;
    Predictions m_prediction;
};


struct Chain {
    std::vector<int> ordered_chain;
    std::unordered_set<int> lookup_list;
};

class MLFind : public PredictionBuilder {
public:
    typedef std::vector<int> Triplet;
    typedef std::vector<Triplet> Triplets;
    typedef std::vector<std::vector<std::pair<int, clipper::Coord_orth>>> TripletsCoord;

    typedef std::vector<std::pair<int, clipper::Coord_orth>> TripletCoordinate;
    typedef std::vector<TripletCoordinate> TripletCoordinates;


    typedef std::map<int, std::vector<NucleicAcidDB::NucleicAcidFull>> PossibleFragments;

public:
    MLFind(Predictions &prediction, clipper::Xmap<float> &xmap_wrk) : PredictionBuilder(prediction) {
        p_xmap = xmap_wrk;
    }

    /// @brief Run finding algorithm, dump phosphate positions by passing in path
    /// @param phosphate_position_path
    /// @return
    clipper::MiniMol find();

    /// @brief Load trinucleotide fragment into library
    /// @param path
    void load_library_from_file(const std::string &path);

    static clipper::MiniMol remove_bases(clipper::MiniMol &mol);

private:

    static MLFind::TripletCoordinates find_phosphate_triplets(clipper::MiniMol& mol);

    static float score_fragment(NucleicAcidDB::ChainFull &fragment, clipper::Xmap<float> &xmap);

    clipper::MiniMol filter_and_form_chain(PossibleFragments &fragments);

    static clipper::Coord_orth calculate_com(NucleicAcidDB::ChainFull& chain);

    NucleicAcidDB::ChainFull
    refine_fragment(NucleicAcidDB::ChainFull &original_fragment, float translation_range, float translation_step);

    static clipper::MiniMol
    organise_to_chains(clipper::MiniMol &mol);

    static void find_chains(int current_index, std::map<int, std::vector<int>>& connections, Chain& chain);

    NucleicAcidDB::ChainFull na_db;

};
#endif //NAUTILUS_NAUTILUS_MLFIND_H
