//
// Created by Jordan Dialpuri on 06/03/2024.
//

#ifndef NAUTILUS_FINDML_H
#define NAUTILUS_FINDML_H
#include <clipper/clipper-minimol.h>
#include <optional>
#include <unordered_set>
#include <set>
#include "nucleicacid_db.h"
//
struct ChainData {
    std::vector<int> ordered_chain;
    std::unordered_set<int> lookup_list;
};
struct PlacedFragmentResult;

class PredictedMaps {
public:
    PredictedMaps(const clipper::Xmap<float> &phosphate_map, const clipper::Xmap<float> &sugar_map,
                  const clipper::Xmap<float> &base_map)
            : phosphate(phosphate_map), sugar(sugar_map), base(base_map) {}

    [[nodiscard]] std::optional<clipper::Xmap<float>> get_phosphate_map() const {
        if (!phosphate.cell().is_null()) {return phosphate;}
        return std::nullopt;
    }

    [[nodiscard]] std::optional<clipper::Xmap<float>> get_sugar_map() const {
        if (!sugar.cell().is_null()) {return sugar;}
        return std::nullopt;
    }

    [[nodiscard]] const clipper::Xmap<float>* get_sugar_map_ptr() const {
        if (!sugar.cell().is_null()) {return &sugar;}
        return nullptr;
    }

    [[nodiscard]] const clipper::Xmap<float>* get_base_map_ptr() const {
        if (!base.cell().is_null()) {return &base;}
        return nullptr;
    }


    [[nodiscard]] std::optional<clipper::Xmap<float>> get_base_map() const {
        if (!base.cell().is_null()) {return base;}
        return std::nullopt;
    }

private:
    clipper::Xmap<float> phosphate;
    clipper::Xmap<float> sugar;
    clipper::Xmap<float> base;
};

class FindML {
public:
    explicit FindML(const clipper::MiniMol &mol, const clipper::Xmap<float> &xwrk, PredictedMaps predictions);

    void load_library_from_file(const std::string& path) {nadb.add_pdb(path);}

    clipper::MiniMol find();

    void set_resolution(double resolution) { m_resolution = resolution;}

public:
    typedef std::vector<std::pair<int, clipper::Coord_orth>> TripletCoordinate;
    typedef std::vector<TripletCoordinate> TripletCoordinates;
    typedef std::map<std::pair<int, int>, std::vector<NucleicAcidDB::NucleicAcidFull>> PossibleFragments;
    typedef std::vector<std::pair<std::vector<int>, std::vector<int>>> PairedChainIndices;
private:

    /*
     * PREDICTIONS TO POINTS START
     */

    clipper::MiniMol
    generate_molecule_from_gridpoints(clipper::Xmap<float> &predicted_map, double value_threshold);

    clipper::MiniMol calculate_phosphate_peaks(double value_threshold);

    clipper::MiniMol calculate_sugar_peaks(double value_threshold);

    clipper::MiniMol calculate_base_peaks(double value_threshold);


    static clipper::Coord_grid ascend_grid_gradient(const clipper::Coord_grid &grid_point, const clipper::Xmap<float> &xmap);

    [[nodiscard]] clipper::MiniMol
    find_peaks(const clipper::Xmap<float> &predicted_map, const clipper::MiniMol &mol) const;

    [[nodiscard]] clipper::MiniMol assimilate_peaks(clipper::MiniMol& peaks, float radius, const std::string& name) const;

    [[nodiscard]] clipper::MiniMol refine_phosphate_peaks(const clipper::MiniMol& phosphate_peaks) const;

    [[nodiscard]] TripletCoordinates symmetrise_phosphate_peaks(TripletCoordinates&triplet_coordinates, clipper::MiniMol&phosphate_peaks) const;

    /*
     * PREDICTIONS TO POINTS END
     */

    /*
     * PREDICTED POINTS TO MOLECULE START
     */

    static FindML::TripletCoordinates
    find_triplet_coordinates(const clipper::MiniMol &phosphate_peaks, const clipper::MiniMol &sugar_peaks);

    NucleicAcidDB::ChainFull
    refine_fragment(NucleicAcidDB::ChainFull &original_fragment, float translation_range, float translation_step);

    NucleicAcidDB::ChainFull refine_fragment_coordinates(NucleicAcidDB::ChainFull& original_fragment);

    static clipper::Coord_orth calculate_com(NucleicAcidDB::ChainFull &chain);

    float score_fragment(NucleicAcidDB::ChainFull &fragment, clipper::Xmap<float> &xmap);

    static float score_density(NucleicAcidDB::NucleicAcidFull &chain, clipper::Xmap<float> &xmap, bool terminal);

    float score_sugar(NucleicAcidDB::NucleicAcidFull &chain);

    float score_base(NucleicAcidDB::NucleicAcidFull &chain);

    [[ nodiscard ]] clipper::MiniMol remove_bases(clipper::MiniMol &mol);

    [[ nodiscard ]] clipper::MiniMol filter_and_form_single_chain(PossibleFragments &fragments) const;

    [[ nodiscard ]] clipper::MiniMol filter_and_form_bidirectional_chain(PossibleFragments& fragments) const;

    [[ nodiscard ]] clipper::MiniMol form_chain(PossibleFragments& fragments) const;

    [[ nodiscard ]] clipper::MiniMol form_organised_chains(PossibleFragments& fragments, std::vector<std::vector<int>>& fragment_indices) const;

    [[ nodiscard ]] clipper::MiniMol organise_to_chains(clipper::MiniMol &mol);

    [[ nodiscard ]] clipper::MiniMol remove_clashing_protein(clipper::MiniMol& na_chain);

    [[ nodiscard ]] clipper::MiniMol remove_low_confidence(clipper::MiniMol& mol);

    PlacedFragmentResult place_fragments(const clipper::MiniMol& phosphate_peaks, const std::vector<int>& positions);

    static void find_chains(int current_index, std::map<int, std::vector<int>>& connections, ChainData& chain);

    [[ nodiscard ]] PairedChainIndices organise_triplets_to_chains(TripletCoordinates&triplets);

    /*
     * PREDICTED POINTS TO MOLECULE END
     */

    /*
     * DEBUG -> TO BE REMOVED
     */

    void draw_triplets(TripletCoordinates &triplets, clipper::MiniMol &phosphate_positions, const std::string& path) {

        clipper::MiniMol drawn_mol = {phosphate_positions.spacegroup(), phosphate_positions.cell()};

        clipper::MPolymer mp;
        mp.set_id("D");
        for (int i = 0; i < triplets.size(); i++) {
            clipper::MMonomer monomer;
            monomer.set_type("P");
            monomer.set_id(i);
            clipper::Coord_orth p1 = triplets[i][0].second;
            clipper::Coord_orth p2 = triplets[i][1].second;
            clipper::Coord_orth p3 = triplets[i][2].second;

            clipper::Vec3<> p1p2 = p1 - p2;
            clipper::Vec3<> p2p3 = p2 - p3;

            clipper::Vec3<double> p1p2step = {p1p2[0] / 10, p1p2[1] / 10, p1p2[2] / 10};
            clipper::Vec3<double> p2p3step = {p2p3[0] / 10, p2p3[1] / 10, p2p3[2] / 10};

            //        std::cout << p1.format() << " " << p2.format() << " " << p1p2.format() << " " << p1p2step.format() << std::endl;

            for (int j = 1; j <= 10; j++) {

                clipper::Vec3<> p1p2_stepped = p2 + (double(j) * p1p2step);
                clipper::Vec3<> p2p3_stepped = p3 + (double(j) * p2p3step);

                clipper::Coord_orth p1p2_stepped_orth = clipper::Coord_orth(p1p2_stepped);
                clipper::Coord_orth p2p3_stepped_orth = clipper::Coord_orth(p2p3_stepped);

                clipper::MAtom atom_p1 = NautilusUtil::create_atom(p1p2_stepped_orth, i + j, "P");
                monomer.insert(atom_p1);

                clipper::MAtom atom_p2 = NautilusUtil::create_atom(p2p3_stepped_orth, i + j + 3, "P");
                monomer.insert(atom_p2);
            }
            mp.insert(monomer);
        }
        drawn_mol.model().insert(mp);
        NautilusUtil::save_minimol(drawn_mol, path);

    }

private:
    clipper::Xmap<float> xwrk;
    PredictedMaps predictions;
    clipper::Xmap<float> xphospred;
    clipper::Xmap<float> xsugarpred;
    clipper::Xmap<float> xbasepred;

    clipper::MiniMol mol;
    NucleicAcidDB::ChainFull nadb;

    double m_resolution = 2;

};

struct PlacedFragmentResult {
    float score;
    FindML::PossibleFragments fragments;

    PlacedFragmentResult(const float score, const FindML::PossibleFragments& fragments) { this->score = score; this->fragments = fragments;}
};


#endif //NAUTILUS_FINDML_H
