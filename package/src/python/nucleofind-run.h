#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>

#include <gemmi/fourier.hpp>  // for get_size_for_hkl, get_f_phi_on_grid, ...
#include "../cpp/nautilus-tools.h"
#include "../cpp/nautilus-target.h"
#include "../cpp/nautilus-join.h"
#include "../cpp/nautilus-sequence.h"
#include "../cpp/nautilus-rebuild-bases.h"
#include "../cpp/nautilus-tidy.h"
#include "../cpp/nautilus-util.h"
#include "../cpp/nautilus-mlfind.h"
#include "../cpp/nautilus-findml.h"
#include "../cpp/nucleofind/nucleofind.h"
#include "../cpp/nucleofind/predicted-maps.h"
#include "../cpp/nucleofind/input.h"

inline clipper::MiniMol &
run_cycle(const int verbose, const NucleicAcidTargets &natools, const clipper::MMoleculeSequence &seq_wrk,
          clipper::MiniMol &mol_wrk, const clipper::Xmap<float> &xwrk, NautilusLog &log, const bool strong_prune) {
    // grow chains

    int nas_found = NautilusUtil::count_nas(mol_wrk);
    if (nas_found == 0) { return mol_wrk;}


    ModelTidy::chain_renumber(mol_wrk, seq_wrk);
    NucleicAcidTools::chain_sort(mol_wrk);
    NucleicAcidTools::chain_label(mol_wrk, clipper::MMDBManager::CIF);
    NucleicAcidTools::residue_label(mol_wrk);

    mol_wrk = natools.grow(xwrk, mol_wrk, 25, 0.001);
    log.log("GROW", mol_wrk, verbose >= 5);

    // join
    if (!strong_prune) {
        NucleicAcidJoin na_join;
        mol_wrk = na_join.join(mol_wrk);
        log.log("JOIN", mol_wrk, verbose >= 5);
    }

    // link
    mol_wrk = natools.link(xwrk, mol_wrk);
    log.log("LINK", mol_wrk, verbose >= 5);

    // prune
    if (strong_prune) {
        mol_wrk = natools.strong_prune(mol_wrk);
    } else {
        mol_wrk = natools.prune(mol_wrk);
    }
    log.log("PRUNE", mol_wrk, verbose >= 5);

    // rebuild
    mol_wrk = natools.rebuild_chain(xwrk, mol_wrk);
    log.log("CHAIN", mol_wrk, verbose >= 5);

    // sequence
    NucleicAcidSequence na_seqnc;
    mol_wrk = na_seqnc.sequence(xwrk, mol_wrk, seq_wrk);
    log.log("SEQNC", mol_wrk, verbose >= 5);

    // rebuild
    NucleicAcidRebuildBases na_bases;
    mol_wrk = na_bases.rebuild_bases(xwrk, mol_wrk);
    log.log("BASES", mol_wrk, verbose >= 5);

    return mol_wrk;
}


void run(NucleoFind::IO::Input &input, NucleoFind::IO::Output &output, int cycles) {

    CCP4Program prog = NucleoFind::IO::initialise_ccp4_program("1.2");

    int nhit = 100;
    double res_in = 2.0;         // Resolution limit
    double srchst = 18.0;        // Search angle step
    int verbose = 0;
    clipper::String title;
    clipper::String ippdb_ref = "NONE";
    clipper::MMDBManager::TYPE cifflag = clipper::MMDBManager::Default;


    // Get reference model
    NautilusUtil::set_reference(ippdb_ref);
    NucleicAcidTargets natools;
    natools.add_pdb(ippdb_ref);
    NucleicAcidTools tools;

    // setup
    auto mtz_file = NucleoFind::IO::MTZ(input, res_in);
    auto prediction_files = NucleoFind::IO::Predictions(input);
    auto sequence_file = NucleoFind::IO::Sequence(input);
    clipper::Xmap<float> xwrk = mtz_file.get_xmap();
    auto model_file = NucleoFind::IO::Model(input, xwrk.spacegroup(), xwrk.cell());
    clipper::MiniMol mol_wrk = model_file.get_model();
    clipper::MMoleculeSequence seq_wrk = sequence_file.get_sequence();

    // store copy of input model
    clipper::MiniMol mol_wrk_in = mol_wrk;

    // map stats
    natools.init_stats(xwrk);
    NautilusLog log(title); // edited
    std::cout << std::endl;

    mol_wrk = NucleicAcidTools::flag_chains(mol_wrk);
    clipper::MiniMol mol_wrk_original = mol_wrk;

    // find with NucleoFind predictions
    NucleoFind::PredictedMaps predicted_maps = prediction_files.get_predictions();
    NucleoFind::Find find = {xwrk, predicted_maps};
    mol_wrk = find.find(mol_wrk, !(input.is_em()));
    log.log("FIND ML", mol_wrk, verbose >= 5);

    int nas_found = NautilusUtil::count_nas(mol_wrk);
    if (nas_found > 0) {
        mol_wrk = natools.link(xwrk, mol_wrk);
        log.log("FIND ML LINK", mol_wrk, verbose >= 5);

        ModelTidy::chain_renumber(mol_wrk, seq_wrk);
        NucleicAcidTools::chain_sort(mol_wrk);
        if (mol_wrk_in.size() > 0) {
            NucleicAcidTools::symm_match(mol_wrk, mol_wrk_in);
        }
        NucleicAcidTools::residue_label(mol_wrk);

        for (int cyc = 0; cyc < cycles; cyc++) {
            std::cout << "ML Based cycle " << clipper::String(cyc + 1, 3) << std::endl << std::endl;
            mol_wrk = run_cycle(verbose, natools, seq_wrk, mol_wrk, xwrk, log, true);
            NucleicAcidTools::chain_label(mol_wrk, clipper::MMDBManager::CIF);
            mol_wrk = NucleicAcidTools::chain_sort(mol_wrk);
            NucleicAcidTools::residue_label(mol_wrk);
        }
    }

    clipper::MiniMol best_model = mol_wrk;
    clipper::MiniMol ml_model = mol_wrk;

    int best_na_count = NautilusUtil::count_well_modelled_nas(mol_wrk, predicted_maps);
    std::cout << "Built " << best_na_count << " residues with positive predicted density" << std::endl;

    for (int cyc = 0; cyc < cycles; cyc++) {
        std::cout << "Internal cycle " << clipper::String(cyc + 1, 3) << std::endl << std::endl; // edited

        mol_wrk = NucleicAcidTools::flag_chains(mol_wrk);

        clipper::MiniMol nautilus_find = natools.find(xwrk, mol_wrk, nhit / 2, nhit / 2, srchst);
        mol_wrk = NucleoFind::Find::aggregate_nucleic_nucleic(ml_model, nautilus_find);
        log.log("FIND", mol_wrk, verbose >= 5);

        mol_wrk = run_cycle(verbose, natools, seq_wrk, mol_wrk, xwrk, log, false);

        int pred_na_count = NautilusUtil::count_well_modelled_nas(mol_wrk, predicted_maps);
        std::cout << "Cycle "<< cyc+1 << " built " << pred_na_count << " residues with positive predicted density" << std::endl;

        if (pred_na_count > best_na_count) {
            std::cout << "Taking model from old cycle " << cyc + 1 << "\n";
            best_model = mol_wrk;
            best_na_count = pred_na_count;
        }
    }

    std::cout << "Taking best model from all cycles with " << best_na_count << " nucleic acids residues with positive built." << std::endl;
    mol_wrk = best_model;

    // move to match input model
    if (mol_wrk_in.size() > 0) {
        NucleicAcidTools::symm_match(mol_wrk, mol_wrk_in);
        log.log("SYMMA", mol_wrk, verbose >= 5);
    }

    // set up residue types
    mol_wrk = ModelTidy::setup_residue_types(xwrk, mol_wrk);
    log.log("SETRES", mol_wrk, verbose >= 5);

    // add true sequence numbers
    ModelTidy::chain_renumber(mol_wrk, seq_wrk);
    log.log("TIDY", mol_wrk, verbose >= 5);

    //new chain labelling routine, for 2-char label, SWH
    NucleicAcidTools::chain_label(mol_wrk, cifflag);
    log.log("LABEL", mol_wrk, verbose >= 5);

    // final file output

    clipper::MMDBfile pdbfile;
    pdbfile.export_minimol(mol_wrk);
    pdbfile.write_file(output.get_pdb_out(), cifflag);

    if (output.get_xml_out().has_value()) {
        std::string xmlpath = output.get_xml_out().value();
        log.xml(xmlpath);
    }

    log.profile();
    prog.set_termination_message("Normal termination");

};

