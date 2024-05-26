#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <algorithm>

#include "../cpp/nucleicacid_db.h"
#include "../cpp/nautilus-tools.h"
#include "../cpp/nautilus-ss-find.h"
#include "../cpp/nautilus-target.h"
#include "../cpp/nautilus-join.h"
#include "../cpp/nautilus-sequence.h"
#include "../cpp/nautilus-rebuild-bases.h"
#include "../cpp/nautilus-tidy.h"
#include "../cpp/nautilus-util.h"
#include "../cpp/nautilus-mlfind.h"
#include "../cpp/nautilus-findml.h"

clipper::MiniMol &
run_cycle(int nhit, double srchst, int verbose, NucleicAcidTargets &natools, const clipper::MMoleculeSequence &seq_wrk,
          clipper::MiniMol &mol_wrk, const clipper::Xmap<float> &xwrk, NautilusLog &log, PredictedMaps& predictions) {
    // grow chains

    int nas_found = NautilusUtil::count_nas(mol_wrk);
    if (nas_found == 0) { return mol_wrk;}

    mol_wrk = natools.grow(xwrk, mol_wrk, 25, 0.001, predictions);
    log.log("GROW", mol_wrk, verbose >= 5);

//    NautilusUtil::save_minimol(mol_wrk, "grow.pdb");

    // join
    NucleicAcidJoin na_join;
    mol_wrk = na_join.join(mol_wrk);
    log.log("JOIN", mol_wrk, verbose >= 5);

//    NautilusUtil::save_minimol(mol_wrk, "join.pdb");


    // link
    mol_wrk = natools.link(xwrk, mol_wrk);
    log.log("LINK", mol_wrk, verbose >= 5);

//    NautilusUtil::save_minimol(mol_wrk, "link.pdb");

    // prune
    mol_wrk = natools.prune(mol_wrk);
    log.log("PRUNE", mol_wrk, verbose >= 5);

//    NautilusUtil::save_minimol(mol_wrk, "prune.pdb");


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

void run(NautilusInput &input, NautilusOutput &output, int cycles) {

    CCP4Program prog("nucleofind-build", "0.6.1", "$Date: 2024/02/26");
    prog.set_termination_message("Failed");

    std::cout << std::endl << "Copyright (c) 2024 Jordan Dialpuri, Jon Agirre, Kevin Cowtan, Paul Bond and University of York. All rights reserved" << std::endl
              << std::endl;
    prog.summary_beg();
    std::cout << "$TEXT:Reference: $$ Please reference $$" << std::endl << std::endl
              << " 'Automated nucleic acid chain tracing in real time'" << std::endl
              << " Cowtan K. (2014). IUCrJ 1, doi:10.1107/S2052252514019290. " << std::endl << std::endl << "$$"
              << std::endl;
    prog.summary_end();

    bool doanis = false;
    int nhit = 100;
    double res_in = 2.0;         // Resolution limit
    double srchst = 18.0;        // Search angle step
    int verbose = 0;
    clipper::String title;
    clipper::String ippdb_ref = "NONE";
    clipper::MMDBManager::TYPE cifflag = clipper::MMDBManager::Default;

    typedef clipper::HKL_data_base::HKL_reference_index HRI;
    using clipper::data32::Compute_fphi_from_fsigf_phifom;
    using clipper::data32::Compute_scale_u_aniso_fphi;
    clipper::Resolution resol;
    clipper::CCP4MTZfile mtzfile;
    mtzfile.set_column_label_mode(clipper::CCP4MTZfile::Legacy);
    const int mmdbflags =
            ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors |
            ::mmdb::MMDBF_IgnoreRemarks;
    typedef clipper::HKL_data_base::HKL_reference_index HRI;

    // Get work reflection data
    clipper::HKL_info hkls;
    mtzfile.open_read(input.get_mtz_path());
    double res = clipper::Util::max(mtzfile.resolution().limit(), res_in);
    resol = clipper::Resolution(res);
    hkls.init(mtzfile.spacegroup(), mtzfile.cell(), resol, true);
    clipper::HKL_data<clipper::data32::F_sigF> wrk_f(hkls);
    clipper::HKL_data<clipper::data32::ABCD> wrk_hl(hkls);
    clipper::HKL_data<clipper::data32::Phi_fom> wrk_pw(hkls);
    clipper::HKL_data<clipper::data32::F_phi> fphi(hkls);
    clipper::HKL_data<clipper::data32::Flag> flag(hkls);
    if (input.get_fobs().has_value()) mtzfile.import_hkl_data(wrk_f, input.get_fobs().value());
    if (input.get_hl().has_value()) mtzfile.import_hkl_data(wrk_hl, input.get_hl().value());
    if (input.get_phifom().has_value()) mtzfile.import_hkl_data(wrk_pw, input.get_phifom().value());
    if (input.get_fc().has_value()) mtzfile.import_hkl_data(fphi, input.get_fc().value());
    if (input.get_free().has_value()) mtzfile.import_hkl_data(flag, input.get_free().value());
    mtzfile.close_read();

    // do anisotropy correction
    clipper::U_aniso_orth uaniso(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    if (input.get_fobs().has_value()) for (HRI ih = hkls.first(); !ih.last(); ih.next()) wrk_f[ih].f() = fphi[ih].f();
    // scale obs data
    typedef clipper::SFscale_aniso<float> SFscale;
    SFscale sfscl(3.0, SFscale::SHARPEN);
    sfscl(wrk_f);
    uaniso = sfscl.u_aniso_orth(SFscale::F);
    // scale map coeffs
    Compute_scale_u_aniso_fphi compute_aniso(1.0, -uaniso);
    if (input.get_fc().has_value()) fphi.compute(fphi, compute_aniso);
    // output
    std::cout << std::endl << "Applying anisotropy correction:" << std::endl << uaniso.format() << std::endl
              << std::endl;


    // apply free flag
    clipper::HKL_data<clipper::data32::F_sigF> wrk_f1 = wrk_f;
    //wrk_f1.mask( flag != 0 );
    for (HRI ih = hkls.first(); !ih.last(); ih.next())
        if (flag[ih].flag() == 0)
            wrk_f1[ih] = clipper::data32::F_sigF();  //ugly hack for broken SGI compilers

    // Get reference model
    NautilusUtil::set_reference(ippdb_ref);
    NucleicAcidTargets natools;
    natools.add_pdb(ippdb_ref);
    NucleicAcidTools tools;

    // Get sequence
    clipper::MMoleculeSequence seq_wrk;
    clipper::SEQfile seqf_wrk;
    seqf_wrk.read_file(input.get_seq_path());
    seqf_wrk.import_molecule_sequence(seq_wrk);


    // Get model
    clipper::MiniMol mol_wrk(hkls.spacegroup(), hkls.cell());
    if (input.get_pdb_path().has_value()) {
        clipper::MiniMol mol_tmp;
        clipper::MMDBfile mmdb;
        mmdb.SetFlag(mmdbflags);
        mmdb.read_file(input.get_pdb_path().value());
        mmdb.import_minimol(mol_tmp);
        std::cout << mol_tmp.spacegroup().symbol_hm() << " " << mol_tmp.cell().format() << " "
                  << mol_tmp.atom_list().size() << std::endl;
        for (int c = 0; c < mol_tmp.size(); c++) mol_wrk.insert(mol_tmp[c]);
    }


    // work map
    if (!input.get_hl().has_value())
        wrk_hl.compute(wrk_pw, clipper::data32::Compute_abcd_from_phifom());
    if (!input.get_phifom().has_value())
        wrk_pw.compute(wrk_hl, clipper::data32::Compute_phifom_from_abcd());
    if (!input.get_fc().has_value())
        fphi.compute(wrk_f1, wrk_pw, Compute_fphi_from_fsigf_phifom());

    clipper::Spacegroup cspg = hkls.spacegroup();
    clipper::Cell cxtl = hkls.cell();
    clipper::Grid_sampling grid(cspg, cxtl, hkls.resolution());
    clipper::Xmap<float> xwrk(cspg, cxtl, grid);
    xwrk.fft_from(fphi);

    // output some statistics
    std::cout << std::endl;
    std::cout << " Spgr " << hkls.spacegroup().symbol_xhm() << std::endl;
    std::cout << hkls.cell().format() << std::endl;
    std::cout << " Nref " << hkls.num_reflections() << " " << fphi.num_obs() << std::endl;
    double smax = 0.0;
    for (HRI ih = fphi.first(); !ih.last(); ih.next())
        if (!fphi[ih].missing())
            if (fphi[ih].f() > 0.0)
                smax = std::max(smax, double(hkls.invresolsq(ih.index())));
    std::cout << " Reso " << hkls.resolution().limit() << " " << 1.0 / std::max(sqrt(smax), 1.0e-3) << std::endl;
    if (input.get_fobs().has_value()) {
        double sf(0.0), sw(0.0);
        for (HRI ih = fphi.first(); !ih.last(); ih.next())
            if (!wrk_f1[ih].missing() && !fphi[ih].missing()) {
                sf += wrk_f1[ih].f();
                sw += fphi[ih].f();
            }
        std::cout << " Fw/Fo " << sw / sf << std::endl;
    }
    if (smax == 0.0) { std::cerr << "No density provided." << std::endl; }

    // store copy of input model
    clipper::MiniMol mol_wrk_in = mol_wrk;

    // map stats
    natools.init_stats(xwrk);
    NautilusLog log(title); // edited
    std::cout << std::endl;

    // initial phosphates from PDB
    // if ( mol_pho.size() > 0 ) mol_wrk = natools.phosphate( xwrk, mol_wrk, mol_pho );

    clipper::Xmap<float> xphospred;
    clipper::Xmap<float> xsugarpred;
    clipper::Xmap<float> xbasepred;

    if (input.get_phosphate_prediction_path().has_value()) {
        clipper::CCP4MAPfile mapfile;
        mapfile.open_read(input.get_phosphate_prediction_path().value());
        mapfile.import_xmap(xphospred);
    }

    if (input.get_sugar_prediction_path().has_value()) {
        clipper::CCP4MAPfile mapfile;
        mapfile.open_read(input.get_sugar_prediction_path().value());
        mapfile.import_xmap(xsugarpred);
        mapfile.close_read();
    }

    if (input.get_base_prediction_path().has_value()) {
        clipper::CCP4MAPfile mapfile;
        mapfile.open_read(input.get_base_prediction_path().value());
        mapfile.import_xmap(xbasepred);
        mapfile.close_read();
    }

    PredictedMaps predictions = {xphospred, xsugarpred, xbasepred};

    mol_wrk = NucleicAcidTools::flag_chains(mol_wrk);
    clipper::MiniMol mol_wrk_original = mol_wrk;

    FindML find_ml = FindML(mol_wrk, xwrk, predictions);
    find_ml.load_library_from_file(ippdb_ref);
    find_ml.set_resolution(hkls.resolution().limit()); // Needed for the RSRZ calculation, but if not set defaults to 2
    mol_wrk = find_ml.find();
    log.log("FIND ML", mol_wrk, verbose >= 5);
    mol_wrk = natools.link(xwrk, mol_wrk);
    log.log("FIND ML LINK", mol_wrk, verbose >= 5);

    ModelTidy::chain_renumber(mol_wrk, seq_wrk);
    NucleicAcidTools::chain_sort(mol_wrk);

    // NautilusUtil::save_minimol(mol_wrk, "find.pdb");

    int nas_found = NautilusUtil::count_nas(mol_wrk);

    if (nas_found > 0) {
        for (int cyc = 0; cyc < cycles; cyc++) {
            std::cout << "ML Based cycle " << clipper::String(cyc + 1, 3) << std::endl << std::endl;
            mol_wrk = run_cycle(nhit, srchst, verbose, natools, seq_wrk, mol_wrk, xwrk, log, predictions);
            NucleicAcidTools::chain_label(mol_wrk, clipper::MMDBManager::CIF);
            mol_wrk = NucleicAcidTools::chain_sort(mol_wrk);
            NucleicAcidTools::residue_label(mol_wrk);
        }
    }

//    NautilusUtil::save_minimol(mol_wrk, "mlbuiltmodel.pdb");
    clipper::MiniMol best_model = mol_wrk;

    int best_na_count = NautilusUtil::count_well_modelled_nas(mol_wrk, xwrk, hkls.resolution().limit());

    std::cout << "Built " << best_na_count << " residues with RSRZ >= -1" << std::endl;

    for (int cyc = 0; cyc < cycles; cyc++) {
        std::cout << "Internal cycle " << clipper::String(cyc + 1, 3) << std::endl << std::endl; // edited

        mol_wrk = NucleicAcidTools::flag_chains(mol_wrk);

        mol_wrk = natools.find(xwrk, mol_wrk, nhit / 2, nhit / 2, srchst);
        log.log("FIND", mol_wrk, verbose >= 5);

        mol_wrk = run_cycle(nhit, srchst, verbose, natools, seq_wrk, mol_wrk, xwrk, log, predictions);

        int current_count = NautilusUtil::count_well_modelled_nas(mol_wrk, xwrk, hkls.resolution().limit());
        std::cout << "Cycle "<< cyc+1 << " built " << best_na_count << " residues with RSRZ >= -1" << std::endl;

        if (current_count > best_na_count) {
            std::cout << "Taking model from old cycle " << cyc + 1 << "\n";
            best_model = mol_wrk;
            best_na_count = current_count;
        }
    }

    std::cout << "Taking best model from all cycles with " << best_na_count << " nucleic acids residues with RSRZ >= -1 built." << std::endl;
    mol_wrk = best_model;

    // move to match input model
    if (mol_wrk_in.size() > 0) {
        NucleicAcidTools::symm_match(mol_wrk, mol_wrk_in);
        log.log("SYMMA", mol_wrk, verbose >= 5);
    }

    // set up residue types
    const clipper::String basetypes = "ACGTU";
    clipper::MiniMol mol_new(xwrk.spacegroup(), xwrk.cell());
    for (int c = 0; c < mol_wrk.size(); c++) {
        clipper::MPolymer mpx = mol_wrk[c];
        if (!mpx.exists_property("NON-NA")) {
            bool dna = false;
            for (int r = 0; r < mpx.size(); r++) {
                const clipper::String type = mpx[r].type().trim() + " ";
                const char ctype = type[0];
                int t = NucleicAcidTools::base_index(ctype);
                if (t >= 0) mpx[r].set_type("  " + basetypes.substr(t, 1));
                else mpx[r].set_type("  U");
                if (mpx[r].type().trim() == "T") dna = true;
            }
            // for DNA, prefix type with D and remove O2'
            if (dna) {
                for (int r = 0; r < mpx.size(); r++) {
                    const clipper::String type = " D" + mpx[r].type().trim();
                    mpx[r].set_type(type);
                    const int io2 = mpx[r].lookup(" O2'", clipper::MM::ANY);
                    clipper::MMonomer mm;
                    mm.copy(mpx[r], clipper::MM::COPY_MP);
                    for (int a = 0; a < mpx[r].size(); a++) {
                        if (a != io2) mm.insert(mpx[r][a]);
                    }
                    mpx[r] = mm;
                }
            }
        }
        mol_new.insert(mpx);
    }
    log.log("SETRES", mol_wrk, verbose >= 5);

    // add true sequence numbers
    ModelTidy::chain_renumber(mol_new, seq_wrk);
    log.log("TIDY", mol_new, verbose >= 5);
    //new chain labelling routine, for 2-char label, SWH
    NucleicAcidTools::chain_label(mol_new, cifflag);
    log.log("LABEL", mol_new, verbose >= 5);
    // final file output
    clipper::MMDBfile pdbfile;
    pdbfile.export_minimol(mol_new);
    pdbfile.write_file(output.get_pdb_out(), cifflag);
    if (output.get_xml_out().has_value()) { 
        std::string xmlpath = output.get_xml_out().value();
        log.xml(xmlpath);
    }
//   msg = log.log_info( mol_new, true );	// added by SWH
//   std::cout << "$TEXT:Result: $$ $$" << std::endl << msg << "\n$$" << std::endl; // added by SWH
    log.profile();
    prog.set_termination_message("Normal termination");

};