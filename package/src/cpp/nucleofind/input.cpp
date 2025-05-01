//
// Created by Jordan Dialpuri on 01/05/2025.
//

#include "input.h"


CCP4Program NucleoFind::IO::initialise_ccp4_program(const std::string &version) {
    CCP4Program prog("nucleofind-build", version.c_str(), "$Date: 2024/02/26");
    prog.set_termination_message("Failed");

    std::cout << std::endl << "Copyright Jordan Dialpuri, Kathryn Cowtan, Jon Agirre, Paul Bond and the University of York." << std::endl
              << std::endl;
    prog.summary_beg();
    std::cout << "$TEXT:Reference: $$ Please reference $$" << std::endl << std::endl
              << " 'NucleoFind : a deep-learning network for interpreting nucleic acid electron density'" << std::endl
              << " J. S. Dialpuri, J. Agirre, K. D. Cowtan, and P. S. Bond, Nucleic Acids Research, 2024, 52, e84 " << std::endl
              << " https://doi.org/10.1093/nar/gkae715" << std::endl << std::endl << "$$"
              << std::endl;
    prog.summary_end();
    return prog;
}

void NucleoFind::IO::MTZ::load_work_map() {
    typedef clipper::HKL_data_base::HKL_reference_index HRI;
    using clipper::data32::Compute_fphi_from_fsigf_phifom;
    using clipper::data32::Compute_scale_u_aniso_fphi;
    clipper::Resolution resol;
    clipper::CCP4MTZfile mtzfile;
    mtzfile.set_column_label_mode(clipper::CCP4MTZfile::Legacy);

    // Get work reflection data
    mtzfile.open_read(input.get_mtz_path());
    double res = clipper::Util::max(mtzfile.resolution().limit(), res_in);
    resol = clipper::Resolution(res);
    hkls.init(mtzfile.spacegroup(), mtzfile.cell(), resol, true);
    clipper::HKL_data<clipper::data32::F_sigF> wrk_f(hkls);
    clipper::HKL_data<clipper::data32::ABCD> wrk_hl(hkls);
    clipper::HKL_data<clipper::data32::Phi_fom> wrk_pw(hkls);
    fphi = clipper::HKL_data<clipper::data32::F_phi>(hkls);
    clipper::HKL_data<clipper::data32::Flag> flag(hkls);
    if (input.get_fobs().has_value()) mtzfile.import_hkl_data(wrk_f, input.get_fobs().value());
    if (input.get_hl().has_value()) mtzfile.import_hkl_data(wrk_hl, input.get_hl().value());
    if (input.get_phifom().has_value()) mtzfile.import_hkl_data(wrk_pw, input.get_phifom().value());
    if (input.get_fc().has_value()) mtzfile.import_hkl_data(fphi, input.get_fc().value());
    if (input.get_free().has_value()) mtzfile.import_hkl_data(flag, input.get_free().value());
    mtzfile.close_read();

    // apply free flag
    wrk_f1 = wrk_f;
    //wrk_f1.mask( flag != 0 );
    for (HRI ih = hkls.first(); !ih.last(); ih.next())
        if (flag[ih].flag() == 0)
            wrk_f1[ih] = clipper::data32::F_sigF();  //ugly hack for broken SGI compilers
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
    xwrk = clipper::Xmap<float>(cspg, cxtl, grid);
    xwrk.fft_from(fphi);
}

void NucleoFind::IO::MTZ::print_stats() {
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
}

void NucleoFind::IO::Predictions::load_predicted_maps() {
    if (input.get_phosphate_prediction_path().has_value()) {
        clipper::CCP4MAPfile mapfile;
        mapfile.open_read(input.get_phosphate_prediction_path().value());
        mapfile.import_xmap(xphospred);
        mapfile.close_read();
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
}

void NucleoFind::IO::Sequence::load_sequence() {
    clipper::SEQfile seqf_wrk;
    seqf_wrk.read_file(input.get_seq_path());
    seqf_wrk.import_molecule_sequence(seq_wrk);
}

void NucleoFind::IO::Model::load_model() {
    if (input.get_pdb_path().has_value()) {
        constexpr int mmdbflags =
        ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors |
        ::mmdb::MMDBF_IgnoreRemarks;

        clipper::MiniMol mol_tmp;
        clipper::MMDBfile mmdb;
        mmdb.SetFlag(mmdbflags);
        mmdb.read_file(input.get_pdb_path().value());
        mmdb.import_minimol(mol_tmp);
        std::cout << mol_tmp.spacegroup().symbol_hm() << " " << mol_tmp.cell().format() << " "
                  << mol_tmp.atom_list().size() << std::endl;
        for (int c = 0; c < mol_tmp.size(); c++) mol_wrk.insert(mol_tmp[c]);
    }
}




