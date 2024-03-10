// Clipper ssfind
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */

//#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <algorithm>

#include "nucleicacid_db.h"
#include "nautilus-tools.h"
#include "nautilus-ss-find.h"
#include "nautilus-target.h"
#include "nautilus-join.h"
#include "nautilus-sequence.h"
#include "nautilus-rebuild-bases.h"
#include "nautilus-tidy.h"
#include "nautilus-util.h"
// #include "nautilus-mlfind.h"



int main( int argc, char** argv )
{
  CCP4Program prog( "cnautilus", "0.5.4", "$Date: 2021/09/21" );
  prog.set_termination_message( "Failed" );

  std::cout << std::endl << "Copyright 2011-2017 Kevin Cowtan and University of York." << std::endl << std::endl;
  prog.summary_beg();
  std::cout << "$TEXT:Reference: $$ Please reference $$" << std::endl << std::endl << " 'Automated nucleic acid chain tracing in real time'" << std::endl << " Cowtan K. (2014). IUCrJ 1, doi:10.1107/S2052252514019290. " << std::endl << std::endl << "$$" << std::endl;
  prog.summary_end();


  // defaults
  clipper::String title;
  clipper::String ipmtz = "NONE";
  clipper::String ipcol_fo = "NONE";
  clipper::String ipcol_hl = "NONE";
  clipper::String ipcol_pw = "NONE";
  clipper::String ipcol_fc = "NONE";
  clipper::String ipcol_fr = "NONE";
  clipper::String ipseq = "NONE";
  clipper::String ippdb = "NONE";
  clipper::String ippdb_ref = "NONE";
  clipper::String ippdb_pho = "NONE";
  clipper::String ippredicted_phos_map = "NONE"; // added by JSD

  clipper::String oppdb = "nautilus.pdb";
  clipper::String opmap = "NONE";
  clipper::String opxml = "NONE";
  clipper::String msg; // added by SWH
  clipper::MMDBManager::TYPE cifflag = clipper::MMDBManager::Default;
  int ncyc = 3;
  bool doanis = false;
  int nhit = 100;
  double res_in = 2.0;         // Resolution limit
  double srchst = 18.0;        // Search angle step
  int verbose = 0;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if        ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-pdbin-ref" ) {
      if ( ++arg < args.size() ) ippdb_ref = args[arg];
    } else if ( args[arg] == "-pdbin-phosphate" ) {
      if ( ++arg < args.size() ) ippdb_pho = args[arg];
    } else if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipmtz = args[arg];
    } else if ( args[arg] == "-seqin" ) {
      if ( ++arg < args.size() ) ipseq = args[arg];
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) ippdb = args[arg];
    } else if ( args[arg] == "-pdbout" ) {
      if ( ++arg < args.size() ) oppdb = args[arg];
    } else if (args[arg] == "-predicted-phos-map") {
        if (++arg < args.size()) ippredicted_phos_map = args[arg];
    } else if ( args[arg] == "-mapout" ) {
      if ( ++arg < args.size() ) opmap  = args[arg];
    } else if ( args[arg] == "-xmlout" ) {
      if ( ++arg < args.size() ) opxml  = args[arg];
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcol_fo = args[arg];
    } else if ( args[arg] == "-colin-hl" ) {
      if ( ++arg < args.size() ) ipcol_hl = args[arg];
    } else if ( args[arg] == "-colin-phifom" ) {
      if ( ++arg < args.size() ) ipcol_pw = args[arg];
    } else if ( args[arg] == "-colin-fc" ) {
      if ( ++arg < args.size() ) ipcol_fc = args[arg];
    } else if ( args[arg] == "-colin-free" ) {
      if ( ++arg < args.size() ) ipcol_fr = args[arg];
    } else if ( args[arg] == "-cycles" ) {
      if ( ++arg < args.size() ) ncyc  = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-anisotropy-correction" ) {
      doanis = true;
    } else if ( args[arg] == "-fragments" ) {
      if ( ++arg < args.size() ) nhit = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) res_in = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-cif" ){      // for cif output
      cifflag = clipper::MMDBManager::CIF;
    } else if ( args[arg] == "-search-step" ) {
      if ( ++arg < args.size() ) srchst = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "\nUsage: cnautilus\n\t-mtzin <filename>\t\tCOMPULSORY\n\t-seqin <filename>\n\t-pdbin <filename>\n\t-pdbout <filename>\n\t-predicted-phos-map <filename>\n\t-xmlout <filename>\n\t-colin-fo <colpath>\n\t-colin-hl <colpath> or -colin-phifom <colpath>\n\t-colin-fc <colpath>\n\t-colin-free <colpath>\n\t-cycles <number>\n\t-anisotropy-correction\n\t-fragments <number>\n\t-resolution <resolution/A>\n\t-pdbin-ref <filename>\n\t-cif\t\t*will only output model in cif format\n\t-search-step <float>\t*search angle step\n\t-verbose <number>\n";
    std::cout << "\nAn input mtz file for the work structure is required. Chains will be located and \ngrown for the work structure and written to the output pdb/cif file.\n";
    std::cout << "This involves following steps:\n finding, growing, joining, linking, pruning, rebuilding chains, sequencing, and rebuilding bases.\nIf the optional input pdb file is provided for the work structure, then the input model is extended.\n";
    return 1;
  }

  // check data present
  if ( ipcol_fc == "NONE" && ipcol_fo == "NONE" )
    { std::cerr << "No F's provided." << std::endl; return 1; }
  if ( ipcol_fc == "NONE" && ipcol_hl == "NONE" && ipcol_pw == "NONE" )
    { std::cerr << "No phases provided." << std::endl; return 1; }

  // other initialisations
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  using clipper::data32::Compute_fphi_from_fsigf_phifom;
  using clipper::data32::Compute_scale_u_aniso_fphi;
  clipper::Resolution resol;
  clipper::CCP4MTZfile mtzfile;
  mtzfile.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
  const int mmdbflags = ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors | ::mmdb::MMDBF_IgnoreRemarks;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;

  // Get work reflection data
  clipper::HKL_info hkls;
  mtzfile.open_read( ipmtz );
  double res = clipper::Util::max( mtzfile.resolution().limit(), res_in );
  resol = clipper::Resolution( res );
  hkls.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  clipper::HKL_data<clipper::data32::F_sigF>  wrk_f ( hkls );
  clipper::HKL_data<clipper::data32::ABCD>    wrk_hl( hkls );
  clipper::HKL_data<clipper::data32::Phi_fom> wrk_pw( hkls );
  clipper::HKL_data<clipper::data32::F_phi>   fphi( hkls );
  clipper::HKL_data<clipper::data32::Flag>    flag( hkls );
  if ( ipcol_fo != "NONE" ) mtzfile.import_hkl_data( wrk_f ,ipcol_fo );
  if ( ipcol_hl != "NONE" ) mtzfile.import_hkl_data( wrk_hl,ipcol_hl );
  if ( ipcol_pw != "NONE" ) mtzfile.import_hkl_data( wrk_pw,ipcol_pw );
  if ( ipcol_fc != "NONE" ) mtzfile.import_hkl_data( fphi,  ipcol_fc );
  if ( ipcol_fr != "NONE" ) mtzfile.import_hkl_data( flag,  ipcol_fr );
  mtzfile.close_read();

  // do anisotropy correction
  clipper::U_aniso_orth uaniso( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
  if ( doanis ) {
    if ( ipcol_fo == "NONE" ) for ( HRI ih = hkls.first(); !ih.last(); ih.next() ) wrk_f[ih].f() = fphi[ih].f();
    // scale obs data
    typedef clipper::SFscale_aniso<float> SFscale;
    SFscale sfscl( 3.0, SFscale::SHARPEN );
    sfscl( wrk_f );
    uaniso = sfscl.u_aniso_orth( SFscale::F );
    // scale map coeffs
    Compute_scale_u_aniso_fphi compute_aniso( 1.0, -uaniso );
    if ( ipcol_fc != "NONE" ) fphi.compute( fphi, compute_aniso );
    // output
    std::cout << std::endl << "Applying anisotropy correction:"
              << std::endl << uaniso.format() << std::endl << std::endl;
  }

  // apply free flag
  clipper::HKL_data<clipper::data32::F_sigF> wrk_f1 = wrk_f;
  //wrk_f1.mask( flag != 0 );
  for ( HRI ih = hkls.first(); !ih.last(); ih.next() ) if ( flag[ih].flag() == 0 ) wrk_f1[ih] = clipper::data32::F_sigF();  //ugly hack for broken SGI compilers

  // Get reference model
  if ( ippdb_ref == "NONE" ) NautilusUtil::set_reference( ippdb_ref );
  NucleicAcidTargets natools;
  natools.add_pdb( ippdb_ref );
  NucleicAcidTools tools;

  // Get sequence
  clipper::MMoleculeSequence seq_wrk;
  if ( ipseq != "NONE" ) {
    clipper::SEQfile seqf_wrk;
    seqf_wrk.read_file( ipseq );
    seqf_wrk.import_molecule_sequence( seq_wrk );
  }

  // Get model
  clipper::MiniMol mol_wrk( hkls.spacegroup(), hkls.cell() );
  if ( ippdb != "NONE" ) {
    clipper::MiniMol mol_tmp;
    clipper::MMDBfile mmdb;
    mmdb.SetFlag( mmdbflags );
    mmdb.read_file( ippdb );
    mmdb.import_minimol( mol_tmp );
    std::cout << mol_tmp.spacegroup().symbol_hm() << " " << mol_tmp.cell().format() << " " << mol_tmp.atom_list().size() << std::endl;
    for ( int c = 0; c < mol_tmp.size(); c++ ) mol_wrk.insert( mol_tmp[c] );
  }

  // Get model phosphates
  clipper::MiniMol mol_pho( hkls.spacegroup(), hkls.cell() );
  if ( ippdb_pho != "NONE" ) {
    clipper::MiniMol mol_tmp;
    clipper::MMDBfile mmdb;
    mmdb.SetFlag( mmdbflags );
    mmdb.read_file( ippdb_pho );
    mmdb.import_minimol( mol_tmp );
    std::cout << mol_tmp.spacegroup().symbol_hm() << " " << mol_tmp.cell().format() << " " << mol_tmp.atom_list().size() << std::endl;
    for ( int c = 0; c < mol_tmp.size(); c++ ) mol_pho.insert( mol_tmp[c] );
  }

  // work map
  if ( ipcol_hl == "NONE" )
    wrk_hl.compute( wrk_pw, clipper::data32::Compute_abcd_from_phifom() );
  if ( ipcol_pw == "NONE" )
    wrk_pw.compute( wrk_hl, clipper::data32::Compute_phifom_from_abcd() );
  if ( ipcol_fc == "NONE" )
    fphi.compute( wrk_f1, wrk_pw, Compute_fphi_from_fsigf_phifom() );
  clipper::Spacegroup cspg = hkls.spacegroup();
  clipper::Cell       cxtl = hkls.cell();
  clipper::Grid_sampling grid( cspg, cxtl, hkls.resolution() );
  clipper::Xmap<float>   xwrk( cspg, cxtl, grid );
  xwrk.fft_from( fphi );

  // output some statistics
  std::cout << std::endl;
  std::cout << " Spgr " << hkls.spacegroup().symbol_xhm() << std::endl;
  std::cout << hkls.cell().format() << std::endl;
  std::cout << " Nref " << hkls.num_reflections() << " " << fphi.num_obs() << std::endl;
  double smax = 0.0;
  for ( HRI ih = fphi.first(); !ih.last(); ih.next() )
    if ( !fphi[ih].missing() )
      if ( fphi[ih].f() > 0.0 )
        smax = std::max( smax, double(hkls.invresolsq(ih.index())) );
  std::cout << " Reso " << hkls.resolution().limit() << " " << 1.0/std::max(sqrt(smax),1.0e-3) << std::endl;
  if ( ipcol_fo != "NONE" ) {
    double sf(0.0), sw(0.0);
    for ( HRI ih = fphi.first(); !ih.last(); ih.next() )
      if ( !wrk_f1[ih].missing() && !fphi[ih].missing() ) {
        sf += wrk_f1[ih].f();
        sw += fphi[ih].f();
      }
    std::cout << " Fw/Fo " << sw/sf << std::endl;
  }
  if ( smax == 0.0 )
    { std::cerr << "No density provided." << std::endl; return 1; }

  // store copy of input model
  clipper::MiniMol mol_wrk_in = mol_wrk;

  // map stats
  natools.init_stats( xwrk );
  NautilusLog log( title ); // edited
  std::cout << std::endl;

  // initial phosphates from PDB
  if ( mol_pho.size() > 0 ) mol_wrk = natools.phosphate( xwrk, mol_wrk, mol_pho );

  for ( int cyc = 0; cyc < ncyc; cyc++ ) {
    std::cout << "Internal cycle " << clipper::String( cyc+1, 3 ) << std::endl << std::endl; // edited

    // adjust labels and label non-NA chains to keep
    mol_wrk = NucleicAcidTools::flag_chains( mol_wrk );

    // find chains
    mol_wrk = natools.find( xwrk, mol_wrk, nhit/2, nhit/2, srchst );
    log.log( "FIND", mol_wrk, verbose >= 5 );

    if (cyc == 0) {
        // Predictions predictions(ippredicted_phos_map);
        // MLFind ml_find = MLFind(predictions, xwrk);
        // ml_find.load_library_from_file(ippdb_ref);
        // clipper::MiniMol mol_find = ml_find.find();
        // log.log("FIND ML", mol_find, verbose >= 5);
        // mol_wrk = mol_find;
    }

    // grow chains
    mol_wrk = natools.grow( xwrk, mol_wrk, 25, 0.001 );
    log.log( "GROW", mol_wrk, verbose >= 5 );

    // join
    NucleicAcidJoin na_join;
    mol_wrk = na_join.join( mol_wrk );
    log.log( "JOIN", mol_wrk, verbose >= 5 );
    //for ( int c = 0; c < mol_wrk.size(); c++ ) { for ( int r = 0; r < mol_wrk[c].size(); r++ ) std::cout << mol_wrk[c][r].type().trim(); std::cout << std::endl; }

    // link
    mol_wrk = natools.link( xwrk, mol_wrk );
    log.log( "LINK", mol_wrk, verbose >= 5 );
    //for ( int c = 0; c < mol_wrk.size(); c++ ) { for ( int r = 0; r < mol_wrk[c].size(); r++ ) std::cout << mol_wrk[c][r].type().trim(); std::cout << std::endl; }

    // prune
    mol_wrk = natools.prune( mol_wrk );
    log.log( "PRUNE", mol_wrk, verbose >= 5 );
    //for ( int c = 0; c < mol_wrk.size(); c++ ) { for ( int r = 0; r < mol_wrk[c].size(); r++ ) std::cout << mol_wrk[c][r].type().trim(); std::cout << std::endl; }

    // rebuild
    mol_wrk = natools.rebuild_chain( xwrk, mol_wrk );
    log.log( "CHAIN", mol_wrk, verbose >= 5 );
    //for ( int c = 0; c < mol_wrk.size(); c++ ) { for ( int r = 0; r < mol_wrk[c].size(); r++ ) std::cout << mol_wrk[c][r].type().trim(); std::cout << std::endl; }

    // sequence
    NucleicAcidSequence na_seqnc;
    mol_wrk = na_seqnc.sequence( xwrk, mol_wrk, seq_wrk );
    log.log( "SEQNC", mol_wrk, verbose >= 5 );
    //for ( int c = 0; c < mol_wrk.size(); c++ ) { for ( int r = 0; r < mol_wrk[c].size(); r++ ) std::cout << mol_wrk[c][r].type().trim(); std::cout << std::endl; }

    // rebuild
    NucleicAcidRebuildBases na_bases;
    mol_wrk = na_bases.rebuild_bases( xwrk, mol_wrk );
    log.log( "BASES", mol_wrk, verbose >= 5 );
    //for ( int c = 0; c < mol_wrk.size(); c++ ) { for ( int r = 0; r < mol_wrk[c].size(); r++ ) std::cout << mol_wrk[c][r].type().trim(); std::cout << std::endl; }

    prog.summary_beg();
    msg = log.log_info( mol_wrk, true ); // edited
    std::cout << "Internal cycle " << clipper::String( cyc+1, 3 ) << std::endl << msg << std::endl ;
    prog.summary_end();

    // file output edited SWH Nov'17
    if ( opxml != "NONE" ) log.xml( opxml ); //, mol_wrk );
  }

  // move to match input model
  if ( mol_wrk_in.size() > 0 ){
    NucleicAcidTools::symm_match( mol_wrk, mol_wrk_in );
    log.log( "SYMMA", mol_wrk, verbose >=5 );
  }

  // set up residue types
  const clipper::String basetypes = "ACGTU";
  clipper::MiniMol mol_new( xwrk.spacegroup(), xwrk.cell() );
  for ( int c = 0; c < mol_wrk.size(); c++ ) {
    clipper::MPolymer mpx = mol_wrk[c];
    if ( !mpx.exists_property( "NON-NA" ) ) {
      bool dna = false;
      for ( int r = 0; r < mpx.size(); r++ ) {
        const clipper::String type = mpx[r].type().trim()+" ";
        const char ctype = type[0];
        int t = NucleicAcidTools::base_index( ctype );
        if ( t >= 0 ) mpx[r].set_type( "  "+basetypes.substr(t,1) );
        else          mpx[r].set_type( "  U" );
        if ( mpx[r].type().trim() == "T" ) dna = true;
      }
      // for DNA, prefix type with D and remove O2'
      if ( dna ) {
        for ( int r = 0; r < mpx.size(); r++ ) {
          const clipper::String type = " D" + mpx[r].type().trim();
          mpx[r].set_type( type );
          const int io2 = mpx[r].lookup( " O2'", clipper::MM::ANY );
          clipper::MMonomer mm;
          mm.copy( mpx[r], clipper::MM::COPY_MP );
          for ( int a = 0; a < mpx[r].size(); a++ ) {
            if ( a != io2 ) mm.insert( mpx[r][a] );
          }
          mpx[r] = mm;
        }
      }
    }
    mol_new.insert( mpx );
  }
  log.log( "SETRES", mol_wrk, verbose >=5 );

  // add true sequence numbers
  ModelTidy::chain_renumber( mol_new, seq_wrk );
  log.log( "TIDY", mol_new, verbose >= 5 );
  //new chain labelling routine, for 2-char label, SWH
  NucleicAcidTools::chain_label( mol_new, cifflag );
  log.log( "LABE", mol_new, verbose >= 5 );
  // final file output
  clipper::MMDBfile pdbfile;
  pdbfile.export_minimol( mol_new );
  pdbfile.write_file( oppdb, cifflag );
  msg = log.log_info( mol_new, true );	// added by SWH
  std::cout << "$TEXT:Result: $$ $$" << std::endl << msg << "\n$$" << std::endl; // added by SWH
  log.profile();
  prog.set_termination_message( "Normal termination" );
}
