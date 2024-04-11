/*! \file nautilus-util.cpp nautilus library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */

#include "nautilus-util.h"

#include <fstream>
#include <clipper/contrib/edcalc.h>
#include "nautilus-tools.h"


extern "C" {
#include <stdlib.h>
}


void NautilusUtil::set_reference( clipper::String& pdb )
{
  const char* clibdptr = getenv( "CLIBD" );
  const char* ccp4ptr = getenv( "CCP4" );
  if ( clibdptr != NULL ) {
    clipper::String clibd( clibdptr );
    clipper::String ccp4( ccp4ptr );
    clipper::String path;
    std::ifstream file;
    if ( pdb == "NONE" ) {
      path = clibd+"/nautilus_lib.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) pdb = path;
    }
    if ( pdb == "NONE" ) {
      path = clibd+"\\nautilus_lib.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) pdb = path;
    }
    if ( pdb == "NONE" ) {
      path = ccp4+"/share/nautilus_lib.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
           if ( !file.fail() ) pdb = path;
    }
    if ( pdb == "NONE" ) {
      path = ccp4+"\\share\\nautilus_lib.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) pdb = path;
    }
    if ( pdb == "NONE" ) 
      clipper::Message::message( clipper::Message_fatal( "No reference data specified and not in $CLIBD" ) );
  } else {
    clipper::Message::message( clipper::Message_fatal( "No reference data specified and $CLIBD not found" ) );
  }
}


#ifdef NAUTILUS_PROFILE
#include <sys/times.h>
#include <unistd.h>
void NautilusLog::log( const clipper::String& id )
{
  int i;
  tms tmst;
  times( &tmst );
  long ticks = sysconf(_SC_CLK_TCK);
  double cpu = double( tmst.tms_utime ) / double( ticks );
  double elapsed = cpu - currentcpu;
  if ( id != "" ) {
    for ( i = 0; i < prof.size(); i++ )
      if ( id == prof[i].first ) break;
    if ( i < prof.size() )
      prof[i].second += elapsed;
    else
      prof.push_back( std::pair<std::string,double>( id, elapsed ) );
  }
  currentcpu = cpu;
}
#else
void NautilusLog::log( const clipper::String& id ) {}
#endif


void NautilusLog::log( const clipper::String& id, const clipper::MiniMol& mol, bool view )
{
  clipper::String msg = log_info( mol, false );
  std::cout << id << ": " << msg << std::endl;
  if ( view ) {
    for ( int c = 0; c < mol.size(); c++ )
      std::cout << mol[c].size() << " ";
    std::cout << std::endl;
    for ( int c = 0; c < mol.size(); c++ ) {
      if ( !mol[c].exists_property( "NON-NA" ) ) {
        for ( int r1 = 0; r1 < mol[c].size()-1; r1++ ) {
          int r2 = r1 + 1;
          int a1 = mol[c][r1].lookup( " O3'", clipper::MM::ANY );
          int a2 = mol[c][r2].lookup( " O5'", clipper::MM::ANY );
          if ( a1 >= 0 && a2 >= 0 ) {
            double r = sqrt( ( mol[c][r1][a1].coord_orth() -
                               mol[c][r2][a2].coord_orth() ).lengthsq() );
            if ( r > 5.0 ) std::cout << "BREAK " << c << " " << r1 << " " << r2 << " " << r << std::endl;
          }
        }
      }
    }
  }
  log( id );
}


clipper::String NautilusLog::log_info( const clipper::MiniMol& mol, bool summary )
{
  //int nnc(0), nna(0);
  // added by me
  //clipper::MiniMol mol_wrk = mol;
  int nmax, nseq, nnc, nna;
  nnc = nna = nmax = nseq = 0;
  // count aa and chains
  // edited by SWH
  for ( int c = 0; c < mol.size(); c++ )
    if ( !mol[c].exists_property( "NON-NA" ) ) {
      if ( mol[c].size() > nmax ) nmax = mol[c].size();
      nnc += 1;
      nna += mol[c].size();     //count all built with or without neucleobases - gives sugar-phosphate backbone
      for ( int r = 0; r < mol[c].size(); r++ )
      {
          if ( mol[c][r].lookup( " C4 ", clipper::MM::ANY ) >= 0 ) // base ring
          {
              if ( mol[c][r].lookup( " O2 ", clipper::MM::ANY ) >= 0 ||
                      mol[c][r].lookup( " N6 ", clipper::MM::ANY ) >= 0 ||
                      mol[c][r].lookup( " O6 ", clipper::MM::ANY ) >= 0 ) nseq++;  // properly built nucleobases
          }
      }
    }
  // store every end of cycle, added SWH Nov'17
    cycdat dat;
    dat.nmax = nmax;
    dat.nchns = nnc;
    dat.nres = nna;
    dat.nseq = nseq;
    data.push_back( dat );

  if ( summary )
  {
    cycdat dat;
    dat.nmax = nmax;
    dat.nchns = nnc;
    dat.nres = nna;
    dat.nseq = nseq;
    data.push_back( dat );
    char s[1000];
    sprintf( s, " %6i sugar-phosphate residues built in %3i chains, the longest having %4i residues.\n "
            "%6i nucleic acids were sequenced.\n", nna, nnc, nmax, nseq);
    return clipper::String(s);
  }

  return clipper::String(nna,4) + " nucleic acids built in " +
         clipper::String(nnc,3) + " chains.";
}

int NautilusUtil::count_na( const clipper::MiniMol& mol )
{
  //int nnc(0), nna(0);
  // added by me
  //clipper::MiniMol mol_wrk = mol;
  int nmax, nseq, nnc, nna;
  nnc = nna = nmax = nseq = 0;
  // count aa and chains
  // edited by SWH
  for ( int c = 0; c < mol.size(); c++ )
    if ( !mol[c].exists_property( "NON-NA" ) ) {
      if ( mol[c].size() > nmax ) nmax = mol[c].size();
      nnc += 1;
      nna += mol[c].size();     //count all built with or without neucleobases - gives sugar-phosphate backbone
      for ( int r = 0; r < mol[c].size(); r++ )
      {
        if ( mol[c][r].lookup( " C4 ", clipper::MM::ANY ) >= 0 ) // base ring
        {
          if ( mol[c][r].lookup( " O2 ", clipper::MM::ANY ) >= 0 ||
                  mol[c][r].lookup( " N6 ", clipper::MM::ANY ) >= 0 ||
                  mol[c][r].lookup( " O6 ", clipper::MM::ANY ) >= 0 ) nseq++;  // properly built nucleobases
        }
      }
    }
  // store every end of cycle, added SWH Nov'17

  return nna;
}

float NautilusUtil::per_residue_rscc(clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, float res) {
    clipper::Cell cell = xmap.cell();
    clipper::Spacegroup spg = xmap.spacegroup();
    mol = NucleicAcidTools::flag_chains(mol);
    std::vector<clipper::Atom> atom_list = {};

    clipper::MiniMol na_only = {spg, cell};

    for (int p = 0; p < mol.size(); p++) {
        clipper::MPolymer mp;
        mp.set_id(mol[p].id());
        if (!mol[p].exists_property("NON-NA")) {
            for (int m = 0; m < mol[p].size(); m++) {
                mp.insert(mol[p][m]);
                for (int a = 0; a < mol[p][m].size(); a++) {
                    atom_list.emplace_back(mol[p][m][a]);
                }
            }
        }
        na_only.insert(mp);
    }

    clipper::MAtomNonBond neighbour_search = clipper::MAtomNonBond(na_only, 1.5);
    clipper::Coord_frac cf0( 0,0,0 );
    clipper::Coord_frac cf1( 1,1,1 );

    clipper::Resolution reso( res);
    clipper::Grid_sampling grid( spg, cell, reso );
    clipper::EDcalc_iso<float> maskcalc( res );
    clipper::Xmap<float> calc_map = {spg, cell, grid };
    maskcalc(calc_map, atom_list);

    clipper::Coord_grid g0( cf0.coord_grid(grid).u(), cf0.coord_grid(grid).v(), cf0.coord_grid(grid).w() );
    clipper::Coord_grid g1( cf1.coord_grid(grid).u(), cf1.coord_grid(grid).v(), cf1.coord_grid(grid).w() );

    clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;
    i0 = clipper::Xmap_base::Map_reference_coord( xmap, g0 );

    std::map<std::pair<int, int>, std::vector<std::pair<double, double>>> residue_pairs = {};

    for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
        for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
            for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
                clipper::Coord_orth co = iw.coord().coord_frac(xmap.grid_sampling()).coord_orth(cell);
                auto nearby = neighbour_search.atoms_near(co, 1.5);

                double min_distance = 1e5;
                clipper::MAtomIndexSymmetry min_atom;
                for (const auto& a: nearby) {
                    double distance = (co - na_only[a.polymer()][a.monomer()][a.atom()].coord_orth()).lengthsq();
                    if (distance < min_distance) {
                        min_distance = distance;
                        min_atom = a;
                    }
                }

                if (!nearby.empty()) {
                    std::pair<int, int> residue_key = std::make_pair(min_atom.polymer(), min_atom.monomer());
                    std::pair<double, double> point_pair = std::make_pair(xmap[iw], calc_map[iw]);

                    residue_pairs[residue_key].emplace_back(point_pair);
                }
            }


    for (const auto &residue_pair: residue_pairs) {
        int polymer = residue_pair.first.first;
        int monomer = residue_pair.first.second;
        const std::vector<std::pair<double, double>> &points = residue_pair.second;

        double sum_obs_residue = 0;
        double sum_calc_residue = 0;
        double count_obs_residue = 0;
        double count_calc_residue = 0;

        for (const auto &point_pair: points) {
            double obs_value = point_pair.first;
            double calc_value = point_pair.second;

            sum_obs_residue += obs_value;
            sum_calc_residue += calc_value;

            count_obs_residue++;
            count_calc_residue++;
        }

        double avg_obs = sum_obs_residue/count_obs_residue;
        double avg_calc = sum_calc_residue/count_calc_residue;

        double sum_delta = 0;
        double sum_delta_obs_sq = 0;
        double sum_delta_calc_sq = 0;

        for (const auto &point_pair: points) {
            double obs_delta = point_pair.first-avg_obs;
            double calc_delta = point_pair.second-avg_calc;

            sum_delta += obs_delta*calc_delta;
            sum_delta_obs_sq += pow(obs_delta, 2);
            sum_delta_calc_sq += pow(calc_delta, 2);
        }

        double rscc = sum_delta/sqrt(sum_delta_obs_sq*sum_delta_calc_sq);
        std::cout << polymer << " " << na_only[polymer][monomer].id() << " " << na_only[polymer][monomer].type() << " " << rscc << std::endl;
    }
}

std::map<std::pair<int, int>, double>
NautilusUtil::per_residue_rsrz(clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, float res) {
    clipper::Cell cell = xmap.cell();
    clipper::Spacegroup spg = xmap.spacegroup();
    mol = NucleicAcidTools::flag_chains(mol);
    std::vector<clipper::Atom> atom_list = {};

    clipper::MiniMol na_only = {spg, cell};

    for (int p = 0; p < mol.size(); p++) {
        clipper::MPolymer mp;
        mp.set_id(mol[p].id());
        if (!mol[p].exists_property("NON-NA")) {
            for (int m = 0; m < mol[p].size(); m++) {
                mp.insert(mol[p][m]);
                for (int a = 0; a < mol[p][m].size(); a++) {
                    atom_list.emplace_back(mol[p][m][a]);
                }
            }
        }
        na_only.insert(mp);
    }

    clipper::MAtomNonBond neighbour_search = clipper::MAtomNonBond(na_only, 1.5);
    clipper::Coord_frac cf0( 0,0,0 );
    clipper::Coord_frac cf1( 1,1,1 );

    clipper::Resolution reso( res);
    clipper::Grid_sampling grid( spg, cell, reso );
    clipper::EDcalc_iso<float> maskcalc( res );
    clipper::Xmap<float> calc_map = {spg, cell, grid };
    maskcalc(calc_map, atom_list);

    clipper::Coord_grid g0( cf0.coord_grid(grid).u(), cf0.coord_grid(grid).v(), cf0.coord_grid(grid).w() );
    clipper::Coord_grid g1( cf1.coord_grid(grid).u(), cf1.coord_grid(grid).v(), cf1.coord_grid(grid).w() );

    clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;
    i0 = clipper::Xmap_base::Map_reference_coord( xmap, g0 );

    std::map<std::pair<int, int>, std::vector<std::pair<double, double>>> residue_pairs = {};

    for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
        for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
            for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
                clipper::Coord_orth co = iw.coord().coord_frac(xmap.grid_sampling()).coord_orth(cell);
                auto nearby = neighbour_search.atoms_near(co, 1.5);

                double min_distance = 1e5;
                clipper::MAtomIndexSymmetry min_atom;
                for (const auto& a: nearby) {
                    double distance = (co - na_only[a.polymer()][a.monomer()][a.atom()].coord_orth()).lengthsq();
                    if (distance < min_distance) {
                        min_distance = distance;
                        min_atom = a;
                    }
                }

                if (!nearby.empty()) {
                    std::pair<int, int> residue_key = std::make_pair(min_atom.polymer(), min_atom.monomer());
                    std::pair<double, double> point_pair = std::make_pair(xmap[iw], calc_map[iw]);

                    residue_pairs[residue_key].emplace_back(point_pair);
                }
            }

    std::map<std::pair<int, int>, double> rsrs = {};

    for (const auto &residue_pair: residue_pairs) {
        const std::vector<std::pair<double, double>> &points = residue_pair.second;

        double numerator = 0;
        double denominator = 0;

        for (const auto &point_pair: points) {
            double obs = point_pair.first;
            double calc = point_pair.second;

            numerator += (obs-calc);
            denominator += (obs+calc);

        }

        double rsr = numerator/denominator;
        rsrs.insert({residue_pair.first, rsr});
    }

    double rsr_sum = 0.0;
    for(const auto& entry: rsrs) {
        rsr_sum += entry.second;
    }
    double rsr_mean = rsr_sum / rsrs.size();

    double squared_rsr_sum = 0.0;
    for(const auto& entry: rsrs) {
        squared_rsr_sum += std::pow(entry.second - rsr_mean, 2);
    }
    double rsr_std_dev = std::sqrt(squared_rsr_sum / rsrs.size());

    std::map<std::pair<int, int>, double> rsrzs = {};

    for (const auto& rsr_pair: rsrs) {
        double rsrz = (rsr_pair.second-rsr_mean)/rsr_std_dev;
        rsrzs[rsr_pair.first] = rsrz;
    }

    return rsrzs;
}

int NautilusUtil::count_well_modelled_nas(clipper::MiniMol &mol, clipper::Xmap<float>& xwrk, float res) {
    auto rsrzs = NautilusUtil::per_residue_rsrz(mol, xwrk, res);
    int count = 0;
    for (const auto& rsr_pair: rsrzs) {
        if (rsr_pair.second >= -1) {count += 1;}
    }
    return count;
}

float NautilusUtil::calculate_rscc(clipper::MiniMol&mol, const clipper::Xmap<float>& xmap, float res) {
  clipper::Cell cell = xmap.cell();
  clipper::Spacegroup spg = xmap.spacegroup();
  // clipper::Atom_list atom_list = mol.atom_list();
  mol = NucleicAcidTools::flag_chains(mol);

  std::vector<clipper::Atom> atom_list = {};
  for (int p = 0; p < mol.size(); p++) {
    if (!mol[p].exists_property("NON-NA")) {
      for (int m = 0; m < mol[p].size(); m++) {
        for (int a = 0; a < mol[p][m].size(); a++) {
          atom_list.emplace_back(mol[p][m][a]);
        }
      }
    }
  }

  if (atom_list.empty()) {
    std::cout << "No nucleic acid found, RSCC = 0" << std::endl;
    return 0;
  }
  clipper::Range<clipper::ftype> urange, vrange, wrange;

  for (const auto & atom : atom_list) {
    clipper::Coord_frac cf = atom.coord_orth().coord_frac( cell );
    urange.include( cf.u() );
    vrange.include( cf.v() );
    wrange.include( cf.w() );
  }
  clipper::Coord_frac cf0( urange.min(), vrange.min(), wrange.min() );
  clipper::Coord_frac cf1( urange.max(), vrange.max(), wrange.max() );

  clipper::Resolution reso( res);
  clipper::Grid_sampling grid( spg, cell, reso );
  clipper::EDcalc_iso<float> maskcalc( res );
  clipper::Xmap<float> calc_map = {spg, cell, grid };
  maskcalc(calc_map, atom_list);

  clipper::Coord_grid g0( cf0.coord_grid(grid).u(), cf0.coord_grid(grid).v(), cf0.coord_grid(grid).w() );
  clipper::Coord_grid g1( cf1.coord_grid(grid).u(), cf1.coord_grid(grid).v(), cf1.coord_grid(grid).w() );

  clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;
  i0 = clipper::Xmap_base::Map_reference_coord( xmap, g0 );

  double sum_obs = 0;
  double sum_calc = 0;
  double count_obs = 0;
  double count_calc = 0;

  for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
    for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
      for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
        sum_obs += xmap[iw];
        count_obs += 1;

        sum_calc += calc_map[iw];
        count_calc += 1;
      }

  double avg_obs = sum_obs/count_obs;
  double avg_calc = sum_calc/count_calc;

  double sum_delta = 0;
  double sum_delta_obs_sqrd = 0;
  double sum_delta_calc_sqrd = 0;

  for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
    for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
      for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
        sum_delta += (xmap[iw]-avg_obs)*(calc_map[iw]-avg_calc);

        sum_delta_obs_sqrd += pow(xmap[iw]-avg_obs, 2);
        sum_delta_calc_sqrd += pow(calc_map[iw]-avg_calc, 2);
      }

  double rscc = (sum_delta)/(sqrt(sum_delta_obs_sqrd*sum_delta_calc_sqrd));
  std::cout << "RSCC = " << rscc << std::endl;
  return rscc;
}


void NautilusLog::xml( const clipper::String& file ) const //, const clipper::MiniMol& mol )
{
  /*
  int nres, nseq, nchn, nmax;
  nchn = mol.size();
  nres = nseq = nmax = 0;
  for ( int c = 0; c < mol.size(); c++ ) {
    if ( !mol[c].exists_property( "NON-NA" ) ) {
      if ( mol[c].size() > nmax ) nmax = mol[c].size();
      for ( int r = 0; r < mol[c].size(); r++ ) {
        if ( mol[c][r].lookup( " C4 ", clipper::MM::ANY ) >= 0 ) nres++;    //no C4 means no nitrogenous base
        if ( mol[c][r].lookup( " C4 ", clipper::MM::ANY ) >= 0 && mol[c][r].lookup( " O4 ", clipper::MM::ANY ) >= 0 ) nseq++;
      }
    }
  }*/
  // xml output each cycle and summary, added SWH Nov'17
  std::ofstream f;
  f.open( file.c_str(), std::ios::out );
  f << "<NautilusResult>" << std::endl;
  f << " <Title>" << title_.c_str() << "</Title>" << std::endl;
  f << " <Cycles>" << std::endl;
  for ( int c = 0; c < data.size() ; c++ )
  {
      f << "  <Cycle>" << std::endl;
      f << "   <Number>" << c+1 << "</Number>" << std::endl;
      f << "   <FragmentsBuilt>" << data[c].nchns << "</FragmentsBuilt>" << std::endl;
      f << "   <ResiduesBuilt>" << data[c].nres << "</ResiduesBuilt>" << std::endl;
      f << "   <ResiduesSequenced>" << data[c].nseq << "</ResiduesSequenced>" << std::endl;
      f << "   <ResiduesLongestFragment>" << data[c].nmax << "</ResiduesLongestFragment>" << std::endl;
      f << "  </Cycle>" << std::endl;
  }
  f << " </Cycles>" << std::endl;
  f << " <Final>" << std::endl;
  int c = data.size()-1;
  f << "  <Number>" << c+1 << "</Number>" << std::endl;
  f << "  <FragmentsBuilt>" << data[c].nchns << "</FragmentsBuilt>" << std::endl;
  f << "  <ResiduesBuilt>" << data[c].nres << "</ResiduesBuilt>" << std::endl;
  f << "  <ResiduesSequenced>" << data[c].nseq << "</ResiduesSequenced>" << std::endl;
  f << "  <ResiduesLongestFragment>" << data[c].nmax << "</ResiduesLongestFragment>" << std::endl;
  f << " </Final>" << std::endl;
  f << "</NautilusResult>" << std::endl;
  f.close();
}


void NautilusLog::profile()
{
  if ( prof.size() > 0 ) {
    std::cout << std::endl << "Profile:" << std::endl;
    for ( int i = 0; i < prof.size(); i++ )
      std::cout << prof[i].first << ": " << clipper::String( prof[i].second, 8 ) << " s" << std::endl;
  }
}

