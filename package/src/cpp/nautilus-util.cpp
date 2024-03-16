/*! \file nautilus-util.cpp nautilus library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */

#include "nautilus-util.h"

#include <fstream>


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

