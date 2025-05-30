/*! Nautilus tools main chain database */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#include "nautilus-tools.h"

#include <algorithm>


int NucleicAcidTools::bindex[256], NucleicAcidTools::bindext[256];


NucleicAcidTools::NucleicAcidTools()
{
  for ( int i = 0; i < 256; i++ ) {
    char c = char(i);
    int t = -1;
    if ( c == 'A' ) t = 0;
    if ( c == 'C' ) t = 1;
    if ( c == 'G' ) t = 2;
    if ( c == 'T' ) t = 3;
    if ( c == 'U' ) t = 4;
    bindex[i]  = t;
    bindext[i] = std::min(t,3);
  }
}


// static functions

clipper::MiniMol NucleicAcidTools::flag_chains( const clipper::MiniMol& mol )
{
  // flag any chain containing at least one non-NA
  clipper::MiniMol mol_new = mol;
  for ( int c = 0; c < mol_new.size(); c++ ) {
    bool flg = false;
    for ( int r = 0; r < mol_new[c].size(); r++ ) {
      NucleicAcidDB::NucleicAcid na( mol_new[c][r] );
      if ( na.flag() == NucleicAcidDB::NucleicAcid::NONE ) {
        flg = true;
        break;
      }
    }
    if ( flg )
      mol_new[c].set_property( "NON-NA", clipper::Property<bool>( true ) );
  }
  // for the remaining chains, adjust names to a single character
  for ( int c = 0; c < mol_new.size(); c++ ) {
    if ( !mol_new[c].exists_property( "NON-NA" ) ) {
      std::vector<int> flag( mol_new[c].size(), 1 );  // flag NAs as good
      for ( int r = 0; r < mol_new[c].size(); r++ ) {
        clipper::String type = mol_new[c][r].type().trim();
        while ( type.length() > 1 ) type = type.substr(1);
        mol_new[c][r].set_type( type );
        if ( type == "U" ) {  // detect unknown NAs
          if ( mol_new[c][r].lookup( " O4 ", clipper::MM::ANY ) < 0 ) {
            if ( mol_new[c][r].lookup( " C4 ", clipper::MM::ANY ) < 0 )
              flag[r] =  0;  // not base atoms: maybe unknown
            else
              flag[r] = -1;  // C4 but no O4: definitely unknown
          }
        }
      }
      if ( flag[0] == 0 )             flag[0] = -1;
      if ( flag[flag.size()-1] == 0 ) flag[flag.size()-1] = -1;
      for ( int r = 1; r < flag.size()-1; r++ )
        if ( flag[r] == 0 && flag[r-1] == -1 ) flag[r] = -1;
      for ( int r = flag.size()-2; r > 0; r-- )
        if ( flag[r] == 0 && flag[r+1] == -1 ) flag[r] = -1;
      for ( int r = 0; r < mol_new[c].size(); r++ )
        if ( flag[r] == -1 ) mol_new[c][r].set_type( "?" );
    }
  }
  return mol_new;
}


clipper::RTop_orth NucleicAcidTools::symmetry_rtop( const std::vector<clipper::Coord_orth>& cowrk, clipper::Coord_orth& coref, const clipper::Spacegroup& spgr, const clipper::Cell& cell )
{
  std::vector<clipper::Coord_frac> cwrk( cowrk.size() );
  for ( int a = 0; a < cowrk.size(); a++ )
    cwrk[a] = cowrk[a].coord_frac(cell);
  clipper::Coord_frac cref = coref.coord_frac(cell);
  clipper::Coord_frac c1, c2;
  double d2, d2min(1.0e12);
  int smin(0);
  clipper::Coord_frac dmin(0.0,0.0,0.0);
  for ( int s = 0; s < spgr.num_symops(); s++ )
    for ( int a = 0; a < cwrk.size(); a++ ) {
      c1 = ( spgr.symop(s) * cwrk[a] );
      c2 = c1.lattice_copy_near( cref );
      d2 = ( c2 - cref ).lengthsq( cell );
      if ( d2 < d2min ) {
        d2min = d2;
        smin = s;
        dmin = c2 - c1;
      }
    }
  clipper::RTop_frac rf( spgr.symop(smin).rot(), spgr.symop(smin).trn()+dmin );
  return rf.rtop_orth( cell );
}


clipper::MiniMol NucleicAcidTools::chain_sort( const clipper::MiniMol& mol )
{
  std::vector<std::pair<int,int> > chnsiz( mol.size() );
  for ( int chn = 0; chn < mol.size(); chn++ )
    chnsiz[chn] = std::pair<int,int>( -mol[chn].size(), chn );
  std::sort( chnsiz.begin(), chnsiz.end() );
  clipper::MiniMol molnew( mol.spacegroup(), mol.cell() );
  for ( int chn = 0; chn < mol.size(); chn++ )
    molnew.insert( mol[chnsiz[chn].second] );
  return molnew;
}


clipper::Coord_orth NucleicAcidTools::coord_adjust( const clipper::Coord_orth& co, const clipper::Coord_orth& cc3, const clipper::Coord_orth& cf3, const clipper::Coord_orth& cc4, const clipper::Coord_orth& cf4, double rad )
{
  if ( co.is_null() ) return co;
  clipper::Coord_orth result = co;
  double w3 = 1.0 - sqrt( ( co - cf3 ).lengthsq() ) / rad;
  double w4 = 1.0 - sqrt( ( co - cf4 ).lengthsq() ) / rad;
  if ( w3 > 0.0 ) result += w3 * ( cc3 - cf3 );
  if ( w4 > 0.0 ) result += w4 * ( cc4 - cf4 );
  return result;
}


bool NucleicAcidTools::symm_match( clipper::MiniMol& molwrk, const clipper::MiniMol& molref )
{
  clipper::Spacegroup spg1 = clipper::Spacegroup(clipper::Spacegroup::P1);
  clipper::Spacegroup spgr = molwrk.spacegroup();
  clipper::Cell       cell = molwrk.cell();

  // calculate extent of model
  clipper::Atom_list atomr = molref.atom_list();
  clipper::Range<clipper::ftype> urange, vrange, wrange;
  clipper::Coord_frac cfr( 0.0, 0.0, 0.0 );
  for ( int i = 0; i < atomr.size(); i++ ) {
    clipper::Coord_frac cf = atomr[i].coord_orth().coord_frac( cell );
    cfr += cf;
    urange.include( cf.u() );
    vrange.include( cf.v() );
    wrange.include( cf.w() );
  }
  clipper::Coord_frac cf0( urange.min(), vrange.min(), wrange.min() );
  clipper::Coord_frac cf1( urange.max(), vrange.max(), wrange.max() );
  cfr = (1.0/double(atomr.size())) * cfr;

  // calculate mask using wrk cell and ref atoms
  clipper::Resolution reso( 5.0 );
  clipper::Grid_sampling grid( spg1, cell, reso );
  clipper::Grid_range    grng( grid,  cf0,  cf1 );
  grng.add_border(4);
  clipper::NXmap<float> nxmap( cell, grid, grng ), nxflt( cell, grid, grng );
  clipper::EDcalc_mask<float> maskcalc( 2.0 );
  nxmap = 0.0;
  maskcalc( nxmap, atomr );
  MapFilterFn_g5 fn;
  clipper::MapFilter_fft<float>
    fltr( fn, 1.0, clipper::MapFilter_fft<float>::Relative );
  fltr( nxflt, nxmap );

  // now score each chain, symmetry and offset in turn
  for ( int c = 0; c < molwrk.size(); c++ ) {
    double              bestscr = 0.0;
    int                 bestsym = 0;
    clipper::Coord_frac bestoff( 0.0, 0.0, 0.0 );
    const clipper::Coord_frac cfh( 0.5, 0.5, 0.5 );
    for ( int sym = 0; sym < spgr.num_symops(); sym++ ) {
      clipper::Atom_list atomw = molwrk[c].atom_list();
      clipper::RTop_orth rtop = spgr.symop(sym).rtop_orth( cell );
      clipper::Coord_orth cow( 0.0, 0.0, 0.0 );
      for ( int a = 0; a < atomw.size(); a++ ) {
        atomw[a].transform( rtop );
        cow += atomw[a].coord_orth();
      }
      if ( atomw.size() > 0 ) cow = (1.0/double(atomw.size())) * cow;
      clipper::Coord_frac cfw = cow.coord_frac( cell );
      clipper::Coord_frac cfwt = cfw.lattice_copy_near( cfr - cfh );
      clipper::Coord_frac off0 = cfwt - cfw;

      // try offsets
      for ( double du = 0.0; du <= 1.01; du += 1.0 )
        for ( double dv = 0.0; dv < 1.01; dv += 1.0 )
          for ( double dw = 0.0; dw < 1.01; dw += 1.0 ) {
            clipper::Coord_frac off( rint( off0.u() ) + du,
                                     rint( off0.v() ) + dv,
                                     rint( off0.w() ) + dw );
            clipper::Coord_orth ofo = off.coord_orth( cell );
            double scr = 0.0;
            for ( int a = 0; a < atomw.size(); a++ ) {
              clipper::Coord_orth coa = atomw[a].coord_orth() + ofo;
              clipper::Coord_grid cga = nxflt.coord_map( coa ).coord_grid();
              if ( nxflt.in_map( cga ) ) scr += nxflt.get_data( cga );
            }
            if ( scr > bestscr ) {
              bestscr = scr;
              bestsym = sym;
              bestoff = off;
            }
          }
    }
    // now transform using the best operator
    clipper::Coord_orth cot = bestoff.coord_orth( cell );
    clipper::RTop_orth rtop = spgr.symop(bestsym).rtop_orth( cell );
    rtop = clipper::RTop_orth( rtop.rot(), rtop.trn()+cot );
    molwrk[c].transform( rtop );
  }

  return true;
}


bool NucleicAcidTools::chain_label( clipper::MiniMol& mol, clipper::MMDBManager::TYPE cifflag )
{
  // set up default chain labels
  std::vector<clipper::String> labels;
  labels.push_back( "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz" );
  labels.push_back( "0123456789" );
  bool lbl_alpha[53][52] = {{0}};
  bool lbl_num[10][10]   = {{0}};

  // get existing labels
  for ( int chn = 0; chn < mol.size(); chn++ )
  {
    std::pair<int, int> index;
    index = get_usedlabels(mol[chn].id(), labels);
    if ( index.second - 52 < 0 )
      lbl_alpha[index.first][index.second] = 1;
    else
      lbl_num[index.first-1][index.second] = 1;
  }
    
  // label chains
  int label = 0;
  for ( int chn = 0; chn < mol.size(); chn++ ) {
    if ( mol[chn].id() == "" ) {
      bool newlabelled = false;
      do{
        newlabelled = false;
        int r, c;
        // label chains with letters
        if ( label < labels[0].length() ) {
          if (!lbl_alpha[0][label])
          {
            mol[chn].set_id( labels[0].substr( label, 1 ) );
            newlabelled = true;
            label++;
          }
          else label++;
        }
        else
        {
          if ( cifflag == clipper::MMDBManager::CIF )
          {
            if (label < 2756 )
            {
              c = ( label - labels[0].length() ) % labels[0].length();
              r = ( label - labels[0].length() ) / labels[0].length();
              if (!lbl_alpha[r][c])
              {
                mol[chn].set_id( labels[0].substr( r, 1 ) + labels[0].substr( c, 1 ) ); 
                newlabelled = true;
                label++;
              }
              else label++;
            }
            else
            {
              if ( label < 2856 )
              {
                c = ( label - 2756 ) % labels[1].length();
                r = ( label - 2756 ) / labels[1].length();
              }
              else
              {
                c = ( label - 2756 ) % labels[1].length();
                r = ( label - 2756 ) / labels[1].length() % labels[1].length();
              }
              if (!lbl_num[r][c])
              {
                clipper::String newid="";
                if (r==0)
                  newid = labels[1].substr( c, 1 ); 
                else
                  newid = labels[1].substr( r, 1 ) + labels[1].substr( c, 1 );
                int rmax = 1;
                for ( int f = 0; f < mol.size(); f++)
                  if ( mol[f].id() == newid )
                    rmax = std::max(rmax, mol[f][mol[f].size()-1].seqnum()+5);
                mol[chn].set_id( newid );
                for ( int res = 0; res < mol[chn].size(); res++)
                  mol[chn][res].set_seqnum( rmax + res );
                newlabelled = true;
                label++;
              }
              else label++;
            }
          }
          else
          {
            c = ( label - labels[0].length() ) % labels[1].length();
            if (!lbl_num[0][c])
            {
              // pack remaining residues into numbered chains
              clipper::String newid = labels[1].substr( c, 1 );
              int rmax = 1;
              for ( int f = 0; f < mol.size(); f++)
                if ( mol[f].id() == newid )
                  rmax = std::max(rmax, mol[f][mol[f].size()-1].seqnum()+5);
              mol[chn].set_id( newid );
              for ( int res = 0; res < mol[chn].size(); res++)
                mol[chn][res].set_seqnum( rmax + res );
              newlabelled = true;
              label++;
            }
            else label++;
          }
        }
      }while(!newlabelled);
    }
  }
  return true;

}


std::pair<int, int> NucleicAcidTools::get_usedlabels(clipper::String chainid, std::vector<clipper::String> labels)
{
    int ind[2] = {-1, -1};
    int row = -1, column = -1;
    // check chars in chain id
    for (int i=0; i<chainid.length(); i++)
      for ( int v=0; v<1; v++)
        for ( int j =0; j<labels[v].length();j++)
          if ( chainid[i] == labels[v][j]) ind[i] = j;

    if( ind[1] == -1) // single char chain id
    {
      row = 0;
      column = ind[0];
    }
    else // double char chain id
    {
      row = ind[0] + 1; // offset the first row of single char
      column = ind[1];
    }
    return std::make_pair(row, column);
}


void NucleicAcidTools::residue_label(clipper::MiniMol &mol) {
    for (int p = 0; p < mol.size(); p++) {
        for (int m = 0; m < mol[p].size(); m++) {
            mol[p][m].set_id(m);
        }
    }
}

/*
bool NucleicAcidTools::chain_label( clipper::MiniMol& mol )
{
  // set up default chain labels
  std::vector<clipper::String> labels;
  labels.push_back( "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz" );
  labels.push_back( "0123456789" );

  // get existing labels
  clipper::String chainids = "";
  for(int chn=0; chn<mol.size(); chn++)
     chainids = chainids + mol[chn].id();

  // eliminate used labels
  for (int i = 0; i < 1; i++)
  {
    clipper::String newlabels;
    for ( int j = 0; j < labels[i].length(); j++)
    {
     if ( chainids.find( labels[i].substr(j,1) ) == clipper::String::npos)
        newlabels += labels[i].substr(j,1);
    }
    labels[i] = newlabels;
  }
  // label chains
  int label = 0;
  std::vector<int> nresc( labels[1].length(), 0 );
  for ( int chn = 0; chn < mol.size(); chn++ ) {
    if ( label < labels[0].length() ) { // label chains with letters
      if ( mol[chn].id() == "" ) {
        mol[chn].set_id( labels[0].substr( label, 1 ) );
        label++;
      }
    } else { // pack remaining into numbered chains
      int c = label - labels[0].length();
      int c1 = c % labels[1].length();
      mol[chn].set_id( labels[1].substr( c1, 1 ) );
      for ( int res = 0; res < mol[chn].size(); res++ )
        mol[chn][res].set_seqnum( res + nresc[c1] + 1 );
      nresc[c1] += mol[chn].size()+5;
      label++;
    }
  }
  return true;
}
*/
