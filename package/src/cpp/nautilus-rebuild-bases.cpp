/*! Nautilus rebuild bases */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#include "nautilus-rebuild-bases.h"

#include "nautilus-tools.h"
#include "nautilus-sequence.h"
#include "nucleicacid_db.h"


NucleicAcidRebuildBases::NucleicAcidRebuildBases()
{
  const char baseatoms[6][11][3] = {
    {"N1","C2","N3","C4","C5","C6","N6","N7","C8","N9",""  },
    {"N1","C2","O2","N3","C4","N4","C5","C6",""  ,""  ,""  },
    {"N1","C2","N2","N3","C4","C5","C6","O6","N7","C8","N9"},
    {"N1","C2","O2","N3","C4","O4","C5","C6","C7",""  ,""  },
    {"N1","C2","O2","N3","C4","O4","C5","C6",""  ,""  ,""  },
    {"N1","C2","N3","C4","C5","C6",""  ,""  ,""  ,""  ,""  },
  };

  const double basecoords[6][11][3] = {
    {{ 2.987, 0.000,-2.697},{ 1.706, 0.000,-3.099},{ 0.590, 0.000,-2.367},{ 0.883, 0.000,-1.055},{ 2.145, 0.000,-0.488},{ 3.244, 0.000,-1.369},{ 4.516, 0.000,-0.964},{ 2.061, 0.000, 0.897},{ 0.772, 0.000, 1.134},{ 0.000, 0.000, 0.000},{ 0.000, 0.000, 0.000}},
    {{-0.000, 0.000,-0.000},{ 0.755, 0.000,-1.182},{ 0.165, 0.000,-2.272},{ 2.106, 0.000,-1.103},{ 2.705, 0.000, 0.089},{ 4.039, 0.000, 0.116},{ 1.963, 0.000, 1.306},{ 0.625, 0.000, 1.216},{ 0.000, 0.000, 0.000},{ 0.000, 0.000, 0.000},{ 0.000, 0.000, 0.000}},
    {{ 2.931, 0.000,-2.645},{ 1.648, 0.000,-3.118},{ 1.521, 0.000,-4.453},{ 0.569, 0.000,-2.346},{ 0.886, 0.000,-1.044},{ 2.138, 0.000,-0.463},{ 3.285, 0.000,-1.305},{ 4.472, 0.000,-0.991},{ 2.037, 0.000, 0.917},{ 0.753, 0.000, 1.147},{ 0.000, 0.000, 0.000}},
    {{ 0.000, 0.000, 0.000},{ 0.716, 0.000,-1.178},{ 0.186, 0.000,-2.270},{ 2.079, 0.000,-1.029},{ 2.784, 0.000, 0.154},{ 4.016, 0.000, 0.123},{ 1.971, 0.000, 1.334},{ 0.639, 0.000, 1.218},{ 2.800, 0.000, 2.300},{ 0.000, 0.000, 0.000},{ 0.000, 0.000, 0.000}},
    {{ 0.000, 0.000, 0.000},{ 0.716, 0.000,-1.178},{ 0.186, 0.000,-2.270},{ 2.079, 0.000,-1.029},{ 2.784, 0.000, 0.154},{ 4.016, 0.000, 0.123},{ 1.971, 0.000, 1.334},{ 0.639, 0.000, 1.218},{ 0.000, 0.000, 0.000},{ 0.000, 0.000, 0.000},{ 0.000, 0.000, 0.000}},
    {{ 0.000, 0.000, 0.000},{ 0.716, 0.000,-1.178},{ 2.079, 0.000,-1.029},{ 2.784, 0.000, 0.154},{ 1.971, 0.000, 1.334},{ 0.639, 0.000, 1.218},{ 0.000, 0.000, 0.000},{ 0.000, 0.000, 0.000},{ 0.000, 0.000, 0.000},{ 0.000, 0.000, 0.000},{ 0.000, 0.000, 0.000}},
  };

  for ( int i = 0; i < 6; i++ ) {
    std::vector<std::pair<clipper::String,clipper::Coord_orth> > atoms;
    for ( int j = 0; j < 11; j++ )
      if ( clipper::String(baseatoms[i][j]) != "" )
        atoms.push_back( std::pair<clipper::String,clipper::Coord_orth>( baseatoms[i][j], clipper::Coord_orth( basecoords[i][j][0], basecoords[i][j][1], basecoords[i][j][2] ) ) );
    basedata.push_back( atoms );
  }
}


clipper::MiniMol NucleicAcidRebuildBases::rebuild_bases( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol ) const
{
  clipper::MiniMol mol_new = mol;

  // build bases in default position
  for ( int c = 0; c < mol_new.size(); c++ ) {
    if ( !mol_new[c].exists_property( "NON-NA" ) ) {
      for ( int r = 0; r < mol_new[c].size(); r++ ) {
        NucleicAcidDB::NucleicAcid na( mol_new[c][r] );
        if ( na.flag() == NucleicAcidDB::NucleicAcid::COMPLETE ) {
          const clipper::String type = mol_new[c][r].type()+"?";
          const char ctype = type.trim()[0];
          const int t = (ctype!='?') ? NucleicAcidTools::base_index(ctype) : 5;
          if ( t >= 0 ) {
            const clipper::RTop_orth rtop =
              NucleicAcidSequence::base_rtop( mol_new[c][r], xmap );
            if ( !rtop.is_null() ) {
              clipper::MMonomer mm0, mm1;
              mm0 = na.mmonomer();
              mm1.set_type( mm0.type() );
              for ( int a = 0; a < mm0.size(); a++ )
                if ( mm0[a].id().trim() != "N1" ) mm1.insert( mm0[a] );
              clipper::MAtom ma = clipper::MAtom::null();
              ma.set_occupancy( 1.0 ); ma.set_u_iso( 0.25 );
              for ( int a = 0; a < basedata[t].size(); a++ ) {
                ma.set_id( " "+basedata[t][a].first+" " );
                ma.set_element( basedata[t][a].first.substr(0,1) );
                ma.set_coord_orth( rtop * basedata[t][a].second );
                mm1.insert( ma );
              }
              mol_new[c][r] = mm1;
            }
          }
        }
      }
    }
  }
  // prune any clashing bases
  const char clashatoms[5][5] = {" C2 "," N3 "," C4 "," C5 "," C6 "};
  const int nclashatoms = sizeof(clashatoms)/sizeof(clashatoms[0]);
  clipper::MAtomNonBond nb( mol_new, 4.0 );
  for ( int c1 = 0; c1 < mol_new.size(); c1++ ) {
    if ( !mol_new[c1].exists_property( "NON-NA" ) ) {
      for ( int r1 = 0; r1 < mol_new[c1].size(); r1++ ) {
        bool clash = false;
        for ( int a1 = 0; a1 < mol_new[c1][r1].size(); a1++ ) {
          bool test = false;
          for ( int i = 0; i < nclashatoms; i++ )
            if ( mol_new[c1][r1][a1].id() == clashatoms[i] )
              test = true;
          if ( test ) {
            const clipper::Coord_orth co = mol_new[c1][r1][a1].coord_orth();
            std::vector<clipper::MAtomIndexSymmetry> atoms = nb( co, 2.0 );
            for ( int i = 0; i < atoms.size(); i++ ) {
              int c2 = atoms[i].polymer();
              int r2 = atoms[i].monomer();
              if ( c1 != c2 || r1 != r2 ) clash = true;
            }
          }
        }
        if ( clash )
          mol_new[c1][r1] =
            NucleicAcidDB::NucleicAcid( mol_new[c1][r1] ).mmonomer();
      }
    }
  }

  // now add the remaining atoms
  for ( int c = 0; c < mol_new.size(); c++ ) {
    if ( !mol_new[c].exists_property( "NON-NA" ) ) {
      // insert OP1, OP2
      clipper::MAtom ma = clipper::MAtom::null();
      ma.set_occupancy( 1.0 ); ma.set_u_iso( 0.25 ); ma.set_element( "O" );
      for ( int r = 1; r < mol_new[c].size(); r++ ) {
        const int i1 = mol_new[c][r].lookup(" OP1",clipper::MM::ANY);
        const int i2 = mol_new[c][r].lookup(" OP2",clipper::MM::ANY);
        const int i3 = mol_new[c][r-1].lookup(" O3'",clipper::MM::ANY);
        const int ip = mol_new[c][r].lookup(" P  ",clipper::MM::ANY);
        const int i5 = mol_new[c][r].lookup(" O5'",clipper::MM::ANY);
        if ( i1 < 0 && i2 < 0 && i3 >= 0 && ip >= 0 && i5 >= 0 ) {
          const clipper::Coord_orth c3 = mol_new[c][r-1][i3].coord_orth();
          const clipper::Coord_orth cp = mol_new[c][r][ip].coord_orth();
          const clipper::Coord_orth c5 = mol_new[c][r][i5].coord_orth();
          const clipper::Coord_orth u = c3 - cp;
          const clipper::Coord_orth v = c5 - cp;
          if ( u.lengthsq() < 3.0 && v.lengthsq() < 3.0 ) {
            const clipper::Coord_orth x( clipper::Vec3<>::cross(u,v).unit() );
            const clipper::Coord_orth y( (u.unit()+v.unit()).unit() );
            const clipper::Coord_orth co1 = cp - 0.80*y + 1.25*x;
            const clipper::Coord_orth co2 = cp - 0.80*y - 1.25*x;
            ma.set_coord_orth( co2 ); ma.set_id( " OP2" );
            mol_new[c][r].insert( ma, i5 );
            ma.set_coord_orth( co1 ); ma.set_id( " OP1" );
            mol_new[c][r].insert( ma, i5 );
          }
        }
      }
      // insert O2'
      for ( int r = 0; r < mol_new[c].size(); r++ ) {
        const int io = mol_new[c][r].lookup(" O2'",clipper::MM::ANY);
        const int i1 = mol_new[c][r].lookup(" C1'",clipper::MM::ANY);
        const int i2 = mol_new[c][r].lookup(" C2'",clipper::MM::ANY);
        const int i3 = mol_new[c][r].lookup(" C3'",clipper::MM::ANY);
        if ( io < 0 && i1 >= 0 && i2 >= 0 && i3 >= 0 ) {
          const clipper::Coord_orth c1 = mol_new[c][r][i1].coord_orth();
          const clipper::Coord_orth c2 = mol_new[c][r][i2].coord_orth();
          const clipper::Coord_orth c3 = mol_new[c][r][i3].coord_orth();
          const clipper::Coord_orth u = c1 - c2;
          const clipper::Coord_orth v = c3 - c2;
          if ( u.lengthsq() < 3.0 && v.lengthsq() < 3.0 ) {
            const clipper::Coord_orth x( clipper::Vec3<>::cross(u,v).unit() );
            const clipper::Coord_orth y( (u.unit()+v.unit()).unit() );
            const clipper::Coord_orth co = c2 - 0.80*y - 1.20*x;
            ma.set_coord_orth( co ); ma.set_id( " O2'" );
            mol_new[c][r].insert( ma, i1 );
          }
        }
      }
    }
  }

  return mol_new;
}
