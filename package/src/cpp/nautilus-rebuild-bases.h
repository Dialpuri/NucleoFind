/*! \file nautilus-rebuild-bases.h nautilus library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#ifndef NAUTILUS_REBUILD_BASES_H
#define NAUTILUS_REBUILD_BASES_H


#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>


class NucleicAcidRebuildBases {
 public:
  NucleicAcidRebuildBases();
  clipper::MiniMol rebuild_bases( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol ) const;

 private:
  std::vector<std::vector<std::pair<clipper::String,clipper::Coord_orth> > > basedata;
};


#endif
