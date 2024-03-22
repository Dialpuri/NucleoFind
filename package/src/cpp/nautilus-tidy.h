/*! \file buccaneer-tidy.h buccaneer library */
/* (C) 2002-2010 Kevin Cowtan & University of York all rights reserved */


#ifndef NAUTILUS_TIDY
#define NAUTILUS_TIDY

#include "nautilus-util.h"


//! Class for tidying model
class ModelTidy {
 public:
  static std::vector<int> chain_renumber( clipper::MiniMol& mol, const clipper::MMoleculeSequence& seq );
 private:
  static int chain_renumber( clipper::MPolymer& mp, const clipper::MMoleculeSequence& seq );
  static clipper::String chain_sequence( const clipper::MPolymer& mp );
};


#endif
