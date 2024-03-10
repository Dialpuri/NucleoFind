/*! \file nautilus-sequence.h nautilus library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#ifndef NAUTILUS_SEQUENCE_H
#define NAUTILUS_SEQUENCE_H


#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>


class NucleicAcidSequence {
 public:
  static double sequence_similarity( const clipper::String& seq1, const clipper::String& seq2 );
  static clipper::RTop_orth base_rtop( const clipper::MMonomer& mm, const clipper::Xmap<float>& xmap );
  static std::vector<std::vector<double> > score( const clipper::Xmap<float>& xmap, std::vector<clipper::RTop_orth>& rtops );
  static std::pair<double,clipper::String> sequence_score( const std::vector<std::vector<double> >& scores, const clipper::String& subseq );
  static clipper::String sequence_chain( const clipper::Xmap<float>& xmap, const clipper::MChain& chain, const clipper::MMoleculeSequence& seq );
  clipper::MiniMol sequence( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol, const clipper::MMoleculeSequence& seq ) const;
};


#endif
