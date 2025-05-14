/*! \file nautilus-target.h nautilus library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#ifndef NAUTILUS_TARGET_H
#define NAUTILUS_TARGET_H


#include "nucleicacid_db.h"
#include "nautilus-ss-find.h"
#include "nautilus-findml.h"

class NucleicAcidTarget {
 public:
  typedef std::pair<clipper::Coord_orth,clipper::Coord_orth> Coords;
  typedef std::vector<Coords> Target;
  typedef std::vector<clipper::Coord_orth> Standard;
  NucleicAcidTarget() {}
  void init( const float c_hi[][3], const float c_lo[][3], const float c_repr[3][3], const int ncoord );
  void init_stats( const clipper::Xmap<float>& xmap );
  const Target& target() const { return target_; }
  const Standard& standard() const { return standard_; }
  double radius() const;

  float score_min( const clipper::Xmap<float>& xmap, const clipper::RTop_orth& rtop ) const;
  float score_sum( const clipper::Xmap<float>& xmap, const clipper::RTop_orth& rtop ) const;
  float cutoff_min( double p ) const;
  float cutoff_sum( double p ) const;

 private:
  Target target_;
  std::vector<clipper::Coord_orth> standard_;
  std::vector<float> smin, ssum;
};

class NucleicAcidTargets {
 public:
  NucleicAcidTargets();
  void add_pdb( const clipper::String& file );
  void init_stats( const clipper::Xmap<float>& xmap );
  const NucleicAcidDB::Chain& db() const { return nadb; } 
  static void superpose_sugar( NucleicAcidDB::Chain& frag, int posn, const NucleicAcidDB::NucleicAcid& na );
  float score_sugar( const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na ) const;
  static float score_sugar_from_predictions(const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid &na) ;

  float score_phosphate( const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na1, const NucleicAcidDB::NucleicAcid& na2 ) const;
  NucleicAcidDB::NucleicAcid next_na_group(const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na) const;
  NucleicAcidDB::NucleicAcid prev_na_group(const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na ) const;
  const NucleicAcidTarget& target_sugar() const { return target_s; }
  const NucleicAcidTarget& target_phosphate() const { return target_p; }

  const NucleicAcidDB::Chain join_sugars( const clipper::Xmap<float>& xmap, const NucleicAcidDB::NucleicAcid& na1, const NucleicAcidDB::NucleicAcid& na2, int len, double rmsdlim ) const;

  const clipper::MiniMol phosphate( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol, const clipper::MiniMol& mol_pho );
  const clipper::MiniMol find( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol, int nsugar, int nphosp, double step );
  const clipper::MiniMol grow( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol, int ngrow, double fcut) const;
  const clipper::MiniMol link( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol ) const;
  const clipper::MiniMol prune( clipper::MiniMol& mol ) const;
  const clipper::MiniMol strong_prune( clipper::MiniMol& mol ) const;
  const clipper::MiniMol rebuild_chain( const clipper::Xmap<float>& xmap, const clipper::MiniMol& mol ) const;

 private:
  NucleicAcidDB::Chain nadb;
  NucleicAcidDB::NucleicAcid narepr;
  NucleicAcidTarget target_s, target_p;
  std::vector<SearchResult> found_s, found_p;
};

#endif
