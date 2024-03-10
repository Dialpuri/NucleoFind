/*! \file nautilus-tidy.cpp nautilus library */
/* (C) 2013 Kevin Cowtan & University of York all rights reserved */

#include "nautilus-tidy.h"


//#include <algorithm>
#include <string>


std::vector<int> ModelTidy::chain_renumber( clipper::MiniMol& mol, const clipper::MMoleculeSequence& seq )
{
  std::vector<int> result( mol.size(), -1 );
  for ( int c = 0; c < mol.size(); c++ )
    result[c] = chain_renumber( mol[c], seq );
  return result;
}


int ModelTidy::chain_renumber( clipper::MPolymer& mp, const clipper::MMoleculeSequence& seq )
{
  // initial numbering
  for ( int i = 0; i < mp.size(); i++ ) mp[i].set_seqnum( i+1 );
  if ( seq.size() == 0 ) return -1;

  // convert sequences to unique strings
  clipper::String chnseq = ModelTidy::chain_sequence( mp );
  std::vector<clipper::String> seqs( seq.size() );
  for ( int chn = 0; chn < seq.size(); chn++ ) {
    clipper::String s = "";
    for ( int res = 0; res < seq[chn].sequence().length(); res++ )
      s += tolower( seq[chn].sequence()[res] );
    seqs[chn] = s;
  }
  
  // now find best match
  int bestchn = -1;
  int bestscr = 0;
  clipper::MSequenceAlign align( clipper::MSequenceAlign::LOCAL,
                                   1.0, -1000.0, -4.0 );
  std::pair<std::vector<int>,std::vector<int> > result;
  for ( int chn = 0; chn < seqs.size(); chn++ ) {
    const clipper::String& seqseq = seqs[chn];
    result = align( chnseq, seqseq );
    int scr = 0;
    for ( int i = 0; i < result.first.size(); i++ ) {
      if ( result.first[i] >= 0 ) {
        if ( chnseq[i] == seqseq[result.first[i]] )
           scr++;
          else
            scr--;
      }
    }
    if ( scr > bestscr ) {
      bestchn = chn;
      bestscr = scr;
    }
  }
  if ( bestchn < 0 ) return -1;
  
  // now number residues
  clipper::String truseq = seqs[bestchn];
  result = align( chnseq, truseq );
  std::vector<int> nums = result.first;
  
  /*
  std::cout << bestchn << " " << bestscr << std::endl;
  for ( int i = 0; i < chnseq.size(); i++ )
    std::cout << chnseq[i];
  std::cout << std::endl;
  for ( int i = 0; i < chnseq.size(); i++ )
    std::cout << (nums[i]>=0&&nums[i]<truseq.length() ? truseq[nums[i]] : '-');
  std::cout << std::endl;
  */
  
  // find bounds of sequenced region
  int i0, i1, i;
  for ( i0 = 0;    i0 < nums.size(); i0++ ) if ( nums[i0] >= 0 ) break;
  for ( i1 = nums.size()-1; i1 >= 0; i1-- ) if ( nums[i1] >= 0 ) break;
  if ( i0 < nums.size() )
    for ( i = i0 - 1; i >= 0;          i-- ) nums[i] = nums[i+1] - 1;
  if ( i1 >= 0 )
    for ( i = i1 + 1; i < nums.size(); i++ ) nums[i] = nums[i-1] + 1;
  
  // renumber the model, with inscodes if necessary
  const clipper::String inscodes = " ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  const int l = inscodes.length();
  for ( i = 0; i < nums.size(); i++ ) nums[i] = nums[i] * l;
  for ( i = 1; i < nums.size(); i++ )
    if ( nums[i] <= nums[i-1] ) nums[i] = nums[i-1]+1;
  for ( i = 0; i < mp.size(); i++ ) {
    const int nres = (nums[i]/l) + 1;
    const int nins = clipper::Util::mod( nums[i], l );
    if ( nins == 0 ) mp[i].set_seqnum( nres );
    else             mp[i].set_seqnum( nres, inscodes.substr( nins, 1 ) );
  }

  /*
  for ( int i = 0; i < chnseq.size()-1; i++ )
    if ( nums[i+1] != nums[i]+1 ) std::cout << "! " << mp.id() << " " << nums[i] << " " << nums[i+1] << std::endl;
  */
  return bestchn;
}


clipper::String ModelTidy::chain_sequence( const clipper::MPolymer& mp )
{
  const int NTYPE = 27;
  const char rtype1[NTYPE] =
    {  'A',  'R',  'N',  'D',  'C',  'Q',  'E',  'G',  'H',  'I',
       'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V',
       'M',  'a',  'c',  'g',  't',  'u',  '?'};
  const char rtype3[NTYPE][4] =
    {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
     "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
     "MSE","A","C","G","T","U","UNK"};
  clipper::String seq = "";
  for ( int res = 0; res < mp.size(); res++ ) {
    clipper::String typ = mp[res].type().trim();
    char c = ' ';
    for ( int t = 0; t < NTYPE; t++ )
      if ( typ == clipper::String(rtype3[t]) )
        c = rtype1[t];
    if ( c == ' ' && mp[res].type().length() > 0 )
      c = char( mp[res].type()[0] + 128 );  // use dummy symbols for unknowns
    seq += c;
  }
  return seq;
}
