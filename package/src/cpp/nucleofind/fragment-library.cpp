//
// Created by Jordan Dialpuri on 23/04/2025.
//

#include "fragment-library.h"

void NucleoFind::TriNucleotide::setup(clipper::MMonomer &m1, clipper::MMonomer &m2, clipper::MMonomer &m3) {
    int ip1 = m1.lookup( " P  ", clipper::MM::ANY);
    int ip2 = m2.lookup( " P  ", clipper::MM::ANY);
    int ip3 = m3.lookup( " P  ", clipper::MM::ANY);

    if (ip1 == -1 && ip2 == -1 && ip3 == -1) {
        throw std::runtime_error("CriticalError: Library file is missing phosphate atoms");
    }

    P1 = m1[ip1].coord_orth();
    P2 = m2[ip2].coord_orth();
    P3 = m3[ip3].coord_orth();
}

void NucleoFind::TriNucleotideLibrary::add_library(const std::string &library_path) {
    const int mmdbflags = ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors | ::mmdb::MMDBF_IgnoreRemarks;
    clipper::MMDBfile mfile;
    clipper::MiniMol mol;
    mfile.SetFlag( mmdbflags );
    mfile.read_file( library_path );
    mfile.import_minimol( mol );

    for ( int c = 0; c < mol.size(); c++ ) {
        clipper::MMonomer monomer1 = mol[c][0];
        clipper::MMonomer monomer2 = mol[c][1];
        clipper::MMonomer monomer3 = mol[c][2];
        library.emplace_back( monomer1, monomer2, monomer3 );
    }

}
