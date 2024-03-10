/*! NucleicAcidDB Top 500 main chain database */
/* (C) 2008-2009 Kevin Cowtan & University of York all rights reserved */


#include "nucleicacid_db.h"
#include <fstream>


namespace NucleicAcidDB {

    NucleicAcidFull::NucleicAcidFull(const clipper::MMonomer &monomer_1) {
        P =    clipper::Coord_orth(clipper::Vec3<>::null());
        OP1 =  clipper::Coord_orth(clipper::Vec3<>::null());
        OP2 =  clipper::Coord_orth(clipper::Vec3<>::null());
        O5p1 = clipper::Coord_orth(clipper::Vec3<>::null());
        C5p1 = clipper::Coord_orth(clipper::Vec3<>::null());
        C4p1 = clipper::Coord_orth(clipper::Vec3<>::null());
        O4p1 = clipper::Coord_orth(clipper::Vec3<>::null());
        C3p1 = clipper::Coord_orth(clipper::Vec3<>::null());
        O3p1 = clipper::Coord_orth(clipper::Vec3<>::null());
        C2p1 = clipper::Coord_orth(clipper::Vec3<>::null());
        C1p1 = clipper::Coord_orth(clipper::Vec3<>::null());
        C2_1 = clipper::Coord_orth(clipper::Vec3<>::null());
        C4_1 = clipper::Coord_orth(clipper::Vec3<>::null());
        C5_1 = clipper::Coord_orth(clipper::Vec3<>::null());
        C6_1 = clipper::Coord_orth(clipper::Vec3<>::null());
        C8_1 = clipper::Coord_orth(clipper::Vec3<>::null());
        N1_1 = clipper::Coord_orth(clipper::Vec3<>::null());
        N2_1 = clipper::Coord_orth(clipper::Vec3<>::null());
        N3_1 = clipper::Coord_orth(clipper::Vec3<>::null());
        N4_1 = clipper::Coord_orth(clipper::Vec3<>::null());
        N6_1 = clipper::Coord_orth(clipper::Vec3<>::null());
        N7_1 = clipper::Coord_orth(clipper::Vec3<>::null());
        N9_1 = clipper::Coord_orth(clipper::Vec3<>::null());
        O2_1 = clipper::Coord_orth(clipper::Vec3<>::null());
        O6_1 = clipper::Coord_orth(clipper::Vec3<>::null());

        int ip = monomer_1.lookup(" P  ", clipper::MM::ANY);
        int io5p1 = monomer_1.lookup(" O5'", clipper::MM::ANY);
        int ic5p1 = monomer_1.lookup(" C5'", clipper::MM::ANY);
        int ic4p1 = monomer_1.lookup(" C4'", clipper::MM::ANY);
        int io4p1 = monomer_1.lookup(" O4'", clipper::MM::ANY);
        int ic3p1 = monomer_1.lookup(" C3'", clipper::MM::ANY);
        int io3p1 = monomer_1.lookup(" O3'", clipper::MM::ANY);
        int ic2p1 = monomer_1.lookup(" C2'", clipper::MM::ANY);
        int ic1p1 = monomer_1.lookup(" C1'", clipper::MM::ANY);
        int in9_1 = monomer_1.lookup(" N9 ", clipper::MM::ANY);
        int in7_1 = monomer_1.lookup(" N7 ", clipper::MM::ANY);
        int in6_1 = monomer_1.lookup(" N6 ", clipper::MM::ANY);
        int in4_1 = monomer_1.lookup(" N4 ", clipper::MM::ANY);
        int in3_1 = monomer_1.lookup(" N3 ", clipper::MM::ANY);
        int in2_1 = monomer_1.lookup(" N2 ", clipper::MM::ANY);
        int in1_1 = monomer_1.lookup(" N1 ", clipper::MM::ANY);
        int ic8_1 = monomer_1.lookup(" C8 ", clipper::MM::ANY);
        int ic6_1 = monomer_1.lookup(" C6 ", clipper::MM::ANY);
        int ic5_1 = monomer_1.lookup(" C5 ", clipper::MM::ANY);
        int ic4_1 = monomer_1.lookup(" C4 ", clipper::MM::ANY);
        int ic2_1 = monomer_1.lookup(" C2 ", clipper::MM::ANY);
        int io2_1 = monomer_1.lookup(" O2 ", clipper::MM::ANY);
        int io6_1 = monomer_1.lookup(" O6 ", clipper::MM::ANY);

//        std::cout << ic8_1 << " " << in9_1 << std::endl;

        if (ip >= 0) { P = monomer_1[ip].coord_orth(); }
        if (io5p1 >= 0) { O5p1 = monomer_1[io5p1].coord_orth(); }
        if (ic5p1 >= 0) { C5p1 = monomer_1[ic5p1].coord_orth(); }
        if (ic4p1 >= 0) { C4p1 = monomer_1[ic4p1].coord_orth(); }
        if (io4p1 >= 0) { O4p1 = monomer_1[io4p1].coord_orth(); }
        if (io3p1 >= 0) { O3p1 = monomer_1[io3p1].coord_orth(); }
        if (ic3p1 >= 0) { C3p1 = monomer_1[ic3p1].coord_orth(); }
        if (ic2p1 >= 0) { C2p1 = monomer_1[ic2p1].coord_orth(); }
        if (ic1p1 >= 0) { C1p1 = monomer_1[ic1p1].coord_orth(); }
        if (in1_1 >= 0) { N1_1 = monomer_1[in1_1].coord_orth(); }
        if (in2_1 >= 0) { N2_1 = monomer_1[in2_1].coord_orth(); }
        if (in3_1 >= 0) { N3_1 = monomer_1[in3_1].coord_orth(); }
        if (in4_1 >= 0) { N4_1 = monomer_1[in4_1].coord_orth(); }
        if (in6_1 >= 0) { N6_1 = monomer_1[in6_1].coord_orth(); }
        if (in7_1 >= 0) { N7_1 = monomer_1[in7_1].coord_orth(); }
        if (in9_1 >= 0) { N9_1 = monomer_1[in9_1].coord_orth(); }
        if (ic2_1 >= 0) { C2_1 = monomer_1[ic2_1].coord_orth(); }
        if (ic4_1 >= 0) { C4_1 = monomer_1[ic4_1].coord_orth(); }
        if (ic5_1 >= 0) { C5_1 = monomer_1[ic5_1].coord_orth(); }
        if (ic6_1 >= 0) { C6_1 = monomer_1[ic6_1].coord_orth(); }
        if (ic8_1 >= 0) { C8_1 = monomer_1[ic8_1].coord_orth(); }
        if (io2_1 >= 0) { O2_1 = monomer_1[io2_1].coord_orth(); }
        if (io6_1 >= 0) { O6_1 = monomer_1[io6_1].coord_orth(); }

        base_type = monomer_1.type();

    }

    clipper::MMonomer NucleicAcidFull::get_mmonomer(int id) {

        clipper::MMonomer monomer;
        monomer.set_type("A");
        monomer.set_id(id);

        if (!P.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(P, "P", "P");
            monomer.insert(atom);
        }

        if (!C1p1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(C1p1, "C1'", "C");
            monomer.insert(atom);
        }
        if (!C2p1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(C2p1, "C2'", "C");
            monomer.insert(atom);
        }
        if (!O3p1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(O3p1, "O3'", "O");
            monomer.insert(atom);
        }
        if (!C3p1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(C3p1, "C3'", "C");
            monomer.insert(atom);
        }
        if (!O4p1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(O4p1, "O4'", "O");
            monomer.insert(atom);
        }
        if (!C4p1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(C4p1, "C4'", "C");
            monomer.insert(atom);
        }
        if (!C5p1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(C5p1, "C5'", "C");
            monomer.insert(atom);
        }
        if (!O5p1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(O5p1, "O5'", "O");
            monomer.insert(atom);
        }

        if (!C2_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(C2_1, "C2", "C");
            monomer.insert(atom);
        }
        if (!C4_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(C4_1, "C4", "C");
            monomer.insert(atom);
        }
        if (!C5_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(C5_1, "C5", "C");
            monomer.insert(atom);
        }
        if (!C6_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(C6_1, "C6", "C");
            monomer.insert(atom);
        }
        if (!C8_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(C8_1, "C8", "C");
            monomer.insert(atom);
        }
        if (!N1_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(N1_1, "N1", "N");
            monomer.insert(atom);
        }
        if (!N2_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(N2_1, "N2", "N");
            monomer.insert(atom);
        }
        if (!N3_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(N3_1, "N3", "N");
            monomer.insert(atom);
        }
        if (!N4_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(N4_1, "N4", "N");
            monomer.insert(atom);
        }
        if (!N6_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(N6_1, "N6", "N");
            monomer.insert(atom);
        }
        if (!N7_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(N7_1, "N7", "N");
            monomer.insert(atom);
        }
        if (!N9_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(N9_1, "N9", "N");
            monomer.insert(atom);
        }
        if (!O2_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(O2_1, "O2", "O");
            monomer.insert(atom);
        }
        if (!O6_1.is_null()) {
            clipper::MAtom atom = NautilusUtil::create_atom(O6_1, "O6", "O");
            monomer.insert(atom);
        }



        return monomer;

    }

    void NucleicAcidFull::debug() {
        std::cout
                << "P" << P.format() << "\n"
                << "OP1." << OP1.format() << "\n"
                << "OP2." << OP2.format() << "\n"
                << "O5p1" << O5p1.format() << "\n"
                << "C5p1" << C5p1.format() << "\n"
                << "C4p1" << C4p1.format() << "\n"
                << "O4p1" << O4p1.format() << "\n"
                << "C3p1" << C3p1.format() << "\n"
                << "O3p1" << O3p1.format() << "\n"
                << "C2p1" << C2p1.format() << "\n"
                << "C1p1" << C1p1.format() << "\n"
                << "C2_1" << C2_1.format() << "\n"
                << "C4_1" << C4_1.format() << "\n"
                << "C5_1" << C5_1.format() << "\n"
                << "C6_1" << C6_1.format() << "\n"
                << "C8_1" << C8_1.format() << "\n"
                << "N1_1" << N1_1.format() << "\n"
                << "N2_1" << N2_1.format() << "\n"
                << "N3_1" << N3_1.format() << "\n"
                << "N4_1" << N4_1.format() << "\n"
                << "N6_1" << N6_1.format() << "\n"
                << "N7_1" << N7_1.format() << "\n"
                << "N9_1" << N9_1.format() << "\n"
                << "O2_1" << O2_1.format() << "\n"
                << "O6_1" << O6_1.format() << "\n" ;
    }


    bool ChainFull::add_pdb(const clipper::String &file, bool strict) {
        const int mmdbflags =
                ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors |
                ::mmdb::MMDBF_IgnoreRemarks | ::mmdb::MMDBF_AllowDuplChainID | ::mmdb::MMDBF_IgnoreSegID;
        clipper::MMDBfile mfile;
        clipper::MiniMol mol;
        mfile.SetFlag(mmdbflags);
        mfile.read_file(file);
        mfile.import_minimol(mol);
        if (mol.size() == 0) {
            return false;
        }

//        std::cout << mol.size() << std::endl;

        for (int c = 0; c < mol.size(); c++) {
            clipper::MPolymer mp;
            // select monomers by occupancy
//            std::cout << "mol[c].size() = " << mol[c].size() << std::endl;
            for (int r = 0; r < mol[c].size(); r++) {
                if (
                        mol[c][r].lookup(" C1'", clipper::MM::ANY) >= 0 &&
                        mol[c][r].lookup(" C2'", clipper::MM::ANY) >= 0 &&
                        mol[c][r].lookup(" C3'", clipper::MM::ANY) >= 0 &&
                        mol[c][r].lookup(" C4'", clipper::MM::ANY) >= 0 &&
                        mol[c][r].lookup(" C5'", clipper::MM::ANY) >= 0 &&
                        mol[c][r].lookup(" O3'", clipper::MM::ANY) >= 0 &&
                        mol[c][r].lookup(" O4'", clipper::MM::ANY) >= 0 &&
                        mol[c][r].lookup(" O5'", clipper::MM::ANY) >= 0 &&
                        mol[c][r].lookup(" P  ", clipper::MM::ANY) >= 0) {

                    int a = mol[c][r].lookup(" P  ", clipper::MM::ANY);
                    if (mol[c][r][a].occupancy() > 0.01 && mol[c][r][a].u_iso() < clipper::Util::b2u(100.0)) {
                        mp.insert(mol[c][r]);
                    }
                }
            }



//            for (int r = 0; r < mp.size(); r++) {
//                int a = mp[r].lookup(" P  ", clipper::MM::UNIQUE);
//                clipper::RTop_orth rt(clipper::Mat33<>::identity(), -mp[r][a].coord_orth());
//                mp[r].transform(rt);
//            }
//
//            clipper::MPolymer out_mp;
//            out_mp.insert(mp[0]);
//
//            clipper::MModel model;
//            model.insert(out_mp);
//            clipper::MiniMol mol;
//            mol.model() = model;
//
//              clipper::MMDBfile mfile;
//              mfile.export_minimol( mol );
//              mfile.write_file( "debug/pdb/database.pdb" );

            // now add the chain to the db
//            std::cout << mp.size() << std::endl;
            for (int r = 0; r < mp.size(); r++) {
                NucleicAcidFull rp(mp[r]);
//                std::cout << r << "/" << mp.size() << mp[r].type() << "\t" << rp.base_type <<  std::endl;
//                if (rp.flag() == NucleicAcid::COMPLETE && rp.type() != ' ') {
//          std::cout << "Adding monomer to db " << std::endl;
                add_monomer(rp);
//                } else {
//                    std::cout << "MONOMER INCOMPLETE" << std::endl;
//                }
            }
//            std::cout << "299" << std::endl;
        }

        return true;
    }

    ChainFull ChainFull::extract(int offset, int len) const {
        ChainFull dbc;
        for (int i = 0; i < len; i++) dbc.add_monomer(dbmonomers[offset + i]);
        return dbc;
    }

    bool ChainFull::is_continuous() const {
        // go through and find elements where there is a chain break
        const double dmin = 2.0;
        std::vector<bool> cterm(dbmonomers.size(), false);
        for (int i = 0; i < dbmonomers.size() - 1; i++) {
            int j = i + 1;
            const clipper::Coord_orth co1 = dbmonomers[i].O3p1;
            const clipper::Coord_orth co2 = dbmonomers[j].P;
            if (co1.is_null() || co2.is_null()) return false;
            const double d2 = (co1 - co2).lengthsq();
            if (d2 > dmin * dmin) return false;
        }
        return true;
    }

    void NucleicAcidFull::transform(const clipper::RTop_orth &rtop) {
        if (!P.is_null()) P = rtop * P;
        if (!OP1.is_null()) OP1 = rtop * OP1;
        if (!OP2.is_null()) OP2 = rtop * OP2;
        if (!O5p1.is_null()) O5p1 = rtop * O5p1;
        if (!C5p1.is_null()) C5p1 = rtop * C5p1;
        if (!C4p1.is_null()) C4p1 = rtop * C4p1;
        if (!O4p1.is_null()) O4p1 = rtop * O4p1;
        if (!C3p1.is_null()) C3p1 = rtop * C3p1;
        if (!O3p1.is_null()) O3p1 = rtop * O3p1;
        if (!C2p1.is_null()) C2p1 = rtop * C2p1;
        if (!C1p1.is_null()) C1p1 = rtop * C1p1;
        if (!C2_1.is_null()) C2_1 = rtop * C2_1;
        if (!C4_1.is_null()) C4_1 = rtop * C4_1;
        if (!C5_1.is_null()) C5_1 = rtop * C5_1;
        if (!C6_1.is_null()) C6_1 = rtop * C6_1;
        if (!C8_1.is_null()) C8_1 = rtop * C8_1;
        if (!N1_1.is_null()) N1_1 = rtop * N1_1;
        if (!N2_1.is_null()) N2_1 = rtop * N2_1;
        if (!N3_1.is_null()) N3_1 = rtop * N3_1;
        if (!N4_1.is_null()) N4_1 = rtop * N4_1;
        if (!N6_1.is_null()) N6_1 = rtop * N6_1;
        if (!N7_1.is_null()) N7_1 = rtop * N7_1;
        if (!N9_1.is_null()) N9_1 = rtop * N9_1;
        if (!O2_1.is_null()) O2_1 = rtop * O2_1;
        if (!O6_1.is_null()) O6_1 = rtop * O6_1;
    }



    clipper::Atom_list NucleicAcidFull::return_atom_list() {

        std::vector<clipper::Atom> atoms;
        clipper::MMonomer monomer = this->get_mmonomer();
        atoms.reserve(monomer.size());
        for (int a = 0; a < monomer.size(); a++) {
            atoms.emplace_back(monomer[a]);
        }

        return {atoms};
    }


    void ChainFull::transform(const clipper::RTop_orth &rtop) {
        for (int r = 0; r < dbmonomers.size(); r++)
            dbmonomers[r].transform(rtop);
    }

NucleicAcid::NucleicAcid( const clipper::Coord_orth& cp, const clipper::Coord_orth& co5, const clipper::Coord_orth& cc5, const clipper::Coord_orth& cc4, const clipper::Coord_orth& co4, const clipper::Coord_orth& cc3, const clipper::Coord_orth& co3, const clipper::Coord_orth& cc2, const clipper::Coord_orth& cc1, const clipper::Coord_orth& cn, const clipper::String& type )
{
  clipper::String t = type + "?";
  typ = t.trim()[0];
  clipper::Util::set_null( p_x );
  clipper::Util::set_null( o5x );
  clipper::Util::set_null( c5x );
  clipper::Util::set_null( c4x );
  clipper::Util::set_null( o4x );
  clipper::Util::set_null( c3x );
  clipper::Util::set_null( o3x );
  clipper::Util::set_null( c2x );
  clipper::Util::set_null( c1x );
  clipper::Util::set_null( n_x );
  if ( !cp.is_null()  ) { p_x = cp.x(); p_y = cp.y(); p_z = cp.z(); }
  if ( !co5.is_null() ) { o5x = co5.x(); o5y = co5.y(); o5z = co5.z(); }
  if ( !cc5.is_null() ) { c5x = cc5.x(); c5y = cc5.y(); c5z = cc5.z(); }
  if ( !cc4.is_null() ) { c4x = cc4.x(); c4y = cc4.y(); c4z = cc4.z(); }
  if ( !co4.is_null() ) { o4x = co4.x(); o4y = co4.y(); o4z = co4.z(); }
  if ( !cc3.is_null() ) { c3x = cc3.x(); c3y = cc3.y(); c3z = cc3.z(); }
  if ( !co3.is_null() ) { o3x = co3.x(); o3y = co3.y(); o3z = co3.z(); }
  if ( !cc2.is_null() ) { c2x = cc2.x(); c2y = cc2.y(); c2z = cc2.z(); }
  if ( !cc1.is_null() ) { c1x = cc1.x(); c1y = cc1.y(); c1z = cc1.z(); }
  if ( !cn.is_null()  ) { n_x = cn.x(); n_y = cn.y(); n_z = cn.z(); }
  set_flag();
}

NucleicAcid::NucleicAcid( const clipper::MMonomer& mm )
{
  clipper::String t = mm.type() + "?";
  typ = t.trim()[0];
  clipper::Util::set_null( p_x );
  clipper::Util::set_null( o5x );
  clipper::Util::set_null( c5x );
  clipper::Util::set_null( c4x );
  clipper::Util::set_null( o4x );
  clipper::Util::set_null( c3x );
  clipper::Util::set_null( o3x );
  clipper::Util::set_null( c2x );
  clipper::Util::set_null( c1x );
  clipper::Util::set_null( n_x );
  int ip  = mm.lookup( " P  ", clipper::MM::ANY );
  int io5 = mm.lookup( " O5'", clipper::MM::ANY );
  int ic5 = mm.lookup( " C5'", clipper::MM::ANY );
  int ic4 = mm.lookup( " C4'", clipper::MM::ANY );
  int io4 = mm.lookup( " O4'", clipper::MM::ANY );
  int ic3 = mm.lookup( " C3'", clipper::MM::ANY );
  int io3 = mm.lookup( " O3'", clipper::MM::ANY );
  int ic2 = mm.lookup( " C2'", clipper::MM::ANY );
  int ic1 = mm.lookup( " C1'", clipper::MM::ANY );
  int in  = mm.lookup( " N9 ", clipper::MM::ANY );
  if ( in < 0 ) in = mm.lookup( " N1 ", clipper::MM::ANY );
  if ( ip  >= 0 ) {
    p_x = mm[ip].coord_orth().x();
    p_y = mm[ip].coord_orth().y();
    p_z = mm[ip].coord_orth().z();
  }
  if ( io5 >= 0 ) {
    o5x = mm[io5].coord_orth().x();
    o5y = mm[io5].coord_orth().y();
    o5z = mm[io5].coord_orth().z();
  }
  if ( ic5 >= 0 ) {
    c5x = mm[ic5].coord_orth().x();
    c5y = mm[ic5].coord_orth().y();
    c5z = mm[ic5].coord_orth().z();
  }
  if ( ic4 >= 0 ) {
    c4x = mm[ic4].coord_orth().x();
    c4y = mm[ic4].coord_orth().y();
    c4z = mm[ic4].coord_orth().z();
  }
  if ( io4 >= 0 ) {
    o4x = mm[io4].coord_orth().x();
    o4y = mm[io4].coord_orth().y();
    o4z = mm[io4].coord_orth().z();
  }
  if ( ic3 >= 0 ) {
    c3x = mm[ic3].coord_orth().x();
    c3y = mm[ic3].coord_orth().y();
    c3z = mm[ic3].coord_orth().z();
  }
  if ( io3 >= 0 ) {
    o3x = mm[io3].coord_orth().x();
    o3y = mm[io3].coord_orth().y();
    o3z = mm[io3].coord_orth().z();
  }
  if ( ic2 >= 0 ) {
    c2x = mm[ic2].coord_orth().x();
    c2y = mm[ic2].coord_orth().y();
    c2z = mm[ic2].coord_orth().z();
  }
  if ( ic1 >= 0 ) {
    c1x = mm[ic1].coord_orth().x();
    c1y = mm[ic1].coord_orth().y();
    c1z = mm[ic1].coord_orth().z();
  }
  if ( in  >= 0 ) {
    n_x = mm[in].coord_orth().x();
    n_y = mm[in].coord_orth().y();
    n_z = mm[in].coord_orth().z();
  }
  set_flag();
}

clipper::Coord_orth NucleicAcid::coord_p () const
{
  if ( clipper::Util::is_null( p_x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( p_x, p_y, p_z );
}

clipper::Coord_orth NucleicAcid::coord_o5() const
{
  if ( clipper::Util::is_null( o5x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( o5x, o5y, o5z );
}

clipper::Coord_orth NucleicAcid::coord_c5() const
{
  if ( clipper::Util::is_null( c5x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( c5x, c5y, c5z );
}

clipper::Coord_orth NucleicAcid::coord_c4() const
{
  if ( clipper::Util::is_null( c4x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( c4x, c4y, c4z );
}

clipper::Coord_orth NucleicAcid::coord_c3() const
{
  if ( clipper::Util::is_null( c3x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( c3x, c3y, c3z );
}

clipper::Coord_orth NucleicAcid::coord_o3() const
{
  if ( clipper::Util::is_null( o3x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( o3x, o3y, o3z );
}

clipper::Coord_orth NucleicAcid::coord_c2() const
{
  if ( clipper::Util::is_null( c2x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( c2x, c2y, c2z );
}

clipper::Coord_orth NucleicAcid::coord_c1() const
{
  if ( clipper::Util::is_null( c1x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( c1x, c1y, c1z );
}

clipper::Coord_orth NucleicAcid::coord_o4() const
{
  if ( clipper::Util::is_null( o4x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( o4x, o4y, o4z );
}

clipper::Coord_orth NucleicAcid::coord_n () const
{
  if ( clipper::Util::is_null( n_x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( n_x, n_y, n_z );
}

clipper::MMonomer NucleicAcid::mmonomer() const
{
  clipper::MMonomer mm;
  clipper::MAtom ma_p, ma_n, ma_c, ma_o;
  clipper::MAtom mao5, mac5, mac4, mac3, mao3, mac2, mac1, mao4;
  ma_p = ma_n = ma_c = ma_o = clipper::MAtom::null();
  ma_p.set_u_iso ( 0.25 ); ma_p.set_occupancy( 1.0 );
  ma_p.set_id( " P  " ); ma_p.set_element( "P" );
  ma_n.set_u_iso ( 0.25 ); ma_n.set_occupancy( 1.0 );
  ma_n.set_id( " N1 " ); ma_n.set_element( "N" );
  ma_c.set_u_iso ( 0.25 ); ma_c.set_occupancy( 1.0 );
  ma_c.set_id( " C  " ); ma_c.set_element( "C" );
  ma_o.set_u_iso ( 0.25 ); ma_o.set_occupancy( 1.0 );
  ma_o.set_id( " O  " ); ma_o.set_element( "O" );
  mao5 = mao4 = mao3 = ma_o;
  mac5 = mac4 = mac3 = mac2 = mac1 = ma_c;
  mao5.set_id( " O5'" ); mao4.set_id( " O4'" ); mao3.set_id( " O3'" ); 
  mac5.set_id( " C5'" ); mac4.set_id( " C4'" ); mac3.set_id( " C3'" ); 
  mac2.set_id( " C2'" ); mac1.set_id( " C1'" );
  ma_p.set_coord_orth( coord_p() );
  ma_n.set_coord_orth( coord_n() );
  mao5.set_coord_orth( coord_o5() );
  mao4.set_coord_orth( coord_o4() );
  mao3.set_coord_orth( coord_o3() );
  mac5.set_coord_orth( coord_c5() );
  mac4.set_coord_orth( coord_c4() );
  mac3.set_coord_orth( coord_c3() );
  mac2.set_coord_orth( coord_c2() );
  mac1.set_coord_orth( coord_c1() );
  if ( !ma_p.coord_orth().is_null() ) mm.insert( ma_p );
  if ( !mao5.coord_orth().is_null() ) mm.insert( mao5 );
  if ( !mac5.coord_orth().is_null() ) mm.insert( mac5 );
  if ( !mac4.coord_orth().is_null() ) mm.insert( mac4 );
  if ( !mao4.coord_orth().is_null() ) mm.insert( mao4 );
  if ( !mac3.coord_orth().is_null() ) mm.insert( mac3 );
  if ( !mao3.coord_orth().is_null() ) mm.insert( mao3 );
  if ( !mac2.coord_orth().is_null() ) mm.insert( mac2 );
  if ( !mac1.coord_orth().is_null() ) mm.insert( mac1 );
  if ( !ma_n.coord_orth().is_null() ) mm.insert( ma_n );
  mm.set_type( std::string( 1, typ ) );
  return mm;
}


void NucleicAcid::transform( const clipper::RTop_orth& rtop )
{
  if ( !clipper::Util::is_null( p_x ) ) {
    clipper::Coord_orth c = rtop * coord_p();
    p_x = float( c.x() );
    p_y = float( c.y() );
    p_z = float( c.z() );
  }
  if ( !clipper::Util::is_null( o5x ) ) {
    clipper::Coord_orth c = rtop * coord_o5();
    o5x = float( c.x() );
    o5y = float( c.y() );
    o5z = float( c.z() );
  }
  if ( !clipper::Util::is_null( c5x ) ) {
    clipper::Coord_orth c = rtop * coord_c5();
    c5x = float( c.x() );
    c5y = float( c.y() );
    c5z = float( c.z() );
  }
  if ( !clipper::Util::is_null( c4x ) ) {
    clipper::Coord_orth c = rtop * coord_c4();
    c4x = float( c.x() );
    c4y = float( c.y() );
    c4z = float( c.z() );
  }
  if ( !clipper::Util::is_null( o4x ) ) {
    clipper::Coord_orth c = rtop * coord_o4();
    o4x = float( c.x() );
    o4y = float( c.y() );
    o4z = float( c.z() );
  }
  if ( !clipper::Util::is_null( c3x ) ) {
    clipper::Coord_orth c = rtop * coord_c3();
    c3x = float( c.x() );
    c3y = float( c.y() );
    c3z = float( c.z() );
  }
  if ( !clipper::Util::is_null( o3x ) ) {
    clipper::Coord_orth c = rtop * coord_o3();
    o3x = float( c.x() );
    o3y = float( c.y() );
    o3z = float( c.z() );
  }
  if ( !clipper::Util::is_null( c2x ) ) {
    clipper::Coord_orth c = rtop * coord_c2();
    c2x = float( c.x() );
    c2y = float( c.y() );
    c2z = float( c.z() );
  }
  if ( !clipper::Util::is_null( c1x ) ) {
    clipper::Coord_orth c = rtop * coord_c1();
    c1x = float( c.x() );
    c1y = float( c.y() );
    c1z = float( c.z() );
  }
  if ( !clipper::Util::is_null( n_x ) ) {
    clipper::Coord_orth c = rtop * coord_n();
    n_x = float( c.x() );
    n_y = float( c.y() );
    n_z = float( c.z() );
  }
}


void NucleicAcid::set_flag()
{
  if ( !clipper::Util::is_null( c1x ) &&
       !clipper::Util::is_null( c3x ) &&
       !clipper::Util::is_null( c4x ) ) {
    if ( !clipper::Util::is_null( n_x ) &&
         !clipper::Util::is_null( p_x ) &&
         !clipper::Util::is_null( c2x ) &&
         !clipper::Util::is_null( c5x ) &&
         !clipper::Util::is_null( o3x ) &&
         !clipper::Util::is_null( o4x ) &&
         !clipper::Util::is_null( o5x ) ) {
      flg = COMPLETE;
    } else {
      flg = INCOMPLETE;
    }
  }
  else {
    flg = NONE;
  }
}


// Chain classes

bool Chain::add_pdb( const clipper::String file )
{
  const int mmdbflags = ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors | ::mmdb::MMDBF_IgnoreRemarks;
  clipper::MMDBfile mfile;
  clipper::MiniMol mol;
  mfile.SetFlag( mmdbflags );
  mfile.read_file( file );
  mfile.import_minimol( mol );
  if ( mol.size() == 0 ) return false;
  for ( int c = 0; c < mol.size(); c++ ) {
    clipper::MPolymer mp;
    // select monomers by occupancy
    for ( int r = 0; r < mol[c].size(); r++ ) {
      if ( mol[c][r].lookup( " C1'", clipper::MM::ANY ) >= 0 &&
           mol[c][r].lookup( " C2'", clipper::MM::ANY ) >= 0 &&
           mol[c][r].lookup( " C3'", clipper::MM::ANY ) >= 0 &&
           mol[c][r].lookup( " C4'", clipper::MM::ANY ) >= 0 &&
           mol[c][r].lookup( " C5'", clipper::MM::ANY ) >= 0 &&
           mol[c][r].lookup( " O3'", clipper::MM::ANY ) >= 0 &&
           mol[c][r].lookup( " O4'", clipper::MM::ANY ) >= 0 &&
           mol[c][r].lookup( " O5'", clipper::MM::ANY ) >= 0 &&
           mol[c][r].lookup( " P  ", clipper::MM::ANY ) >= 0 ) {
        int a = mol[c][r].lookup( " C4'", clipper::MM::ANY );
        if ( mol[c][r][a].occupancy() > 0.01 &&
             mol[c][r][a].u_iso() < clipper::Util::b2u(100.0) )
          mp.insert( mol[c][r] );
      }
    }
    // shift centre-of-mass of chain to the origin
    clipper::Coord_orth cm( 0.0, 0.0, 0.0 );
    double              sm = 0.0;
    for ( int r = 0; r < mp.size(); r++ ) {
      int a = mp[r].lookup( " C4'", clipper::MM::ANY );
      cm += mp[r][a].coord_orth();
      sm += 1.0;
    }
    cm = (-1.0/sm) * cm;
    clipper::RTop_orth rt( clipper::Mat33<>::identity(), cm );
    mp.transform( rt );
    // now add the chain to the db
    for ( int r = 0; r < mp.size(); r++ ) {
      NucleicAcid rp( mp[r] );
      if ( rp.flag() == NucleicAcid::COMPLETE && rp.type() != ' ' )
        add_monomer( rp );
    }
  }
  return true;
}


bool Chain::save_db( const clipper::String file ) const
{
  /*
  std::ofstream fs( file.c_str(), std::ios::out | std::ios::binary );
  char d[20];
  for ( int r = 0; r < dbmonomers.size(); r++ ) {
    dbmonomers[r].data_export( d );
    fs.write( d, 20 );
  }
  fs.close();
  */
  return true;
}


bool Chain::load_db( const clipper::String file )
{
  /*
  dbmonomers.clear();
  // read whole file (for speed)
  std::ifstream fs( file.c_str(), std::ios::in | std::ios::binary );
  if ( !fs ) return false;
  fs.seekg( 0, std::ios::end );
  int i2 = fs.tellg();
  fs.seekg( 0, std::ios::beg );
  int i1 = fs.tellg();
  int l = i2 - i1;
  char d[l];
  fs.read( d, l );
  fs.close();
  if ( l%20 != 0 ) return false;
  // import file data
  dbmonomers.resize( l/20 );
  for ( int r = 0; r < dbmonomers.size(); r++ ) {
    dbmonomers[r].data_import( d + 20*r );
  }
  */
  return true;
}


bool Chain::merge( const Chain& other, const std::vector<double>& wgt )
{
  /*
  if ( other.size() != size() ) return false;
  if ( wgt.size() != 3*size() ) return false;
  for ( int r = 0; r < dbmonomers.size(); r++ )
    dbmonomers[r].merge( other.dbmonomers[r],
                         wgt[3*r], wgt[3*r+1], wgt[3*r+2] );
  */
  return true;
}


Chain Chain::extract( int offset, int len ) const
{
  Chain dbc;
  for ( int i = 0; i < len; i++ ) dbc.add_monomer( dbmonomers[offset+i] );
  return dbc;
}


bool Chain::is_continuous() const
{
  // go through and find elements where there is a chain break
  const double dmin = 2.0;
  std::vector<bool> cterm( dbmonomers.size(), false );
  for ( int i = 0; i < dbmonomers.size()-1; i++ ) {
    int j = i + 1;
    const clipper::Coord_orth co1 = dbmonomers[i].coord_o3();
    const clipper::Coord_orth co2 = dbmonomers[j].coord_p();
    if ( co1.is_null() || co2.is_null() ) return false;
    const double d2 = ( co1 - co2 ).lengthsq();
    if ( d2 > dmin*dmin ) return false;
  }
  return true;
}


void Chain::transform( const clipper::RTop_orth& rtop )
{
  for ( int r = 0; r < dbmonomers.size(); r++ )
    dbmonomers[r].transform( rtop );
}


}
