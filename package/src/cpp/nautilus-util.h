/*! \file nautilus-util.h nautilus library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#ifndef NAUTILUS_UTIL_H
#define NAUTILUS_UTIL_H

#include <clipper/clipper-minimol.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/core/map_interp.h>

class NautilusUtil {
 public:
  static void set_reference( clipper::String& pdb );
    static void save_minimol(clipper::MiniMol& mol, const std::string& path) {
        clipper::MMDBfile mfile;
        mfile.export_minimol(mol);
        mfile.write_file(path);
    }

    static clipper::Coord_orth div(clipper::Coord_orth& c1, clipper::Coord_orth& c2) {
        return {c1.x()/c2.x(), c1.y()/c2.y(), c1.z()/c2.z()};
    }

    static clipper::Vec3<> coord_to_vec(const clipper::Coord_orth& coord) {
        return {coord.x(), coord.y(), coord.z()};
    }

    static std::vector<double> coord_to_vector(const clipper::Coord_orth& coord) {
        return {coord.x(), coord.y(), coord.z()};
    }

    static clipper::Vec3<float> coord_to_vec_float(const clipper::Coord_orth& coord) {
        return {(float)coord.x(), (float)coord.y(), (float)coord.z()};
    }

    static clipper::Coord_orth vec_to_coord(clipper::Vec3<float>& v) {
        return {v[0], v[1], v[2]};
    }

    static clipper::Coord_orth vec_to_coord(clipper::Vec3<>& v) {
        return {v[0], v[1], v[2]};
    }

    static clipper::MAtom create_atom(const clipper::Coord_orth& coord, int id = 0, std::string element = "C") {
        clipper::MAtom m_atom;
        m_atom.set_u_iso ( 0.25 );
        m_atom.set_occupancy( 1.0 );
        m_atom.set_id(std::to_string(id));
        m_atom.set_element( std::move(element));
        m_atom.set_coord_orth(std::move(coord));
        return m_atom;
    }

    static clipper::MAtom create_atom(const clipper::Coord_orth& coord, std::string id, std::string element = "C") {
        clipper::MAtom m_atom;
        m_atom.set_u_iso ( 0.25 );
        m_atom.set_u_aniso_orth(clipper::U_aniso_orth(clipper::Mat33sym<>::null()));
        m_atom.set_occupancy( 1.0 );
        m_atom.set_id(std::move(id));
        m_atom.set_element( std::move(element));
        m_atom.set_coord_orth(coord);
        return m_atom;
    }

    static clipper::MAtom create_atom(const clipper::Vec3<>& coord, int id = 0, std::string element = "C") {
        clipper::MAtom m_atom;
//        m_atom.set_u_iso ( 0.25 );
        m_atom.set_occupancy( 1.0 );
        m_atom.set_id(id);
        m_atom.set_element( std::move(element));
        m_atom.set_coord_orth(NautilusUtil::vec_to_coord((clipper::Vec3<float> &) coord));
        return m_atom;
    }

    static clipper::MiniMol combine_minimols(clipper::MiniMol& mol1, clipper::MiniMol& mol2) {
        clipper::MiniMol return_mol = {mol1.spacegroup(), mol1.cell()};

        clipper::MModel return_model;
        for (int i = 0; i < mol1.model().size(); i++) {
//            mol1[i].set_id(i);
            return_model.insert(mol1[i]);
        }
        for (int i = 0; i < mol2.model().size(); i++) {
//            mol1[i].set_id(2);
            return_model.insert(mol2[i]);
        }

        return_mol.model() = return_model;
        return return_mol;
    }

    static clipper::MMonomer del_atom(clipper::MMonomer& mon, int atom_index) {

        clipper::MMonomer return_monomer;
        return_monomer.set_id(mon.id());
        return_monomer.set_type(mon.type());
        return_monomer.set_seqnum(mon.seqnum());
        for (int i = 0; i < mon.size(); i++) {
            if (i != atom_index){
                return_monomer.insert(mon[i]);
            }
        }
        return return_monomer;
    }

    static void save_monomer(clipper::MMonomer& mon, const std::string& path) {
        clipper::MPolymer mp;
        mp.set_id(1);
        mp.insert(mon);

        clipper::MModel model;
        model.insert(mp);

        clipper::MiniMol mol;
        mol.model() = model;

        clipper::MMDBfile out_file;
        out_file.export_minimol(mol);
        out_file.write_file(path);
    }

    static std::vector<clipper::Coord_orth> atom_index_to_goal(clipper::MiniMol& mol, std::vector<clipper::MAtomIndexSymmetry>& atom_list) {
        std::vector<clipper::Coord_orth> return_list;

        for (auto& atom: atom_list) {
            clipper::Coord_orth grid = mol[atom.polymer()][atom.monomer()][atom.atom()].coord_orth();
            return_list.emplace_back(grid);
        }

        return return_list;
    }

    static std::vector<clipper::Coord_orth> mol_to_goal(clipper::MiniMol& mol) {
        std::vector<clipper::Coord_orth> return_list;

        for (int poly = 0; poly < mol.model().size(); poly++) {
            for (int mon = 0; mon < mol.model()[poly].size(); mon++) {
                for (int atom = 0; atom < mol.model()[poly][mon].size(); atom++) {
                    clipper::Coord_orth grid = mol[poly][mon][atom].coord_orth();
                    return_list.emplace_back(grid);
                }
            }
        }

        return return_list;

    }


    static std::vector<clipper::Coord_grid> orth_to_grid(clipper::Xmap<float>& xmap, std::vector<clipper::Coord_orth>& grid) {
        std::vector<clipper::Coord_grid> return_list;

        return_list.reserve(grid.size());
        for (auto& point: grid) {
            return_list.emplace_back(point.coord_frac(xmap.cell()).coord_grid(xmap.grid_sampling()));
        }

        return return_list;
    }


    static clipper::Coord_grid calculate_mean(std::vector<clipper::Coord_grid>& grid_points) {

        clipper::Coord_grid centroid_sum = {0,0,0};

        for (auto& point: grid_points) {
            centroid_sum+=point;
        }

        clipper::Coord_grid centroid = {
                centroid_sum.u() / (int)grid_points.size(),
                centroid_sum.v() / (int)grid_points.size(),
                centroid_sum.w() / (int)grid_points.size()
        };

        return centroid;
    }

    static bool compare_grid(const clipper::Coord_grid& g1, const clipper::Coord_grid& g2) {
        return (g1.u() == g2.u() && g1.v() == g2.v() && g1.w() == g2.w());
    }
};


class NautilusLog {
 public:
  NautilusLog( clipper::String& title ) : title_(title), currentcpu(0.0) { log(""); }
  void log( const clipper::String& id );
  void log( const clipper::String& id, const clipper::MiniMol& mol, bool view );
  clipper::String log_info( const clipper::MiniMol& mol, bool summary );
  void xml( const clipper::String& file ) const; //, const clipper::MiniMol& mol ); edited SWH
  void profile();
 private:
  struct cycdat { int nchns, nseq, nres, nmax; }; // added by SWH
  std::vector<cycdat> data; //added by SWH
  std::vector<std::pair<std::string,double> > prof;
  clipper::String title_;
  double currentcpu;
};





#endif
