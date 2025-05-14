//
// Created by Jordan Dialpuri on 01/05/2025.
//

#ifndef IO_H
#define IO_H

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/core/hkl_compute.h>

#include <string>
#include <optional>

#include "predicted-maps.h"

namespace NucleoFind::IO {
    class Input {
    public:
        Input(
            const std::string &mtzin,
            const std::string &seqin,
            const std::string &pdbin,
            const std::string &phospredin,
            const std::string &sugarpredin,
            const std::string &basepredin,
            const std::string &colin_fo,
            const std::string &colin_hl,
            const std::string &colin_phifom,
            const std::string &colin_fc,
            const std::string &colin_free,
            const bool& em,
            const std::string& database
        ) {
            this->mtzin = mtzin;
            this->seqin = seqin;
            this->pdbin = pdbin;
            this->phospredin = phospredin;
            this->sugarpredin = sugarpredin;
            this->basepredin = basepredin;

            this->colin_fo = colin_fo;
            this->colin_hl = colin_hl;
            this->colin_phifom = colin_phifom;
            this->colin_fc = colin_fc;
            this->colin_free = colin_free;
            this->em = em;
            this->database = database;

            if (mtzin.empty()) { throw std::runtime_error("MTZ path must not be empty"); }
            //        if (seqin.empty()) { throw std::runtime_error("SEQ Path must not be empty");}
            // if (pdbin.empty()) { throw std::runtime_error("PDB Path must not be empty");}
            if (phospredin.empty()) { throw std::runtime_error("Phosphate prediction path must not be empty"); }
        };

        [[nodiscard]] std::string get_mtz_path() const { return mtzin; }
        [[nodiscard]] std::string get_seq_path() const { return seqin; }

        [[nodiscard]] std::optional<std::string> get_pdb_path() const {
            if (pdbin.empty()) return std::nullopt;
            return pdbin;
        }

        [[nodiscard]] std::optional<std::string> get_phosphate_prediction_path() const {
            if (phospredin.empty()) return std::nullopt;
            return phospredin;
        }

        [[nodiscard]] std::optional<std::string> get_sugar_prediction_path() const {
            if (sugarpredin.empty()) return std::nullopt;
            return sugarpredin;
        }

        [[nodiscard]] std::optional<std::string> get_base_prediction_path() const {
            if (basepredin.empty()) return std::nullopt;
            return basepredin;
        }

        [[nodiscard]] std::optional<std::string> get_fobs() const {
            if (colin_fo.empty()) {
                return std::nullopt;
            }
            return colin_fo;
        }

        [[nodiscard]] std::optional<std::string> get_hl() const {
            if (colin_hl.empty()) {
                return std::nullopt;
            }
            return colin_hl;
        }

        [[nodiscard]] std::optional<std::string> get_phifom() const {
            if (colin_phifom.empty()) {
                return std::nullopt;
            }
            return colin_phifom;
        }

        [[nodiscard]] std::optional<std::string> get_fc() const {
            if (colin_fc.empty()) {
                return std::nullopt;
            }
            return colin_fc;
        }

        [[nodiscard]] std::optional<std::string> get_free() const {
            if (colin_free.empty()) {
                return std::nullopt;
            }
            return colin_free;
        }

        [[nodiscard]] bool is_em() const {
            return em;
        }

        [[nodiscard]] std::optional<std::string> get_database() const {
            if (database.empty()) {
                return std::nullopt;
            }
            return database;
        }

    private:
        std::string mtzin;
        std::string seqin;
        std::string pdbin;
        std::string phospredin;
        std::string sugarpredin;
        std::string basepredin;
        std::string colin_fo;
        std::string colin_hl;
        std::string colin_phifom;
        std::string colin_fc;
        std::string colin_free;
        bool em;
        std::string database;
    };

    class Output {
    public:
        Output(const std::string &pdbout, const std::string &xmlout) {
            this->pdbout = pdbout;
            this->xmlout = xmlout;
            if (pdbout == "") { throw std::runtime_error("PDB Out must not be empty"); }
        };

        [[ nodiscard ]] std::string get_pdb_out() const { return pdbout; }

        [[ nodiscard ]] std::optional<std::string> get_xml_out() const {
            if (xmlout.empty()) {
                return std::nullopt;
            }
            return xmlout;
        }

    private:
        std::string pdbout;
        std::string xmlout;
    };

    CCP4Program initialise_ccp4_program(const std::string &version);

    class MTZ {
    public:
        typedef clipper::HKL_data_base::HKL_reference_index HRI;

        explicit MTZ(const Input &input, double res_in = 2.0): input(input), res_in(res_in) {
            load_work_map();
            print_stats();
        };

        void load_work_map();

        void print_stats();

        [[nodiscard]] clipper::Xmap<float> get_xmap() const {
            return xwrk;
        }

    private:
        Input input;
        double res_in;
        clipper::HKL_data<clipper::data32::F_phi> fphi;
        clipper::HKL_data<clipper::data32::F_sigF> wrk_f1;
        clipper::HKL_info hkls;
        clipper::Xmap<float> xwrk;
    };

    class Predictions {
    public:
        explicit Predictions(const Input &input): input(input) {
            load_predicted_maps();
        }

        void load_predicted_maps();

        [[nodiscard]] PredictedMaps get_predictions() const {
            return {xphospred, xsugarpred, xbasepred};
        }
    private:
        Input input;
        clipper::Xmap<float> xphospred;
        clipper::Xmap<float> xsugarpred;
        clipper::Xmap<float> xbasepred;

    };

    class Sequence {
    public:
        explicit Sequence(const Input &input): input(input) {
            load_sequence();
        }

        void load_sequence();

        [[nodiscard]] clipper::MMoleculeSequence get_sequence() const {
            return seq_wrk;
        }

    private:
        Input input;
        clipper::MMoleculeSequence seq_wrk;

    };

    class Model {
    public:
        explicit Model(const Input &input, const clipper::Spacegroup& spg, const clipper::Cell &cell): input(input) {
            mol_wrk = clipper::MiniMol(spg, cell);
            load_model();
        }

        void load_model() ;

        [[nodiscard]] clipper::MiniMol get_model() const {
            return mol_wrk;
        }

    private:
        Input input;
        clipper::MiniMol mol_wrk;


    };;

};


#endif //IO_H
