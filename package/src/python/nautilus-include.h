#include <string>
#include <optional>

class NautilusInput { 
public:
    NautilusInput(
        const std::string& mtzin,
        const std::string& seqin,
        const std::string& pdbin,
        const std::string& phospredin,
        const std::string& sugarpredin,
        const std::string& basepredin,
        const std::string& colin_fo,
        const std::string& colin_hl, 
        const std::string& colin_phifom,
        const std::string& colin_fc,
        const std::string& colin_free
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

        if (mtzin.empty()) { throw std::runtime_error("MTZ Path must not be empty");}
        if (seqin.empty()) { throw std::runtime_error("SEQ Path must not be empty");}
        if (pdbin.empty()) { throw std::runtime_error("PDB Path must not be empty");}
        if (phospredin.empty()) { throw std::runtime_error("Phosphate predicition Path must not be empty");}

    };

    [[nodiscard]] std::string get_mtz_path() const { return mtzin; }
    [[nodiscard]] std::string get_seq_path() const { return seqin; }
    [[nodiscard]] std::string get_pdb_path() const { return pdbin; }
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
};

class NautilusOutput { 
public: 
    NautilusOutput(const std::string& pdbout) { 
        this->pdbout = pdbout; 

        if (pdbout == "") {throw std::runtime_error("PDB Out must not be empty");}
    }; 

    std::string get_pdb_out() const { return pdbout; }
private:
    std::string pdbout; 
};