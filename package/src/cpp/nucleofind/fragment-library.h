//
// Created by Jordan Dialpuri on 23/04/2025.
//

#ifndef LIBRARY_H
#define LIBRARY_H
#include <string>
#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>

namespace NucleoFind {

    struct TriNucleotide {
        TriNucleotide(clipper::MMonomer& m1, clipper::MMonomer& m2, clipper::MMonomer& m3): m1(m1), m2(m2), m3(m3) {
            setup(m1, m2, m3);
        }

        std::vector<clipper::Coord_orth> get_phosphates() {
            setup(m1, m2, m3);
            return {P1, P2, P3};
        }

        std::vector<clipper::MMonomer> transform(clipper::RTop_orth& rtop) {
            clipper::MMonomer m1_ = m1;
            clipper::MMonomer m2_ = m2;
            clipper::MMonomer m3_ = m3;
            m1_.transform(rtop);
            m2_.transform(rtop);
            m3_.transform(rtop);
            return {m1_, m2_, m3_};
        }

    private:
        void setup(clipper::MMonomer& m1, clipper::MMonomer& m2, clipper::MMonomer& m3);
        clipper::MMonomer m1, m2, m3;
        clipper::Coord_orth P1, P2, P3;
    };


    class TriNucleotideLibrary {
    public:
        explicit TriNucleotideLibrary() = default;
        explicit TriNucleotideLibrary(const std::string& library_path) {
            add_library(library_path);
            // std::cout << "Created a library with " << library.size() << " trinucleotides" << std::endl;
        }

        size_t size() const { return library.size(); }

        TriNucleotide operator[](int index) {
            return library[index];
        }

    private:
        void add_library(const std::string& library_path);

        std::vector<TriNucleotide> library;

    };
};

#endif //LIBRARY_H
