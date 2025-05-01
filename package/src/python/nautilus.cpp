//
// Created by Jordan Dialpuri on 16/02/2024.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include "../cpp/nucleofind/input.h"
#include "nautilus-run.h"

namespace nb = nanobind;

using namespace nb::literals;

NB_MODULE(nautilus_module, m) {
nb::class_<NucleoFind::IO::Input>(m, "Input")
            .def(nb::init< const std::string&, // mtzin
             const std::string&, // seqin
             const std::string&, // pdbin
             const std::string&, // phospredin
             const std::string&, // sugarpredin
             const std::string&, // basepredin
             const std::string&, // colin_fo
             const std::string&, // colin_hl
             const std::string&, // colin_phifom
             const std::string&,  // colin_fc
             const std::string&  >(), // colin_free
             "mtzin"_a, "seqin"_a, "pdbin"_a, "phospredin"_a, "sugarpredin"_a, "basepredin"_a, "colin_fo"_a,
             "colin_hl"_a, "colin_phifom"_a, "colin_fc"_a, "colinfree"_a, "Input to NucleoFind"
             );

nb::class_<NucleoFind::IO::Output>(m, "Output")
            .def(nb::init< const std::string&, const std::string&  >()); // pdbout, xmlout

m.def("run", &run, "input"_a, "output"_a, "cycles"_a, "Run nucleofind-build");
}