//
// Created by Jordan Dialpuri on 16/02/2024.
//
#include <pybind11/pybind11.h>

#include "nautilus-include.h"
#include "nautilus-run.h"

namespace py = pybind11;

using namespace pybind11::literals;

PYBIND11_MODULE(nautilus_module, m) {
py::class_<NautilusInput>(m, "Input")
            .def(py::init< const std::string&, // mtzin
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
             "colin_hl"_a, "colin_phifom"_a, "colin_fc"_a, "colinfree"_a, ""
             );

py::class_<NautilusOutput>(m, "Output")
            .def(py::init< const std::string& >()); // pdbout

m.def("run", &run, "input"_a, "output"_a, "cycles"_a, "Run nucleofind-build");
}