//
// Created by Jordan Dialpuri on 13/09/2023.
//
#include <string>
#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-ccp4.h>

#ifndef NAUTILUS_NAUTILUS_PREDICT_H
#define NAUTILUS_NAUTILUS_PREDICT_H

class Predictions {
public:
    Predictions(const std::string& phosphate_path) {
        load_predicted_phosphate_map(phosphate_path);
    }

    void load_predicted_phosphate_map(const std::string& phosphate_path) {
        clipper::CCP4MAPfile file;
        file.open_read(phosphate_path);
        file.import_xmap(phosphate_map);
        file.close_read();
    }

    clipper::Xmap<float> phosphate_map;
};

#endif //NAUTILUS_NAUTILUS_PREDICT_H
