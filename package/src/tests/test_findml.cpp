// Copyright (c) 2024 Jordan Dialpuri, Jon Agirre, Kevin Cowtan, Paul Bond and University of York. All rights reserved

//
// Created by jordan on 29/05/24.
//
#include <gtest/gtest.h>
#include "../cpp/nautilus-findml.h"

TEST(FindMLTests, PredictedMapTests) {
    PredictedMaps map = {clipper::Xmap<float>(), clipper::Xmap<float>(), clipper::Xmap<float>()};
    EXPECT_FALSE(map.get_base_map().has_value());
}