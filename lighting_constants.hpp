#pragma once

#include "types.hpp"
#include "color.hpp"
#include "vector_math.hpp"

color camb(50,50,50);
vector<color> clight = {color(0,255,255)};
vector<Vec3> light = {{0.5,0.75,1}};
Vec3 view = {0,0,1};
Vec3 kamb = {0.2,0.2,0.2};
Vec3 kdiff = {0.6,0.6,0.6};
Vec3 kspec = {0.82,0.82,0.82};