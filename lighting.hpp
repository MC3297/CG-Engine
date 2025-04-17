#pragma once

#include "types.hpp"
#include "color.hpp"
#include "vector_math.hpp"

color get_amb(color camb, Vec3 kamb);

color get_diff(color clight, Vec3 light, Vec3 kdiff, Vec3 norm);

color get_spec(color clight, Vec3 light, Vec3 view, Vec3 kspec, Vec3 norm);

color limit(color c);

color get_lighting(color camb, vector<color> clight, vector<Vec3> light, Vec3 view, Vec3 norm, Vec3 kamb, Vec3 kdiff, Vec3 kspec);