#pragma once

#include "types.hpp"
using Vec3 = array<db,3>;

//vector math functions

//normalizes vector
void normalize(Vec3& v);

//returns dot product of a and b
db dot_prod(Vec3 a, Vec3 b);

//returns cross product (s-r)x(t-r)
Vec3 normal_surface(Point r, Point s, Point t);