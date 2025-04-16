#pragma once

#include "types.h"
using Vec3 = array<db,3>;

//vector math functions

//normalizes vector
void normalize(Vec3& v) {
    db mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] /= mag;
    v[1] /= mag;
    v[2] /= mag;
}

//returns dot product of a and b
db dot_prod(Vec3 a, Vec3 b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

//returns cross product (s-r)x(t-r)
Vec3 normal_surface(Point r, Point s, Point t) {
    Vec3 a = {s[0]-r[0], s[1]-r[1], s[2]-r[2]};
    Vec3 b = {t[0]-r[0], t[1]-r[1], t[2]-r[2]};
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}