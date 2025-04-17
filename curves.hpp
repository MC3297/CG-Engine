#pragma once

#include "types.hpp"
#include "matrix.hpp"

void circle(matrix& m, db cx, db cy, db cz, db r);

db polynomial(array<db,4> pts, db x, db y, bool c);

void bezier_curve(matrix& m, db x0, db y0, db x1, db y1, db x2, db y2, db x3, db y3);

void hermite_curve(matrix& m, db x0, db y0, db x1, db y1, db rx0, db ry0, db rx1, db ry1);