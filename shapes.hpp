#pragma once

#include "types.hpp"
#include "matrix.hpp"

void add_box(matrix& m, db x, db y, db z, db w, db h, db d);

vector<vector<vector<db>>> gen_sphere(db cx, db cy, db cz, db r, int step);

void add_sphere(matrix& m, db cx, db cy, db cz, db r, const int step);

vector<vector<vector<db>>> gen_torus(db cx, db cy, db cz, db r1, db r2, int step);

void add_torus(matrix& m, db cx, db cy, db cz, db r1, db r2, int step);