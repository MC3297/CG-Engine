#pragma once

#include "types.hpp"
#include "color.hpp"

//implementation of bresenham's line algorithm
//includes z component to account for z-buffer
void draw_line(vector<vector<color>>& plot, vector<vector<db>>& zb, int x0, int y0, int z0, int x1, int y1, int z1, color p);