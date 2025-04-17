#pragma once

#include "types.hpp"

struct color {
    int r,g,b;
    color(int x, int y, int z);
    color();
    string to_str();
};