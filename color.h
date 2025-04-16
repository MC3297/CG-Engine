#pragma once

#include "types.h"

struct color {
    int r,g,b;
    color(int x, int y, int z) {
        r = x, g = y, b = z;
    }
    color() {
        r = g = b = 0;
    }
    string to_str() {
        return to_string(r)+" "+to_string(g)+" "+to_string(b)+" ";
    }
};