#include "color.hpp"

color::color(int x, int y, int z) : r(x), g(y), b(z) {};

color::color() : r(0), g(0), b(0) {}

string color::to_str() {
    return to_string(r)+" "+to_string(g)+" "+to_string(b)+" ";
}