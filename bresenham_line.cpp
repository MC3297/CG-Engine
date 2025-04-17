#include "bresenham_line.hpp"

void draw_line(vector<vector<color>>& plot, vector<vector<db>>& zb, int x0, int y0, int z0, int x1, int y1, int z1, color p) {
    int dx = abs(x1-x0);
    int dy = abs(y1-y0);
    int d = 0;
    if (abs(x1-x0) > abs(y1-y0)) {
        if (x0 > x1) {
            swap(x0, x1);
            swap(y0, y1);
        }
        int y = y0;
        int incr = (y0 < y1)? 1:-1;
        db z = z0;
        for (int x = x0; x <= x1; x++) {
            if (0 <= x && x < (int)plot.size() && 0 <= y && y < (int)plot[0].size()) {
                if (z > zb[x][y]) {
                    plot[x][y] = p;
                    zb[x][y] = z;
                }
            }
            d += (dy<<1) - dx;
            if (d > 0) {
                d -= dx;
                y += incr;
            }
            else d += dx;
            z = z0 + (x-x0+1)*(db)(z1-z0)/(x1-x0);
            //z += (db)dz/dx;
        }
    }
    else {
        if (y0 > y1) {
            swap(x0,x1);
            swap(y0,y1);
        }
        int x = x0;
        int incr = (x0 < x1)? 1:-1, dz = z1-z0;
        db z = z0;
        for (int y = y0; y <= y1; y++) {
            if (0 <= x && x < (int)plot.size() && 0 <= y && y < (int)plot[0].size()) {
                if (z > zb[x][y]) {
                    plot[x][y] = p;
                    zb[x][y] = z;
                }
            }
            d += dy - (dx<<1);
            if (d > 0) d -= dy;
            else {
                d += dy;
                x += incr;
            }
            z = z0 + (y-y0+1)*(db)(z1-z0)/(y1-y0);
            //z += (db)dz/dy;
        }
    }
}