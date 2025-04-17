#pragma once

#include "types.hpp"
#include "color.hpp"
#include "vector_math.hpp"
#include "lighting.hpp"
#include "bresenham_line.hpp"

db sqr(db a);
db distsqred(Point a, Point b);

struct matrix {
    int vert;
    int horiz;
    vector<Point> ma;
    matrix(int horizontal, int vertical);
    matrix(vector<Point> m);
    matrix();
    void clear();
    void ident();
    Point& operator[](int ind);
    matrix operator*(const matrix& x);
    void print();
    void add_point(db x, db y, db z);
    void add_edge(db x0, db y0, db z0, db x1, db y1, db z1);
    void draw_lines(vector<vector<color>>& plot, vector<vector<db>>& zb, color p);

    void add_poly(vector<db> a, vector<db> b, vector<db> c);

    void draw_poly(vector<vector<color>>& plot, vector<vector<db>>& zb, color camb, vector<color> clight, vector<Vec3> light, Vec3 view, Vec3 kamb, Vec3 kdiff, Vec3 kspec);

    void draw_scanline(db x0, db x1, db z0, db z1, db y, vector<vector<color>>& plot, vector<vector<db>>& zb, color c);

    void scanline_convert(int ind, vector<vector<color>>& plot, vector<vector<db>>& zb, color clr);
};