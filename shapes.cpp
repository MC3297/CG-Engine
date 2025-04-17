#include "shapes.hpp"

void add_box(matrix& m, db x, db y, db z, db w, db h, db d) {
    //top
    m.add_poly({x,y,z},{x+w,y,z},{x+w,y,z-d});
    m.add_poly({x+w,y,z-d},{x,y,z-d},{x,y,z});
    //bottom
    m.add_poly({x+w,y-h,z-d},{x+w,y-h,z},{x,y-h,z});
    m.add_poly({x,y-h,z},{x,y-h,z-d},{x+w,y-h,z-d});
    //left
    m.add_poly({x,y,z},{x,y,z-d},{x,y-h,z-d});
    m.add_poly({x,y-h,z-d},{x,y-h,z},{x,y,z});
    //right
    m.add_poly({x+w,y-h,z-d},{x+w,y,z-d},{x+w,y,z});
    m.add_poly({x+w,y,z},{x+w,y-h,z},{x+w,y-h,z-d});
    //front
    m.add_poly({x+w,y-h,z},{x+w,y,z},{x,y,z});
    m.add_poly({x,y,z},{x,y-h,z},{x+w,y-h,z});
    //back
    m.add_poly({x,y,z-d},{x+w,y,z-d},{x+w,y-h,z-d});
    m.add_poly({x+w,y-h,z-d},{x,y-h,z-d},{x,y,z-d});
}

vector<vector<vector<db>>> gen_sphere(db cx, db cy, db cz, db r, int step) {
    vector<vector<vector<db>>> ans(step+1, vector<vector<db>>(step+1));
    for (int rot = 0; rot <= step; rot++) {
        for (int circ = 0; circ <= step; circ++) {
            db x = r*cos(M_PI*circ/(db)step) + cx;
            db y = r*sin(M_PI*circ/(db)step)*cos(2*M_PI*rot/(db)step)+cy;
            db z = r*sin(M_PI*circ/(db)step)*sin(2*M_PI*rot/(db)step)+cz;
            ans[rot][circ] = {x,y,z};
        }
    }
    return ans;
}

void add_sphere(matrix& m, db cx, db cy, db cz, db r, const int step) {
    vector<vector<vector<db>>> pts = gen_sphere(cx,cy,cz,r,step);
    for (int rot = 0; rot < (int)pts.size()-1; rot++) {
        for (int circ = 0; circ < (int)pts.size()-1; circ++) {
            if (circ != 0) m.add_poly(pts[rot+1][circ+1], pts[rot+1][circ], pts[rot][circ]);
            if (circ != step-1) m.add_poly(pts[rot][circ], pts[rot][circ+1], pts[rot+1][circ+1]);
        }
    }
}

vector<vector<vector<db>>> gen_torus(db cx, db cy, db cz, db r1, db r2, int step) {
    vector<vector<vector<db>>> ans(step+1, vector<vector<db>>(step+1));
    for (int rot = 0; rot <= step; rot++) {//step+1
        for (int circ = 0; circ <= step; circ++) {
            db x = cos(2*M_PI*rot/(db)step)*(r1*cos(2*M_PI*circ/(db)step)+r2)+cx;
            db y = r1*sin(2*M_PI*circ/(db)step)+cy;
            db z = -sin(2*M_PI*rot/(db)step)*(r1*cos(2*M_PI*circ/(db)step)+r2)+cz;
            ans[rot][circ] = {x,y,z};
        }
    }
    return ans;
}

void add_torus(matrix& m, db cx, db cy, db cz, db r1, db r2, int step) {
    auto pts = gen_torus(cx, cy, cz, r1, r2, step);
    for (int rot = 0; rot < (int)pts.size()-1; rot++) {
        for (int circ = 0; circ < (int)pts[0].size()-1; circ++) {
            m.add_poly(pts[rot][circ], pts[rot+1][circ], pts[rot+1][circ+1]);
            m.add_poly(pts[rot+1][circ+1], pts[rot][circ+1], pts[rot][circ]);
        }
    }
}