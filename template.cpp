#include "types.h"
#include "lighting_constants.h"
#include "color.h"
#include "vector_math.h"
#include "bresenham_line.h"
#include "matrix.h"
#include "matrix_transform.h"
#include "lighting.h"


#define INCR 0.0005

void circle(matrix& m, db cx, db cy, db cz, db r) {
    for (db t = INCR; t <= 1+INCR; t += INCR) {
        m.add_edge(r*cos(2*M_PI*t)+cx, r*sin(2*M_PI*t)+cy, cz, r*cos(2*M_PI*(t-INCR))+cx, r*sin(2*M_PI*(t-INCR))+cy, cz);
    }
}

db polynomial(array<db,4> pts, db x, db y, bool c=1) {
    //c is true if we want binomial coefficients
    //pts[0] * x^n + pts[1] * x^(n-1) * y + pts[2] * x^(n-2) * y^2 ...
    db a = 1, b = 1, coeff = 1, ans = 0;
    int n = (int)pts.size();
    for (int i = 0; i < n; i++) {
        pts[n-i-1] *= a;
        pts[i] *= b;
        a *= x;
        b *= y;
        if (c) {
            pts[i] *= coeff;
            coeff *= (db)(n-i-1)/(i+1);
        }
    }
    for (int i = 0; i < n; i++) 
        ans += pts[i];
    return ans;
}

void bezier_curve(matrix& m, db x0, db y0, db x1, db y1, db x2, db y2, db x3, db y3) {
    for (db t = INCR; t <= 1+INCR; t += INCR) {
        db xcurr = polynomial({x0,x1,x2,x3}, 1-t, t);
        db ycurr = polynomial({y0,y1,y2,y3}, 1-t, t);
        db xprev = polynomial({x0,x1,x2,x3}, 1-t+INCR, t-INCR);
        db yprev = polynomial({y0,y1,y2,y3}, 1-t+INCR, t-INCR);
        m.add_edge(xprev, yprev, 0, xcurr, ycurr, 0);
    }
}

void hermite_curve(matrix& m, db x0, db y0, db x1, db y1, db rx0, db ry0, db rx1, db ry1) {
    matrix h_inv({{2,-3,0,1},{-2,3,0,0},{1,-2,1,0},{1,-1,0,0}}),
    gx({{x0,x1,rx0,rx1}}), gy({{y0,y1,ry0,ry1}}),
    xcoeff = h_inv*gx, ycoeff = h_inv*gy;
    for (db t = INCR; t <= 1+INCR; t += INCR) {
        db xcurr = polynomial(xcoeff.ma[0], t, 1, 0);
        db ycurr = polynomial(ycoeff.ma[0], t, 1, 0);
        db xprev = polynomial(xcoeff.ma[0], t+INCR, 1, 0);
        db yprev = polynomial(ycoeff.ma[0], t+INCR, 1, 0);
        m.add_edge(xprev, yprev, 0, xcurr, ycurr, 0);
    }
}

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

/*vector<vector<vector<db>>> gen_bezier(vector<db> vx, vector<db> vy, int step) {
    vector<vector<vector<db>>> ans(step+1);
    matrix m;
    for (int t = 0; t <= step; t++) {
        db xcurr = polynomial(vx, 1-(db)t/step, (db)t/step);
        db ycurr = polynomial(vy, 1-(db)t/step, (db)t/step);
        m.add_point(xcurr, ycurr, 0);
    }
    for (int rot = 0; rot <= step; rot++) {//step+1
        ans[rot] = (mk_rotY(rot*-360/step)*m).ma;
    }
    return ans;
}

void add_bezier_shape(matrix& m, array<db,4> vx, array<db,4> vy, int step) {
    auto pts = gen_bezier(vx, vy, step);
    for (int rot = 0; rot < (int)pts.size()-1; rot++) {
        for (int circ = 0; circ < (int)pts[0].size()-1; circ++) {
            //m.add_poly(pts[rot][circ], pts[rot+1][circ], pts[rot+1][circ+1]);
            //m.add_poly(pts[rot+1][circ+1], pts[rot][circ+1], pts[rot][circ]);
            m.add_poly(pts[rot+1][circ+1], pts[rot+1][circ], pts[rot][circ]);
            m.add_poly(pts[rot][circ], pts[rot][circ+1], pts[rot+1][circ+1]);
        }
    }
}*/