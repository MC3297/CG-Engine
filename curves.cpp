#include "curves.hpp"

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