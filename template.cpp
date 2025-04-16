#include "types.h"
#include "lighting_constants.h"
#include "color.h"


void draw_line(vector<vector<color>>& plot, vector<vector<db>>& zb, int x0, int y0, int z0, int x1, int y1, int z1, color p) {
    int dx = abs(x1-x0), dy = abs(y1-y0), d = 0;
    if (abs(x1-x0) > abs(y1-y0)) {
        if (x0 > x1) {
            int temp = x0;
            x0 = x1;
            x1 = temp;
            temp = y0;
            y0 = y1;
            y1 = temp;
            temp = z0;
            z0 = z1;
            z1 = temp;
        }
        int y = y0, incr = (y0 < y1)? 1:-1;
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
            int temp = x0;
            x0 = x1;
            x1 = temp;
            temp = y0;
            y0 = y1;
            y1 = temp;
            temp = z0;
            z0 = z1;
            z1 = temp;
        }
        int x = x0, incr = (x0 < x1)? 1:-1, dz = z1-z0;
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

//vector math
void normalize(vector<db>& v) {
    db mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] /= mag;
    v[1] /= mag;
    v[2] /= mag;
}
db dot_prod(vector<db> a, vector<db> b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
vector<db> normal_surface(vector<db> r, vector<db> s, vector<db> t) {
    vector<db> a{s[0]-r[0], s[1]-r[1], s[2]-r[2]},
    b{t[0]-r[0], t[1]-r[1], t[2]-r[2]};
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}


color get_amb(color camb, vector<db> kamb) {
    return color(camb.r*kamb[0], camb.g*kamb[1], camb.b*kamb[2]);
}

color get_diff(color clight, vector<db> light, vector<db> kdiff, vector<db> norm) {
    db cosin = dot_prod(norm,light);
    return color(clight.r*kdiff[0]*cosin, clight.g*kdiff[1]*cosin, clight.b*kdiff[2]*cosin);
}

color get_spec(color clight, vector<db> light, vector<db> view, vector<db> kspec, vector<db> norm) {
    db cf = 2*dot_prod(norm,light);
    vector<db> cfn = {norm[0]*cf, norm[1]*cf, norm[2]*cf};
    vector<db> subtr = {cfn[0]-light[0], cfn[1]-light[1], cfn[2]-light[2]};
    db cosin = pow(dot_prod(subtr,view),7);
    return color(clight.r*kspec[0]*cosin, clight.g*kspec[1]*cosin, clight.b*kspec[2]*cosin);
}

color limit(color c) {
    return color(min(max(0,c.r),255), min(max(0,c.g),255), min(max(0,c.b),255));
}

color get_lighting(color camb, vector<color> clight, vector<vector<db>> light, vector<db> view, vector<db> norm, vector<db> kamb, vector<db> kdiff, vector<db> kspec) {
    for (auto& l: light) normalize(l);
    normalize(view), normalize(norm);
    color ans;
    color amb = limit(get_amb(camb, kamb));
    ans = amb;
    for (int i = 0; i < (int)light.size(); i++) {
        color diff = limit(get_diff(clight[i], light[i], kdiff, norm));
        color spec = limit(get_spec(clight[i], light[i], view, kspec, norm));
        ans.r += diff.r + spec.r;
        ans.g += diff.g + spec.g;
        ans.b += diff.b + spec.b;
    }
    //cerr << amb.to_str() << '\n';
    //cerr << diff.to_str() << '\n';
    //cerr << spec.to_str() << '\n';
    //cerr << '\n';
    //ans.r = amb.r+diff.r+spec.r;
    //ans.g = amb.g+diff.g+spec.g;
    //ans.b = amb.b+diff.b+spec.b;
    return limit(ans);
}


db sqr(db a) {return a*a;}
db distsqred(vector<db> a, vector<db> b) {
    return sqr(a[0]-b[0])+sqr(a[1]-b[1])+sqr(a[2]-b[2]);
}
#define EPS 0.00001

template<typename T> struct matrix {
    int vert, horiz;
    vector<vector<T>> ma;
    matrix(int horizontal, int vertical) {
        horiz = horizontal, vert = vertical;
        ma.resize(horiz, vector<T>(vert, 0));
    }
    matrix(vector<vector<T>> m) {
        ma = m;
        horiz = (int)m.size();
        vert = (int)m[0].size();
    }
    matrix() {
        vert = 4;
        horiz = 0;
    }
    void clear() {
        horiz = 0;
        ma.clear();
    }
    void ident() {
        ma.clear();
        ma.resize(4, vector<T>(4, 0));
        vert = horiz = 4;
        for (int i = 0; i < 4; i++) ma[i][i] = 1;
    }
    vector<T>& operator[](int ind) {
        return ma[ind];
    }
    matrix<T> operator*(const matrix<T>& x) {
        assert(horiz == x.vert);
        matrix<T> ans(x.horiz, vert);
        for (int i = 0; i < x.horiz; i++) {
            for (int j = 0; j < vert; j++) {
                for (int k = 0; k < horiz; k++) {
                    ans[i][j] += ma[k][j]*x.ma[i][k];
                }
            }
        }
        return ans;
    }
    void print() {
        for (int j = 0; j < vert; j++) {
            cout << '[';
            for (int i = 0; i < horiz-1; i++) {
                cout << ma[i][j] << ' ';
            }
            cout << ma[horiz-1][j] << "]\n";
        }
        cout << '\n';
    }
    void add_point(T x, T y, T z = 0) {
        ma.push_back({x,y,z,1});
        horiz++;
    }
    void add_edge(T x0, T y0, T z0, T x1, T y1, T z1) {
        add_point(x0,y0,z0);
        add_point(x1,y1,z1);
    }
    void draw_lines(vector<vector<color>>& plot, vector<vector<db>>& zb, color p) {
        for (int i = 0; i < horiz; i+=2) {
            draw_line(plot, zb, ma[i][0], ma[i][1], ma[i][2], ma[i+1][0], ma[i+1][1], ma[i+1][2], p);
        }
    }

    void add_poly(vector<db> a, vector<db> b, vector<db> c) {
        assert(a.size() >= 3);
        assert(b.size() >= 3);
        assert(c.size() >= 3);
        add_point(a[0],a[1],a[2]);
        add_point(b[0],b[1],b[2]);
        add_point(c[0],c[1],c[2]);
    }

    void draw_poly(vector<vector<color>>& plot, vector<vector<db>>& zb, 
    color camb, vector<color> clight, vector<vector<db>> light, vector<db> view, vector<db> kamb, vector<db> kdiff, vector<db> kspec) {
        for (int i = 0; i < horiz; i += 3) {
            bool bad = 0;
            for (int a = 0; a < 3; a++) {
                for (int b = a+1; b < 3; b++) {
                    bad |= (ma[i+a] == ma[i+b]);
                    bad |= (distsqred(ma[i+a],ma[i+b]) < EPS);
                }
            }
            if (bad) continue;
            auto norm = normal_surface(ma[i], ma[i+1], ma[i+2]);
            db d = dot_prod(norm, {0,0,1});
            if (d > 0) {
                color abc(rand()%256, rand()%256, rand()%256);
                color c = get_lighting(camb, clight, light, view, norm, kamb, kdiff, kspec);
                //cerr << c.to_str() << '\n';
                scanline_convert(i, plot, zb, c);
            }
        }
    }

    void draw_scanline(db x0, db x1, db z0, db z1, db y, vector<vector<color>>& plot, vector<vector<db>>& zb, color c) {
        if (x0 > x1) {
            swap(x0,x1);
            swap(z0,z1);
        }
        db z = z0;
        x1 = ceil(x1);
        x0 = x0-1;
        for (int i = x0; i <= x1; i++) {
            if (0 <= i && i < (int)plot.size() && 0 <= y && y < (int)plot[0].size()) {
                if (z-zb[i][y] > 5 && !isnan(z)) {//this works better
                    plot[i][y] = c;
                    zb[i][y] = z;
                }
            }
            z = z0 + (i-(int)(x0)+1)*(z1-z0)/(x1-x0);
        }
    }

    void scanline_convert(int ind, vector<vector<color>>& plot, vector<vector<db>>& zb, color clr) {
        vector<db> t = ma[ind], m = ma[ind+1], b = ma[ind+2];
        vector<vector<db>> srt = {t,m,b};
        sort(srt.begin(), srt.end(), [](const vector<db>& p, const vector<db>& q){return p[1]<q[1];});
        t = srt[2], m = srt[1], b = srt[0];
        b[1]--;
        m[1]--;
        db x0,x1,z0,z1,y,ym,yt;
        x0 = x1 = b[0];
        z0 = z1 = b[2];
        y = b[1], yt = t[1], ym = m[1];
        bool switched = 0;
        while (y <= yt) {
            if (!switched && y >= ym) {
                switched = 1;
                x1 = m[0];
                z1 = m[2];
            }
            draw_scanline(x0,x1,z0,z1,y,plot,zb,clr);

            x0 = b[0] + (y-b[1]+1)*(t[0]-b[0])/(t[1]-b[1]);
            if (!switched) x1 = b[0] + (y-b[1]+1)*(m[0]-b[0])/(m[1]-b[1]);
            else x1 = m[0] + (y-m[1]+1)*(t[0]-m[0])/(t[1]-m[1]);

            z0 = b[2] + (y-b[1]+1)*(t[2]-b[2])/(t[1]-b[1]);
            if (!switched) z1 = b[2] + (y-b[1]+1)*(m[2]-b[2])/(m[1]-b[1]);
            else z1 = m[2] + (y-m[1]+1)*(t[2]-m[2])/(t[1]-m[1]);
            y++;
        }
    }
};

matrix<db> mk_translate(db x, db y, db z) {
    matrix<db> ans;
    ans.ident();
    ans[3][0] = x;
    ans[3][1] = y;
    ans[3][2] = z;
    return ans;
}

matrix<db> mk_scale(db x, db y, db z) {
    matrix<db> ans;
    ans.ident();
    ans[0][0] = x;
    ans[1][1] = y;
    ans[2][2] = z;
    return ans;
}

matrix<db> mk_rotX(db theta) {
    matrix<db> ans;
    ans.ident();
    theta *= M_PI/180;//convert to radians
    ans[1][1] = cos(theta);
    ans[1][2] = sin(theta);
    ans[2][1] = -sin(theta);
    ans[2][2] = cos(theta);
    return ans;
}

matrix<db> mk_rotY(db theta) {
    matrix<db> ans;
    ans.ident();
    theta *= M_PI/180;//convert to radians
    ans[0][0] = cos(theta);
    ans[2][0] = sin(theta);
    ans[0][2] = -sin(theta);
    ans[2][2] = cos(theta);
    return ans;
}

matrix<db> mk_rotZ(db theta) {
    matrix<db> ans;
    ans.ident();
    theta *= M_PI/180;//convert to radians
    ans[0][0] = cos(theta);
    ans[0][1] = sin(theta);
    ans[1][0] = -sin(theta);
    ans[1][1] = cos(theta);
    return ans;
}

#define INCR 0.0005

void circle(matrix<db>& m, db cx, db cy, db cz, db r) {
    for (db t = INCR; t <= 1+INCR; t += INCR) {
        m.add_edge(r*cos(2*M_PI*t)+cx, r*sin(2*M_PI*t)+cy, cz, r*cos(2*M_PI*(t-INCR))+cx, r*sin(2*M_PI*(t-INCR))+cy, cz);
    }
}

db polynomial(vector<db> pts, db x, db y, bool c=1) {
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

void bezier_curve(matrix<db>& m, db x0, db y0, db x1, db y1, db x2, db y2, db x3, db y3) {
    for (db t = INCR; t <= 1+INCR; t += INCR) {
        db xcurr = polynomial({x0,x1,x2,x3}, 1-t, t);
        db ycurr = polynomial({y0,y1,y2,y3}, 1-t, t);
        db xprev = polynomial({x0,x1,x2,x3}, 1-t+INCR, t-INCR);
        db yprev = polynomial({y0,y1,y2,y3}, 1-t+INCR, t-INCR);
        m.add_edge(xprev, yprev, 0, xcurr, ycurr, 0);
    }
}

void general_bezier(matrix<db>& m, vector<db> vx, vector<db> vy) {
    for (db t = INCR; t <= 1+INCR; t += INCR) {
        db xcurr = polynomial(vx, 1-t, t);
        db ycurr = polynomial(vy, 1-t, t);
        db xprev = polynomial(vx, 1-t+INCR, t-INCR);
        db yprev = polynomial(vy, 1-t+INCR, t-INCR);
        m.add_edge(xprev, yprev, 0, xcurr, ycurr, 0);
    }
}

void hermite_curve(matrix<db>& m, db x0, db y0, db x1, db y1, db rx0, db ry0, db rx1, db ry1) {
    matrix<db> h_inv({{2,-3,0,1},{-2,3,0,0},{1,-2,1,0},{1,-1,0,0}}),
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

void add_box(matrix<db>& m, db x, db y, db z, db w, db h, db d) {
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

void add_sphere(matrix<db>& m, db cx, db cy, db cz, db r, const int step) {
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

void add_torus(matrix<db>& m, db cx, db cy, db cz, db r1, db r2, int step) {
    auto pts = gen_torus(cx, cy, cz, r1, r2, step);
    for (int rot = 0; rot < (int)pts.size()-1; rot++) {
        for (int circ = 0; circ < (int)pts[0].size()-1; circ++) {
            m.add_poly(pts[rot][circ], pts[rot+1][circ], pts[rot+1][circ+1]);
            m.add_poly(pts[rot+1][circ+1], pts[rot][circ+1], pts[rot][circ]);
        }
    }
}

vector<vector<vector<db>>> gen_bezier(vector<db> vx, vector<db> vy, int step) {
    vector<vector<vector<db>>> ans(step+1);
    matrix<db> m;
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

void add_bezier_shape(matrix<db>& m, vector<db> vx, vector<db> vy, int step) {
    auto pts = gen_bezier(vx, vy, step);
    for (int rot = 0; rot < (int)pts.size()-1; rot++) {
        for (int circ = 0; circ < (int)pts[0].size()-1; circ++) {
            //m.add_poly(pts[rot][circ], pts[rot+1][circ], pts[rot+1][circ+1]);
            //m.add_poly(pts[rot+1][circ+1], pts[rot][circ+1], pts[rot][circ]);
            m.add_poly(pts[rot+1][circ+1], pts[rot+1][circ], pts[rot][circ]);
            m.add_poly(pts[rot][circ], pts[rot][circ+1], pts[rot+1][circ+1]);
        }
    }
}