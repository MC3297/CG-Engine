#include "matrix.hpp"

db sqr(db a) {return a*a;}
db distsqred(Point a, Point b) {
    return sqr(a[0]-b[0])+sqr(a[1]-b[1])+sqr(a[2]-b[2]);
}

matrix::matrix(int horizontal, int vertical) : horiz(horizontal), vert(vertical) {
    ma.resize(horiz);
}

matrix::matrix(vector<Point> m) : ma(m), horiz(m.size()), vert(m[0].size()) {}

matrix::matrix() : vert(4), horiz(0) {}

void matrix::clear() {
    horiz = 0;
    ma.clear();
}

void matrix::ident() {
    ma.clear();
    ma.resize(4);
    vert = horiz = 4;
    for (int i = 0; i < 4; i++) ma[i][i] = 1;
}

Point& matrix::operator[](int ind) {
    return ma[ind];
}

matrix matrix::operator*(const matrix& x) {
    assert(horiz == x.vert);
    matrix ans(x.horiz, vert);
    for (int i = 0; i < x.horiz; i++) {
        for (int j = 0; j < vert; j++) {
            for (int k = 0; k < horiz; k++) {
                ans[i][j] += ma[k][j]*x.ma[i][k];
            }
        }
    }
    return ans;
}

void matrix::print() {
    for (int j = 0; j < vert; j++) {
        cout << '[';
        for (int i = 0; i < horiz-1; i++) {
            cout << ma[i][j] << ' ';
        }
        cout << ma[horiz-1][j] << "]\n";
    }
    cout << '\n';
}

void matrix::add_point(db x, db y, db z = 0) {
    ma.push_back({x,y,z,1});
    horiz++;
}

void matrix::add_edge(db x0, db y0, db z0, db x1, db y1, db z1) {
    add_point(x0,y0,z0);
    add_point(x1,y1,z1);
}

void matrix::draw_lines(vector<vector<color>>& plot, vector<vector<db>>& zb, color p) {
    for (int i = 0; i < horiz; i+=2) {
        draw_line(plot, zb, ma[i][0], ma[i][1], ma[i][2], ma[i+1][0], ma[i+1][1], ma[i+1][2], p);
    }
}

void matrix::add_poly(vector<db> a, vector<db> b, vector<db> c) {
    assert(a.size() >= 3);
    assert(b.size() >= 3);
    assert(c.size() >= 3);
    add_point(a[0],a[1],a[2]);
    add_point(b[0],b[1],b[2]);
    add_point(c[0],c[1],c[2]);
}

#define EPS 0.00001
void matrix::draw_poly(vector<vector<color>>& plot, vector<vector<db>>& zb, color camb, vector<color> clight, vector<Vec3> light, Vec3 view, Vec3 kamb, Vec3 kdiff, Vec3 kspec) {
    for (int i = 0; i < horiz; i += 3) {
        Vec3 norm = normal_surface(ma[i], ma[i+1], ma[i+2]);
        db d = dot_prod(norm, {0,0,1});
        if (d > 0) {
            color c = get_lighting(camb, clight, light, view, norm, kamb, kdiff, kspec);
            scanline_convert(i, plot, zb, c);
        }
    }
}

void matrix::draw_scanline(db x0, db x1, db z0, db z1, db y, vector<vector<color>>& plot, vector<vector<db>>& zb, color c) {
    if (x0 > x1) {
        swap(x0,x1);
        swap(z0,z1);
    }
    db z = z0;
    x1 = ceil(x1);
    x0 = x0-1;
    for (int i = x0; i <= x1; i++) {
        if (0 <= i && i < (int)plot.size() && 0 <= y && y < (int)plot[0].size()) {
            if (z-zb[i][y] >= 7 && !isnan(z)) {//this works better
                plot[i][y] = c;
                zb[i][y] = z;
            }
        }
        z = z0 + (i-(int)(x0)+1)*(z1-z0)/(x1-x0);
    }
}

void matrix::scanline_convert(int ind, vector<vector<color>>& plot, vector<vector<db>>& zb, color clr) {
    Point t = ma[ind], m = ma[ind+1], b = ma[ind+2];
    vector<Point> srt = {t,m,b};
    sort(srt.begin(), srt.end(), [](const Point& p, const Point& q){return p[1]<q[1];});
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