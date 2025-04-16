#include "template.cpp"
#include "types.h"
#include "color.h"
#include "lighting_constants.h"
#include "vector_math.h"
#include "bresenham_line.h"
#include "matrix.h"
#include "matrix_transform.h"
#include "lighting.h"

template<class T> void read(T &x) {cin >> x;}
template<class H, class... T> void read(H &h, T &...t) { read(h); read(t...); }

const int N = 500;
const int M = 500;
const int resolution = 30;

void parse(vector<vector<color>>& plot, vector<vector<db>>& zb, stack<matrix>& stacc, string s, color clr = color(0,255,0)) {
    if (s[0] == '#') {
        string qwerty;
        getline(cin, qwerty);
    }
    else if (s == "clear") {
        plot.clear();
        plot.resize(N, vector<color>(M, color(255,255,255)));
        zb.clear();
        zb.resize(N, vector<db>(M, -69420/0.0));
    }
    else if (s == "push") {
        stacc.push(stacc.top());
    }
    else if (s == "pop") {
        stacc.pop();
    }
    else if (s == "scale") {
        db sx, sy, sz;
        read(sx,sy,sz);
        stacc.top() = stacc.top()*mk_scale(sx,sy,sz);
    }
    else if (s == "move") {
        db tx, ty, tz;
        read(tx,ty,tz);
        stacc.top() = stacc.top()*mk_translate(tx,ty,tz);
    }
    else if (s == "rotate") {
        char axis; db theta;
        read(axis,theta);
        if (axis == 'x') stacc.top() = stacc.top()*mk_rotX(theta);
        else if (axis == 'y') stacc.top() = stacc.top()*mk_rotY(theta);
        else if (axis == 'z') stacc.top() = stacc.top()*mk_rotZ(theta);
    }
    else if (s == "line") {
        db x0, y0, z0, x1, y1, z1;
        read(x0, y0, z0, x1, y1, z1);
        matrix temp;
        temp.add_edge(x0, y0, z0, x1, y1, z1);
        (stacc.top() * temp).draw_lines(plot, zb, clr);
    }
    else if (s == "circle") {
        db cx, cy, cz, r;
        read(cx, cy, cz, r);
        matrix temp;
        circle(temp, cx, cy, cz, r);
        (stacc.top() * temp).draw_lines(plot, zb, clr);
    }
    else if (s == "hermite") {
        db x0, y0, x1, y1, rx0, ry0, rx1, ry1;
        read(x0, y0, x1, y1, rx0, ry0, rx1, ry1);
        matrix temp;
        hermite_curve(temp, x0, y0, x1, y1, rx0, ry0, rx1, ry1);
        (stacc.top() * temp).draw_lines(plot, zb, clr);
    }
    else if (s == "bezier") {
        db x0, y0, x1, y1, x2, y2, x3, y3;
        read(x0, y0, x1, y1, x2, y2, x3, y3);
        matrix temp;
        bezier_curve(temp, x0, y0, x1, y1, x2, y2, x3, y3);
        (stacc.top() * temp).draw_lines(plot, zb, clr);
    }
    /*else if (s == "beziershape") {
        int n; cin >> n;
        vector<db> vx(n), vy(n);
        for (db& i: vx) cin >> i;
        for (db& i: vy) cin >> i;
        matrix temp;
        add_bezier_shape(temp, vx, vy, resolution);
        (stacc.top() * temp).draw_poly(plot, zb, camb, clight, light, view, kamb, kdiff, kspec);
    }*/
    else if (s == "triangle") {
        db x0, y0, z0, x1, y1, z1, x2, y2, z2;
        read(x0, y0, z0, x1, y1, z1, x2, y2, z2);
        matrix temp;
        temp.add_poly({x0,y0,z0}, {x1,y1,z1}, {x2,y2,z2});
        (stacc.top() * temp).draw_poly(plot, zb, camb, clight, light, view, kamb, kdiff, kspec);
    }
    else if (s == "box") {
        db x, y, z, w, h, d;
        read(x, y, z, w, h, d);
        matrix temp;
        add_box(temp, x, y, z, w, h, d);
        (stacc.top() * temp).draw_poly(plot, zb, camb, clight, light, view, kamb, kdiff, kspec);
    }
    else if (s == "sphere") {
        db x, y, z, r;
        read(x, y, z, r);
        matrix temp;
        add_sphere(temp, x, y, z, r, resolution);
        (stacc.top() * temp).draw_poly(plot, zb, camb, clight, light, view, kamb, kdiff, kspec);
    }
    else if (s == "torus") {
        db x, y, z, r1, r2;
        read(x, y, z, r1, r2);
        matrix temp;
        add_torus(temp, x, y, z, r1, r2, resolution);
        (stacc.top() * temp).draw_poly(plot, zb, camb, clight, light, view, kamb, kdiff, kspec);
    }
    else if (s == "display") {
        FILE *f = popen("display", "w");
        fprintf(f, "P3\n%d %d\n%d\n", N, M, 255);
        for (int j = (int)plot[0].size()-1; j >= 0; j--)  {
            for (int i = 0; i < (int)plot.size(); i++) {
                fprintf(f, "%d %d %d ", plot[i][j].r, plot[i][j].g, plot[i][j].b);
            }
            fprintf(f, "\n");
        }
        pclose(f);
    }
    else if (s == "save") {
        string filename;
        cin >> filename;
        freopen("temp.ppm","w",stdout);
        cout << "P3" << N << ' ' << M << '\n' << 255 << '\n';
        for (int j = (int)plot[0].size()-1; j >= 0; j--)  {
            for (int i = 0; i < (int)plot.size(); i++) {
                cout << plot[i][j].to_str();
            }
            cout << '\n';
        }
        fclose(stdout);
        system(("pnmtopng temp.ppm > "+filename).c_str());//please work
    }
    else {
        cout << "unknown command " << s << '\n';
    }
}

vector<vector<color>> grid(N, vector<color>(M, color(255,255,255)));
vector<vector<db>> zb(N, vector<db>(M, -69420/0.0));
color c(255,0,0);

int main() {
    stack<matrix> stacc;
    stacc.push(matrix());
    stacc.top().ident();
    freopen("script","r",stdin);
    while (cin.peek() != -1) {
        string s; cin >> s;
        if (s == "quit") break;
        parse(grid, zb, stacc, s, c);
    }
    return 0;
}