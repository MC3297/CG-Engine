#pragma once

#include "types.h"
#include "matrix.h"

matrix mk_translate(db x, db y, db z) {
    matrix ans;
    ans.ident();
    ans[3][0] = x;
    ans[3][1] = y;
    ans[3][2] = z;
    return ans;
}

matrix mk_scale(db x, db y, db z) {
    matrix ans;
    ans.ident();
    ans[0][0] = x;
    ans[1][1] = y;
    ans[2][2] = z;
    return ans;
}

matrix mk_rotX(db theta) {
    matrix ans;
    ans.ident();
    theta *= M_PI/180;//convert to radians
    ans[1][1] = cos(theta);
    ans[1][2] = sin(theta);
    ans[2][1] = -sin(theta);
    ans[2][2] = cos(theta);
    return ans;
}

matrix mk_rotY(db theta) {
    matrix ans;
    ans.ident();
    theta *= M_PI/180;//convert to radians
    ans[0][0] = cos(theta);
    ans[2][0] = sin(theta);
    ans[0][2] = -sin(theta);
    ans[2][2] = cos(theta);
    return ans;
}

matrix mk_rotZ(db theta) {
    matrix ans;
    ans.ident();
    theta *= M_PI/180;//convert to radians
    ans[0][0] = cos(theta);
    ans[0][1] = sin(theta);
    ans[1][0] = -sin(theta);
    ans[1][1] = cos(theta);
    return ans;
}