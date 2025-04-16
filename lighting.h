#pragma once

#include "types.h"
#include "color.h"
#include "vector_math.h"

color get_amb(color camb, Vec3 kamb) {
    return color(camb.r*kamb[0], camb.g*kamb[1], camb.b*kamb[2]);
}

color get_diff(color clight, Vec3 light, Vec3 kdiff, Vec3 norm) {
    db cosin = dot_prod(norm,light);
    return color(clight.r*kdiff[0]*cosin, clight.g*kdiff[1]*cosin, clight.b*kdiff[2]*cosin);
}

color get_spec(color clight, Vec3 light, Vec3 view, Vec3 kspec, Vec3 norm) {
    db cf = 2*dot_prod(norm,light);
    Vec3 cfn = {norm[0]*cf, norm[1]*cf, norm[2]*cf};
    Vec3 subtr = {cfn[0]-light[0], cfn[1]-light[1], cfn[2]-light[2]};
    db cosin = pow(dot_prod(subtr,view),7);
    return color(clight.r*kspec[0]*cosin, clight.g*kspec[1]*cosin, clight.b*kspec[2]*cosin);
}

color limit(color c) {
    return color(min(max(0,c.r),255), min(max(0,c.g),255), min(max(0,c.b),255));
}

color get_lighting(color camb, vector<color> clight, vector<Vec3> light, Vec3 view, Vec3 norm, Vec3 kamb, Vec3 kdiff, Vec3 kspec) {
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
    return limit(ans);
}