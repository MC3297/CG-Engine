#pragma once

#include "types.h"
#include "color.h"

//colors
color camb(50,50,50);
vector<color> clight = {color(0,255,255)};
vector<vector<db>> light = {{0.5,0.75,1}};
vector<db> view = {0,0,1};
vector<db> kamb = {0.2,0.2,0.2};
vector<db> kdiff = {0.6,0.6,0.6};
vector<db> kspec = {0.82,0.82,0.82};