#pragma once

#include "types.hpp"
#include "matrix.hpp"

matrix mk_translate(db x, db y, db z);

matrix mk_scale(db x, db y, db z);

matrix mk_rotX(db theta);

matrix mk_rotY(db theta);

matrix mk_rotZ(db theta);