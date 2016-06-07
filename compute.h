#ifndef COMPUTE_H
#define COMPUTE_H

#include "common.h"
#include "sparse.h"
#include "dense.h"

extern void multiply(Sparse* a, Dense* b, Dense* c);

#endif
