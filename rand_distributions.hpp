#ifndef RAND_DISTRIBUTIONS_H
#define RAND_DISTRIBUTIONS_H

#include "pcg/pcg_random.hpp"

double nextUc01o(pcg64 &r); // [0,1)

double nextUo01c(pcg64 &r); // (0,1]

double nextUM1P1_cl(pcg64 &r); // [-1,1]

double nextExp(pcg64 &r);

void next2Norm(pcg64 &r, double *z1, double *z2);

#endif /*RAND_DISTRIBUTIONS_H*/