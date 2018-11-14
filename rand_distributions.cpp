#include <math.h>
#include <stdint.h>
#include "rand_distributions.hpp"
#include "pcg/pcg_random.hpp"

double nextUc01o(pcg64 &r) { // [0,1)
    return (r() >> 11) * (1. / (uint64_t(1) << 53));
}

double nextUo01c(pcg64 &r) { // (0,1]
    return ((r() >> 11) + 1) * (1. / (uint64_t(1) << 53));
}

double nextUM1P1_cl(pcg64 &r) { // [-1,1]
    uint64_t u = r();
    if(u > (1023*((uint64_t(1) << 54) + 1) - 1)) {
        return nextUM1P1_cl(r);
    }
    int64_t mod = u % ((uint64_t(1) << 54) + 1);
    return (mod - (int64_t(1) << 53))*(1. / (uint64_t(1) << 53));
}

double nextExp(pcg64 &r) {
    return -log(nextUo01c(r));
}

void next2Norm(pcg64 &r, double *z1, double *z2) {
    double u = nextUM1P1_cl(r);
    double v = nextUM1P1_cl(r);
    double s = u*u + v*v;
    while((s == 0.0) || (s >= 1.0)) {
        u = nextUM1P1_cl(r);
        v = nextUM1P1_cl(r);
        s = u*u + v*v;
    }
    double sqrtlog = sqrt(-2*log(s)/s);
    *z1 = u*sqrtlog;
    *z2 = v*sqrtlog;
}