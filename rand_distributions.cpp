#include <math.h>
#include <stdint.h>
#include "rand_distributions.hpp"
#include "pcg/pcg_random.hpp"

double nextUc01o(pcg64 &r) { // [0,1)
    //Originally from http://xoshiro.di.unimi.it
    //rationals of the form u/2^53 are produced (u is an integer);
    //if we include 0 but skip 1, there are exactly 53 bits worth
    //of floats of this form.
    return (r() >> 11) * (1. / (uint64_t(1) << 53));
}

double nextUo01c(pcg64 &r) { // (0,1]
    //Almost the same as nextUc01o, but we want to skip 0
    //and include 1. There are still 53 bits worth of floats
    //of the form u/2^53
    return ((r() >> 11) + 1) * (1. / (uint64_t(1) << 53));
}

double nextUM1P1_cl(pcg64 &r) { // [-1,1]
    //This is a little bit trickier. We can easily get [-1,1) or (-1,1]
    //Since for the open intervals there are 54 bits worth of v/2^53
    //if we allow v to be an unsigned integer. Instead we just generate
    //ints with a range of 1 more than 54 bits worth, which gives us
    //the closed ranges. This comes at a cost of ~.1% rejection, a
    //modulus operation, and a comparison operation. Not too much!
    uint64_t u = r();
    //Ensure fairness of sampling integers by rejecting.
    //We'll take the modulus later, so just reject if we have the leftover
    //portion. As in, if we wanted to generate numbers 0-6 out of RVs
    //between 0 and 15, we can take the RV % 7. But then we have a few left
    //over (14 & 15 map to an extra 0 or 1), so we'd reject numbers > 13.
    //Hope I did all those constants right!!!
    if(u > (1023*((uint64_t(1) << 54) + 1) - 1)) {
        return nextUM1P1_cl(r);
    }
    //Take modulus of number
    int64_t mod = u % ((uint64_t(1) << 54) + 1);
    return (mod - (int64_t(1) << 53))*(1. / (uint64_t(1) << 53));
}

double nextExp(pcg64 &r) {
    //Simple inverse transform. We want to use (0,1] since 0 would produce
    //infinity.
    return -log(nextUo01c(r));
}

void next2Norm(pcg64 &r, double *z1, double *z2) {
    //Marsaglia's polar box-muller transform.
    //Need [-1,1]
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