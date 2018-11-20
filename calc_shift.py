#!/usr/bin/env python3

import numpy as np

def ltShift(alpha, beta, tau, totalN, tHold, kappa, tDip):
    return (tau -
            (tHold)/
             (tHold/tau + np.log(
                (3*alpha*(tDip+2*kappa)**2 + 2*beta*totalN*(3*kappa + 2*tDip))
                /
                (3*alpha*(tDip+2*kappa)**2 + 4*beta*totalN*(3*kappa + 2*tDip))
              )
             )
           )

def alphabeta(points):
    sumInvErr = sum([1/(e**2) for (x,y,e) in points])
    sumXInvErr = sum([x/(e**2) for (x,y,e) in points])
    sumXXInvErr = sum([x**2/(e**2) for (x,y,e) in points])
    sumYInvErr = sum([y/(e**2) for (x,y,e) in points])
    sumXYInvErr = sum([x*y/(e**2) for (x,y,e) in points])
    #points is [(x,y,yerr), ...]
    delta = np.linalg.det(np.matrix([[sumInvErr, sumXInvErr],[sumXInvErr, sumXXInvErr]]))
    alpha = (1/delta)*np.linalg.det(np.matrix([[sumYInvErr, sumXInvErr],[sumXYInvErr, sumXXInvErr]]))
    beta = (1/delta)*np.linalg.det(np.matrix([[sumInvErr, sumYInvErr],[sumXInvErr, sumXYInvErr]]))
    alphaErr = np.sqrt((1/delta)*sumXXInvErr)
    betaErr = np.sqrt((1/delta)*sumInvErr)

    return (alpha, alphaErr), (beta, betaErr)

#Example efficiencies, taken from running sim from commit be2de83aaa68b76bc47420fa23c4b49661997d86
effs = [(100.000000, 0.906779, 9.016094e-05),
        (100.000000, 0.906759, 9.030045e-05),
        (5000.000000, 0.889910, 8.864058e-05),
        (5000.000000, 0.889967, 8.864998e-05),
        (10000.000000, 0.872398, 8.723504e-05),
        (10000.000000, 0.872419, 8.701708e-05),]

#Example efficiencies, taken from running sim from commit 77a5cc129181ea397132bf119848ef85a0e32633
#Use the outputs corresponding to the DT corrected output
effs = [(100.000000, 0.907139, 9.575018e-05),
        (100.000000, 0.907119, 9.559649e-05),
        (5000.000000, 0.907895, 1.544868e-03),
        (5000.000000, 0.907949, 1.442664e-03),
        (10000.000000, 0.908696, 2.827326e-03),
        (10000.000000, 0.908724, 2.656485e-03),]

ab = alphabeta(effs)
alpha = ab[0][0]
beta = ab[1][0]
print("((Alpha, err), (Beta, err)):",ab)
print("\n\nLT Shift:", ltShift(ab[0][0], ab[1][0], 877.7, 20000, 1550, 4.77, 10))
print("Rough 1-sigma bounds:", ltShift(ab[0][0]+ab[0][1], ab[1][0]+ab[1][1], 877.7, 20000, 1550, 4.77, 10),
      ltShift(ab[0][0]+ab[0][1], ab[1][0]-ab[1][1], 877.7, 20000, 1550, 4.77, 10),
      ltShift(ab[0][0]-ab[0][1], ab[1][0]+ab[1][1], 877.7, 20000, 1550, 4.77, 10),
      ltShift(ab[0][0]-ab[0][1], ab[1][0]-ab[1][1], 877.7, 20000, 1550, 4.77, 10))