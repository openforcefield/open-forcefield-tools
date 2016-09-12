import pymbar
from pymbar import MBAR
from pymbar import testsystems
import numpy as np
import sys

(x_n, u_kn, N_k, s_n) = testsystems.HarmonicOscillatorsTestCase().sample(mode='u_kn')


print np.average(x_n[0:9])
print np.average(u_kn[0][0:9])

mbar = MBAR(u_kn, N_k)
A_n = x_n

print "Input for expectation of position: %s" % (A_n)

(A_ij, dA_ij) = mbar.computeExpectations(A_n)

A_n2 = u_kn[0,:]
print "Input for expectation of reduced potential: %s" % (A_n2)

(A_ij2, dA_ij2) = mbar.computeExpectations(A_n2)

print A_ij
print A_ij2
