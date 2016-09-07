# Ben Coscia
# Homework 10
# Part II Question 3

import math
import scipy.optimize as sci
import matplotlib.pyplot as plt
import numpy as np

# Given Data

I = [50, 80, 130, 200, 250, 350, 450, 550, 700]
P = [99, 177, 202, 248, 229, 219, 173, 142, 72]

# Fit the given data to the mathematical model:
# P = Pm*(I/Isat)*exp(-(I/Isat) + 1)   where
#   P == Photosynthesis Rate (mg / m^3 d)
#   Pm == Maximum Photosynthesis Rate (mg / m^3 d)
#   I == Solar Radiation (mu E /m^4 s)
#   Isat == Optimal Solar Radiation (mu E /m^4 s)

# The simplest way to do it is to minimize the sum of the residuals so that is what I am going to do
# The residual at each point is calculated as (y - f(a0, a1, x))**2
# for this specific case it will be (P - Pm*(I/Isat)*exp(-(I/Isat) + 1))**2

PmIsat = [100, 10]  # initial guesses at Pm and Isat. Format: [Pm, Isat]

def chi(PmIsat):
    sum_chi = 0
    for i in range(0, len(I)):
        sum_chi += (P[i] - PmIsat[0]*(I[i]/PmIsat[1])*math.exp(-(I[i]/PmIsat[1]) + 1))**2
    return sum_chi

A = sci.minimize(chi, PmIsat)


def goodness_of_fit(P, I, A):
    St = 0
    Sr = 0
    mean = np.mean(P)
    for i in range(0, len(P)):
            Sr += (P[i] - A.x[0]*(I[i]/A.x[1])*math.exp(-(I[i]/A.x[1]) + 1))**2
            St += (P[i] - mean)**2
    s = math.sqrt(Sr/(len(P) - len(A.x)))
    R_squared = 1 - (Sr/St)
    return float(R_squared), s

def plot_fit(A, I):
    y = np.zeros((len(I)))
    for i in range(0, len(I)):
        y[i] = A.x[0]*(I[i]/A.x[1])*math.exp(-(I[i]/A.x[1]) + 1)
    return y

R_squared, s = goodness_of_fit(P, I, A)
P_fit = plot_fit(A, I)
plt.plot(I, P)
plt.plot(I, P_fit)
plt.title('P versus I, R-squared = %s, s = %s' %(R_squared, s))
plt.xlabel('Solar Radiation (mu E /m^4 s)')
plt.ylabel('Photosynthesis Rate (mg / m^3 d)')
plt.show()
