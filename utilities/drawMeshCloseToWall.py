#!/usr/bin/env python3
#-*- coding:utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from math import log10, exp, sqrt

try:
    from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
except ImportError:
    raise ImportError("Error loading PyFoam, download it with the command 'pip install PyFoam' if not installed!")

cs_dict = ParsedParameterFile("constant/caseSettings")

def calcDarcyFrictionFactor( Re, Dh = 1.0, roughness = 0.0, th = 1e-6 ):
    f = 0
    if Re > 0 and Re <= 2320: # laminar flow
        f = 64.0/Re
    elif Re > 2320 and Re < 4000:
        print ("No valid correlation exists for this Reynolds number")
    elif Re >= 4000: # turbulent flow
        # First guess based on Zigrang, Sylvester (1982) explict approximation
        f = (-2.0*log10(roughness/3.7 - 5.02/Re*log10(roughness - 5.02/Re*log10(roughness/3.7 + 13.0/Re))))**-2
        while True:
            fPrevious = f
            f = (-2.0 * log10(roughness/Dh/3.7 + 2.51/Re/sqrt(f)))**-2
            if (abs(f - fPrevious) <= th): break
    else:
        print ("Re must be a positive number different from 0")

    return f

# if __name__ == '__main__':
nu = 1E-5
diameter  = 0.01
Re = 400000

f = calcDarcyFrictionFactor(Re)

U = Re*nu/diameter
## Mesh parameters should be given as inputs
cell_exp_ratio = 500
n = 70
Sum = diameter
## Calculate yPlus based on meshGeometry
# Note: the approach used here is based on geometric progressions
#       where q: ratio
#             a1: first term
#             an: last term
#
# In case cell expansion ratio (an/a1) is different from 1
# an/a1 = cell_exp_ratio -> q = cell_exp_ratio**(1/(n-1))
q = cell_exp_ratio**(1/(n-1))

# Geometric progression sum (Sum) = a1*(q**n-1)/(q-1)
# for this case, a1 is the least value term of the G.P.
# which also means the least distance from the wall for the mesh
# based on the total height of the patch (numeric equal to Sum)
a1 = Sum*(q-1)/(q**n-1)

u_tau = sqrt(f/8.0)*U
yPlus = 1
y_eq_yPlus = nu*yPlus/u_tau

print(f"Least wall distance = {a1}")
print(f"A value of yPlus = {yPlus} corresponds to y = {y_eq_yPlus} m")

#pg = np.arange(0, n-1)
pg = np.linspace(0,n-1,n-1)
for idx in range(len(pg)):
    pg[idx] = a1*q**(idx)

xValues = np.ones(len(pg))

plt.plot(xValues, pg, ls='None', marker='o', ms='.5')
plt.hlines(y_eq_yPlus, 0.0, 2.0)
plt.hlines(2*y_eq_yPlus, 0.0, 2.0)
plt.hlines(3*y_eq_yPlus, 0.0, 2.0)
plt.hlines(5*y_eq_yPlus, 0.0, 2.0)
plt.hlines(30*y_eq_yPlus, 0.0, 2.0)
plt.show()
