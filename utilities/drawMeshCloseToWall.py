#!/usr/bin/env python3
#-*- coding:utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

nu = 1E-5
diameter  = 0.04
Re = 400000

U = Re*nu/diameter
## Mesh parameters should be given as inputs
cell_exp_ratio = 4
n = 70
Sum = 1E-3
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

yPlus = 1
y_eq_yPlus = nu*yPlus/U

print(f"Least wall distance = {a1}")
print(f"A value of yPlus = {yPlus} corresponds to y = {y_eq_yPlus}")

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
