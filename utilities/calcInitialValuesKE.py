#!/usr/bin/env python3
#-*- coding:utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

C_mu = 0.09
nu = 1E-5
diameter  = 0.04
Re = 400000

U = Re*nu/diameter

lenght_scale = 0.16*Re**(-1/8)
k = (3/2)*(lenght_scale*U)**2
epsilon = C_mu**(3/4)*k**(3/2)/lenght_scale

print(f"U = {U}")
print(f"lenght_scale = {lenght_scale}")
print(f"k = {k}")
print(f"epsilon = {epsilon}")

## Mesh parameters should be given as inputs
cell_exp_ratio = 2
n = 70
Sum = 4E-5
## Calculate yPlus based on meshGeometry
# Note: the approach used here is based on geometric progressions
#       where q: ratio
#             a1: first term
#             an: last term
#
# In case cell expansion ratio (an/a1) is different from 1
# an/a1 = cell_exp_ratio -> q = cell_exp_ratio**(1/(n-1))
q = cell_exp_ratio**(1/(n-1))

# Geometric progression sum (Sum) = a0*(q**n-1)/(q-1)
# for this case, a0 is the least value term of the G.P.
# which also means the least distance from the wall for the mesh
# based on the total height of the patch (numeric equal to Sum)
a0 = Sum*(q-1)/(q**n-1)

yPlus = 3
y = nu*yPlus/U

print(f"Least wall distance = {a0}")
print(f"A value of yPlus = {yPlus} corresponds to y = {y}")