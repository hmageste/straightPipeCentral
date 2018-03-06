#!/usr/bin/env python3
#-*- coding:utf-8 -*-

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