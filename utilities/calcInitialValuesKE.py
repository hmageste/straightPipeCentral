#!/usr/bin/env python3
#-*- coding:utf-8 -*-

import argparse

# construct the argument parser and parse the arguments
parser = argparse.ArgumentParser(description="Calculate initial values for k epsilon model", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-re","--reynolds", required=False, type=float, help="Flow Reynolds number")
group.add_argument("-u","--velocity", required=False, type=float, help="Flow velocity")
parser.add_argument("-nu","--kinematic_viscosity", required=False, type=float, default=1E-5, help="Fluid kinematic_viscosity. nu = mu/rho")
parser.add_argument("-l","--characteristic_length", required=False, type=float, default=0.04,
    help="Characteristic length of the flow. For pipes: hidraulic diameter; flat-plate: plate length")
parser.add_argument("-v", "--verbosity", action="count", default=0,
    help="increase output verbosity")

args = parser.parse_args()

if args.verbosity >= 2:
    print("Running '{}'".format(__file__))
if args.verbosity >= 1:
    print("The Reynolds number is given by Re = U*L/nu")

C_mu = 0.09
nu = args.kinematic_viscosity
diameter  = args.characteristic_length
U = args.velocity
Re = args.reynolds

if Re:
    U = Re*nu/diameter
elif U:
    Re = U*diameter/nu
else:
    print("Provide either -re or -u")
    

lenght_scale = 0.16*Re**(-1/8)
k = (3/2)*(lenght_scale*U)**2
epsilon = C_mu**(3/4)*k**(3/2)/lenght_scale

print(f"U = {U}")
print(f"Re = {Re}")
print(f"lenght_scale = {lenght_scale}")
print(f"k = {k}")
print(f"epsilon = {epsilon}")