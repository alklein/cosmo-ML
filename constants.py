#!/usr/bin/env python                                                                 

import os
import sys
import math
import numpy as np

from random import *
from numpy import log10
from optparse import OptionParser, OptionGroup

ID, D_GC, PVMAX, VMAX, RMAX, MTIDAL, RTIDAL, X, Y, Z, Vx, Vy, Vz, M300pc, M600pc = range(15)
Rad_ext, J_nfw, J_approx = range(8, 11)
mm, rr, tt, pp, dd, jj = range(6)
x_earth = 8.5 
J_D = 2.e12

R_vir_VL = 402.
R_vir_AQ = 433.
km_per_kpc = 3.086e16

ID, X, Y, Z, VX, VY, VZ = range(7)
N200h, M200h, R200h, N500h, M500h, R500h = range(7, 13)
N200c, M200c, R200c, N500c, M500c, R500c = range(13, 19)
Np = 128**3
Nh = 62897

scalefactor = 1. # for z = 0
a = scalefactor
box = 200.
Ndm = 512.
h0 = 0.7
om = 0.27
Mpc2cm = 3.08560e+24
H0_cgs = 3.24086e-18

x_unit = box/Ndm
time_unit = 2/(3*H0_cgs*h0*(om**.5))*a**2 #[s?]#
vel_unit  = (x_unit/time_unit)/(1e5) #[km/s]#

if __name__ == '__main__':
    print 'Indices:'
    print ID, X, Y, Z, VX, VY, VZ, \
          N200h, M200h, R200h, N500h, M500h, R500h, \
          N200c, M200c, R200c, N500c, M500c, R500c
    print 'Number of particles:',Np
    print 'Number of halos:',Nh
