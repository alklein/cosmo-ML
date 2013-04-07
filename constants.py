#!/usr/bin/env python                                                                 

import os
import sys
import math
import numpy as np

from random import *
from numpy import log10
from optparse import OptionParser, OptionGroup

# constants for Via Lactea II and Aquarius A:
ID, D_GC, PVMAX, VMAX, RMAX, MTIDAL, RTIDAL, X, Y, Z, Vx, Vy, Vz, M300pc, M600pc = range(15)
Rad_ext, J_nfw, J_approx = range(8, 11)
mm, rr, tt, pp, dd, jj = range(6)
x_earth = 8.5 
J_D = 2.e12
R_vir_VL = 402.
R_vir_AQ = 433.
km_per_kpc = 3.086e16

# Hy Trac cosmological parameters
omega_matter = 0.27
omega_lambda = 0.73
omega_baryon = 0.045
omega_rad    = 8.53E-5
hubble0      = 0.70
ns_init      = 0.96
sigma_8      = 0.80
w_de         = -1.0

# Hy Trac sim indices
ID, X, Y, Z, VX, VY, VZ = range(7) # for halos and particles
H_ID = 7 # for particles only
N200a, M200a, R200a, s200a, N500a, M500a, R500a, s500a = range(7, 15) # for halos only
N200c, M200c, R200c, s200c, N500c, M500c, R500c, s500c = range(15, 23) # for halos only

# Hy Trac particle + halo counts
Np = 128**3 # TODO: verify
Nh = 62897 # TODO: verify

# Hy Trac sim parameters
scalefactor = 1. # for z = 0
a = scalefactor
box = 200.
Ndm = 512.
h0 = 0.7
om = 0.27
Mpc2cm = 3.08560e+24
H0_cgs = 3.24086e-18

# Hy Trac conversion constants
x_unit = box/Ndm
time_unit = 2/(3*H0_cgs*h0*(om**.5))*a**2 #[s?]#
vel_unit  = (x_unit/time_unit)/(1e5) #[km/s]#

if __name__ == '__main__':
    print 'Indices:'
    print ID, X, Y, Z, VX, VY, VZ, \
          N200a, M200a, R200a, s200a, N500a, M500a, R500a, s500a, \
          N200c, M200c, R200c, s200c, N500c, M500c, R500c, s500c
    print 'Number of particles:',Np
    print 'Number of halos:',Nh
