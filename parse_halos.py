#!/usr/bin/env python                                                                                           

"""                                                                                                             

@file parse_halos.py
@brief Script to extract halos from binary simulation output file
@author Andrea Klein       <alklein@alumni.stanford.edu>

Usage: python parse_halos.py > halos.txt
Note: the infile halo.z=00.0000 must be in the same directory.

"""

import numpy as np
import struct

infile = "halo.z=00.0000"
fp = open(infile,"rb")
FileContent = fp.read()

# number of halos
nh, = struct.unpack('q',FileContent[0:8])

# halo ID
ID = np.zeros(nh)

# redshift (blank)
z = np.zeros(nh)
# distance [comoving Mpc/h] (blank)
d = np.zeros(nh)
# phi angle [deg] (blank)
phi = np.zeros(nh)
# theta angle [deg] (blank)
theta = np.zeros(nh)

# 3D position [comoving Mpc/h]
x_p = np.zeros(nh)
y_p = np.zeros(nh)
z_p = np.zeros(nh)

# 3D velocity [proper km/s]
x_v = np.zeros(nh)
y_v = np.zeros(nh)
z_v = np.zeros(nh)

# num of particles s.t. <density> = 200*cosmic density (a)
N200a = np.zeros(nh)
# mass [Msolar/h] corresponding to...
M200a = np.zeros(nh)
# ...radius [proper Mpc/h]
R200a = np.zeros(nh)
# velocity dispersion [proper km/s]
s200a = np.zeros(nh)

# num of particles s.t. <density> = 500*cosmic density (a)
N500a = np.zeros(nh)
M500a = np.zeros(nh)
R500a = np.zeros(nh)
s500a = np.zeros(nh)

# num of particles s.t. <density> = 200*critical density (c)
N200c = np.zeros(nh)
M200c = np.zeros(nh)
R200c = np.zeros(nh)
s200c = np.zeros(nh)

# num of particles s.t. <density> = 500*critical density (c)
N500c = np.zeros(nh)
M500c = np.zeros(nh)
R500c = np.zeros(nh)
s500c = np.zeros(nh)

# max circular velocity [proper km/s]
vcmax = np.zeros(nh)
# mass [Msolar/h] corresponding to vcmax
Mcmax = np.zeros(nh)
# radius [proper Mpc/h] corresponding to vcmax
rcmax = np.zeros(nh)

# 1(8) + 29(4)
halo_struct_size = 124

for i in range(nh):
  first_byte = 8+i*halo_struct_size
  
  """
  ID is allocated 8 bytes in the binary file, but you only need 4,
  so bytes 4-8 of each struct are skipped. The halos are in numerical 
  order, so I output i+1 rather than this value.
  """
  ID[i], = struct.unpack('i',FileContent[first_byte:first_byte+4])

  z[i], = struct.unpack('f',FileContent[first_byte+8:first_byte+12])
  d[i], = struct.unpack('f',FileContent[first_byte+12:first_byte+16])
  phi[i], = struct.unpack('f',FileContent[first_byte+16:first_byte+20])
  theta[i], = struct.unpack('f',FileContent[first_byte+20:first_byte+24])

  x_p[i], = struct.unpack('f',FileContent[first_byte+24:first_byte+28])
  y_p[i], = struct.unpack('f',FileContent[first_byte+28:first_byte+32])
  z_p[i], = struct.unpack('f',FileContent[first_byte+32:first_byte+36])

  x_v[i], = struct.unpack('f',FileContent[first_byte+36:first_byte+40])
  y_v[i], = struct.unpack('f',FileContent[first_byte+40:first_byte+44])
  z_v[i], = struct.unpack('f',FileContent[first_byte+44:first_byte+48])

  N200a[i], = struct.unpack('i',FileContent[first_byte+48:first_byte+52])
  M200a[i], = struct.unpack('f',FileContent[first_byte+52:first_byte+56])
  R200a[i], = struct.unpack('f',FileContent[first_byte+56:first_byte+60])
  s200a[i], = struct.unpack('f',FileContent[first_byte+60:first_byte+64])
  N500a[i], = struct.unpack('i',FileContent[first_byte+64:first_byte+68])
  M500a[i], = struct.unpack('f',FileContent[first_byte+68:first_byte+72])
  R500a[i], = struct.unpack('f',FileContent[first_byte+72:first_byte+76])
  s500a[i], = struct.unpack('f',FileContent[first_byte+76:first_byte+80])

  N200c[i], = struct.unpack('i',FileContent[first_byte+80:first_byte+84])
  M200c[i], = struct.unpack('f',FileContent[first_byte+84:first_byte+88])
  R200c[i], = struct.unpack('f',FileContent[first_byte+88:first_byte+92])
  s200c[i], = struct.unpack('f',FileContent[first_byte+92:first_byte+96])

  N500c[i], = struct.unpack('i',FileContent[first_byte+96:first_byte+100])
  M500c[i], = struct.unpack('f',FileContent[first_byte+100:first_byte+104])
  R500c[i], = struct.unpack('f',FileContent[first_byte+104:first_byte+108])
  s500c[i], = struct.unpack('f',FileContent[first_byte+108:first_byte+112])

  vcmax[i], = struct.unpack('f',FileContent[first_byte+112:first_byte+116])
  Mcmax[i], = struct.unpack('f',FileContent[first_byte+116:first_byte+120])
  rcmax[i], = struct.unpack('f',FileContent[first_byte+120:first_byte+124])

  # currently omitting blank properties z, d, phi, and theta:
  # , '\t', z[i],   '\t', d[i],   '\t', phi[i], '\t', theta[i] \
  print i+1 \
           , '\t', x_p[i], '\t', y_p[i], '\t', z_p[i] \
           , '\t', x_v[i], '\t', y_v[i], '\t', z_v[i] \
           , '\t', N200a[i], '\t', M200a[i], '\t', R200a[i], '\t', s200a[i] \
           , '\t', N500a[i], '\t', M500a[i], '\t', R500a[i], '\t', s500a[i] \
           , '\t', N200c[i], '\t', M200c[i], '\t', R200c[i], '\t', s200c[i] \
           , '\t', N500c[i], '\t', M500c[i], '\t', R500c[i], '\t', s500c[i] \
           , '\t', vcmax[i], '\t', Mcmax[i], '\t', rcmax[i]

