import numpy as np
import struct

infile = "halo.z=00.0000"

fp = open(infile,"rb")
FileContent = fp.read()
nh, = struct.unpack('q',FileContent[0:8])
x_p = np.zeros(nh)
y_p = np.zeros(nh)
z_p = np.zeros(nh)
x_v = np.zeros(nh)
y_v = np.zeros(nh)
z_v = np.zeros(nh)
M200h = np.zeros(nh)
R200h = np.zeros(nh)
M500h = np.zeros(nh)
R500h = np.zeros(nh)
N200h = np.zeros(nh)
N500h = np.zeros(nh)
M200c = np.zeros(nh)
R200c = np.zeros(nh)
M500c = np.zeros(nh)
R500c = np.zeros(nh)
N200c = np.zeros(nh)
N500c = np.zeros(nh)

halo_struct_size = 72

for i in range(nh):
  first_byte = 8+i*halo_struct_size
  x_p[i], = struct.unpack('f',FileContent[first_byte:first_byte+4])
  y_p[i], = struct.unpack('f',FileContent[first_byte+4:first_byte+8])
  z_p[i], = struct.unpack('f',FileContent[first_byte+8:first_byte+12])

  x_v[i], = struct.unpack('f',FileContent[first_byte+12:first_byte+16])
  y_v[i], = struct.unpack('f',FileContent[first_byte+16:first_byte+20])
  z_v[i], = struct.unpack('f',FileContent[first_byte+20:first_byte+24])

  N200h[i], = struct.unpack('i',FileContent[first_byte+24:first_byte+28])
  M200h[i], = struct.unpack('f',FileContent[first_byte+28:first_byte+32])
  R200h[i], = struct.unpack('f',FileContent[first_byte+32:first_byte+36])

  N500h[i], = struct.unpack('i',FileContent[first_byte+36:first_byte+40])
  M500h[i], = struct.unpack('f',FileContent[first_byte+40:first_byte+44])
  R500h[i], = struct.unpack('f',FileContent[first_byte+44:first_byte+48])

  N200c[i], = struct.unpack('i',FileContent[first_byte+48:first_byte+52])
  M200c[i], = struct.unpack('f',FileContent[first_byte+52:first_byte+56])
  R200c[i], = struct.unpack('f',FileContent[first_byte+56:first_byte+60])

  N500c[i], = struct.unpack('i',FileContent[first_byte+60:first_byte+64])
  M500c[i], = struct.unpack('f',FileContent[first_byte+64:first_byte+68])
  R500c[i], = struct.unpack('f',FileContent[first_byte+68:first_byte+72])

  print i, '\t', x_p[i], '\t', y_p[i], '\t', z_p[i] \
         , '\t', x_v[i], '\t', y_v[i], '\t', z_v[i] \
         , '\t', N200h[i], '\t', M200h[i], '\t', R200h[i] \
         , '\t', N500h[i], '\t', M500h[i], '\t', R500h[i] \
         , '\t', N200c[i], '\t', M200c[i], '\t', R200c[i] \
         , '\t', N500c[i], '\t', M500c[i], '\t', R500c[i] 

