#!/usr/bin/env python                                                                                           

"""                                                                                                             

@file velocity.py
@brief Script to learn velocity distributions within DM halos
@author Andrea Klein       <alklein@alumni.stanford.edu>

"""

__author__   = "Andrea Klein"

from PyML import *
from constants import *

#-------------------------------------------------------#
#---------------------- FUNCTIONS ----------------------#
#-------------------------------------------------------#

#----------------------- Plotting ----------------------#

# description: determines bins for use by vis function
# input:
     # data to 'visualize'
     # desired pixel density per axis
# output:
     # x and y values for use in plotting
     # bounds of plotting (to fit all data)
# effects: n/a
# notes: n/a

def bns_ext(data, res=400):
    xmins, xmaxs = map(min, data[0]), map(max, data[0])
    ymins, ymaxs = map(min, data[1]), map(max, data[1])
    xmin, xmax = min(xmins), max(xmaxs)
    ymin, ymax = min(ymins), max(ymaxs)
    xrange, yrange = np.linspace(xmin, xmax, res), np.linspace(ymin, ymax, res)
    return [yrange, xrange], [xmin, xmax, ymin, ymax]

# description: makes colored 'scatter' plots for side-by-side comparison of data
# input:
     # data = [[ x1, x2, ... ],[ y1, y2, ... ]]
     # labels = [ xlabel, ylabel ]
     # titles = [ title1, title2, ... ]
# output: n/a (makes graph)
# effects: n/a
# notes: does not take log of data (must use np.log beforehand for semilog or log-log vis)

def vis(data, labels, titles):
    i = 0
    while i < len(data[0]):
        figure(i)
        xs, ys = data[0][i], data[1][i]
        bns, ext = bns_ext(data)
        H, yedges, xedges = np.histogram2d(ys, xs, bins=bns)
        imshow(H, extent=ext, aspect='equal', cmap = cm.jet, interpolation='nearest', origin='lower', vmin=0.0, vmax=10.0)
        colorbar()
        rc('text', usetex=True)
        xlabel(labels[0], fontsize=20)
        ylabel(labels[1], fontsize=20)
        title(titles[i], fontsize=24)
        i += 1
    show()

#-------------------- Pre-Processing -------------------#

def separation(particle, halo):
    px, py, pz = particle[X], particle[Y], particle[Z]
    hx, hy, hz = halo[X], halo[Y], halo[Z]
    return ((px - hx)**2 + (py - hy)**2 + (pz - hz)**2)**.5

def get_halo_ID(particle, halos):
    for halo in halos:
        sep = separation(particle, halo) 
        if separation(particle, halo) < 2*halo[R200a]: return halo[ID]
    return -1

#--------------- Velocities (all in km/s) --------------#

def Vnet(VX, VY, VZ):
    return (VX**2 + VY**2 + VZ**2)**(.5)

def R_dot(X, Y, Z, VX, VY, VZ):
    return (X*VX + Y*VY + Z*VZ) / ((X**2 + Y**2 + Z**2)**.5)

def Theta_dot(X, Y, Z, VX, VY, VZ):
    return (VX*Y - X*VY) / ((X**2 + Y**2)**.5)

def Phi_dot(X, Y, Z, VX, VY, VZ):
    num = Z*(X*VX + Y*VY) - (X**2 + Y**2)*VZ
    denom = (X**2 + Y**2 + Z**2)*((X**2 + Y**2)**.5)
    return num/denom

#------------------- Particle Extraction -------------------#

def get_particles(particle_file, halo_ID):
    P = []
    for line in open(particle_file):
        particle = [float(val) for val in line.split()]
        cur_ID = particle[1]
        if cur_ID > halo_ID: 
            P = np.array(P)
            if len(P) > 0: return np.column_stack((P[:,0], P[:,2:], P[:,1]))
            else: return P
        if cur_ID == halo_ID: P.append(particle)
    if len(P) > 0: return np.column_stack((P[:,0], P[:,2:], P[:,1]))
    else: P = np.array(P)


#-----------------------------------------------------------#
#---------------------- PARSE OPTIONS ----------------------#
#-----------------------------------------------------------#

usage = 'Usage: %prog [options]'
description = 'This program learns the velocity distribution of dark matter subhalos.'
parser = OptionParser(usage=usage, description=description)

parser.add_option('-v','--verbose', default = True, action = 'store_true', help = 'display additional output (default: True)', dest='verbose')
parser.add_option('-p','--plot', default = False, action = 'store_true', help = 'display some relevant plots from AQA & VLII (default: False)', dest='make_plots')
parser.add_option('-s','--hist', default = False, action = 'store_true', help = 'display some relevant histograms from AQA & VLII (default: False)', dest='make_hists')
parser.add_option('-y','--hyplot', default = False, action = 'store_true', help = 'display some relevant plots from Hy\'s sims (default: False)', dest='hy_plots')

(opts, arguments) = parser.parse_args()
are_args = len(arguments)
if (not are_args):
    print 
    parser.print_help()

#-------------------------------------------------------#
#-------------------- GET AQ/VL DATA -------------------#
#-------------------------------------------------------#

if opts.verbose: print '\n ... Loading Aquarius and Via Lactea II ...'
AQ = np.loadtxt('AqA_NoMainHalo.txt')
VL = np.loadtxt('VL2_NoMainHalo.txt')
#--------------- Cut at Virial Radius ----------------#
AQ = np.array([sh for sh in AQ if sh[D_GC] < R_vir_AQ])
VL = np.array([sh for sh in VL if sh[D_GC] < R_vir_VL])
#--------------- X, Y, Z: kpc/s -> km/s --------------#
AQ[:,X:Z+1] = AQ[:,X:Z+1]*km_per_kpc
VL[:,X:Z+1] = VL[:,X:Z+1]*km_per_kpc
#---------------- Compute Velocities -----------------#
AQ_Vnet = [Vnet(sh[Vx], sh[Vy], sh[Vz]) for sh in AQ]
VL_Vnet = [Vnet(sh[Vx], sh[Vy], sh[Vz]) for sh in VL]
AQ_R_dot = [R_dot(sh[X], sh[Y], sh[Z], sh[Vx], sh[Vy], sh[Vz]) for sh in AQ]
VL_R_dot = [R_dot(sh[X], sh[Y], sh[Z], sh[Vx], sh[Vy], sh[Vz]) for sh in VL]
AQ_Theta_dot = [Theta_dot(sh[X], sh[Y], sh[Z], sh[Vx], sh[Vy], sh[Vz]) for sh in AQ]
VL_Theta_dot = [Theta_dot(sh[X], sh[Y], sh[Z], sh[Vx], sh[Vy], sh[Vz]) for sh in VL]
AQ_Phi_dot = [Phi_dot(sh[X], sh[Y], sh[Z], sh[Vx], sh[Vy], sh[Vz]) for sh in AQ]
VL_Phi_dot = [Phi_dot(sh[X], sh[Y], sh[Z], sh[Vx], sh[Vy], sh[Vz]) for sh in VL]

#-------------------------------------------------------#
#--------------------- GET HY DATA ---------------------#
#-------------------------------------------------------#

if opts.verbose: print '\n ... Loading Hy Trac\'s Simulations ... \n'
halos = np.loadtxt('halos.txt')
# TEMP
#num_part = 0
#for line in open('particles.txt'): num_part += 1
#print 'number of particles:',num_part
# TEMP: only considering particles assigned to 1st halo
particles = get_particles('particles.txt', halos[0][ID]) 
#--------------- Order and Cut by Mass -----------------#
halos = halos[halos[:,M200a].argsort()][::-1]

print '\t range of halo Xs:',min(halos[:,X]),max(halos[:,X])
print '\t range of particle Xs:',min(particles[:,X]),max(particles[:,X])
print '\t range of halo VXs:',min(halos[:,VX]), max(halos[:,VX])
print '\t range of particle VXs:',min(particles[:,VX]), max(particles[:,VX])


#-------------------------------------------------------#
#---------------------- MAKE PLOTS ---------------------#
#-------------------------------------------------------#

if opts.make_plots:

#-------------------- Aquarius A ---------------------#

    figure(0)
    semilogx(AQ[:,D_GC], AQ_Vnet, '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel('V (km/s)', fontsize=24)
    title('AQ', fontsize=30)
    savefig('cosmo_plot_0.pdf', format='pdf')
    
    figure(1)
    semilogx(AQ[:,D_GC], AQ_R_dot, '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel(r'$\mathrm{\dot{R}}$ (km/s)',fontsize=24)
    axes = gca()
    axes.set_ylim(-450, 450)
    title('AQ', fontsize=30)
    savefig('cosmo_plot_1.pdf', format='pdf')

    figure(2)
    semilogx(AQ[:,D_GC], AQ_Theta_dot, '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel(r'$\mathrm{\dot{\theta}}$ (rad/s)', fontsize=24)
    axes = gca()
    axes.set_ylim(-450, 450)
    title('AQ', fontsize=30)
    savefig('cosmo_plot_2.pdf', format='pdf')

    figure(3)
    semilogx(AQ[:,D_GC], AQ_Phi_dot, '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel(r'$\mathrm{\dot{\phi}}$ (rad/s)', fontsize=24)
    axes = gca()
    axes.set_ylim(-2e-15, 2e-15)
    title('AQ', fontsize=30)
    savefig('cosmo_plot_3.pdf', format='pdf')

    figure(4)
    semilogx(AQ[:,D_GC], AQ[:,Vx], '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel(r'$\mathrm{V_x}$ (km/s)', fontsize=24)
    axes = gca()
    axes.set_xlim(10**16, 10**20)
    axes.set_ylim(-600, 600)
    title('AQ', fontsize=30)
    savefig('cosmo_plot_4.pdf', format='pdf')

    figure(5)
    semilogx(AQ[:,D_GC], AQ[:,Vy], '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel(r'$\mathrm{V_y}$ (km/s)', fontsize=24)
    axes = gca()
    axes.set_xlim(10**16, 10**20)
    axes.set_ylim(-600, 600)
    title('AQ', fontsize=30)
    savefig('cosmo_plot_5.pdf', format='pdf')

    figure(6)
    semilogx(AQ[:,D_GC], AQ[:,Vz], '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel(r'$\mathrm{V_z}$ (km/s)', fontsize=24)
    axes = gca()
    axes.set_xlim(10**16, 10**20)
    axes.set_ylim(-600, 600)
    title('AQ', fontsize=30)
    savefig('cosmo_plot_6.pdf', format='pdf')
    
#-------------------- Via Lactea ---------------------#
    
    figure(10)
    semilogx(VL[:,D_GC], VL_Vnet, '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel('V (km/s)', fontsize=24)
    title('VL', fontsize=30)
    savefig('cosmo_plot_10.pdf', format='pdf')

    figure(11)
    semilogx(VL[:,D_GC], VL_R_dot, '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel(r'$\mathrm{\dot{R}}$ (km/s)',fontsize=24)
    axes = gca()
    axes.set_ylim(-450, 450)
    title('VL', fontsize=30)
    savefig('cosmo_plot_11.pdf', format='pdf')

    figure(12)
    semilogx(VL[:,D_GC], VL_Theta_dot, '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel(r'$\mathrm{\dot{\theta}}$ (rad/s)', fontsize=24)
    axes = gca()
    axes.set_ylim(-450, 450)
    title('VL', fontsize=30)
    savefig('cosmo_plot_12.pdf', format='pdf')

    figure(13)
    semilogx(VL[:,D_GC], VL_Phi_dot, '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel(r'$\mathrm{\dot{\phi}}$ (rad/s)', fontsize=24)
    axes = gca()
    axes.set_ylim(-2e-15, 2e-15)
    title('VL', fontsize=30)
    savefig('cosmo_plot_13.pdf', format='pdf')

    figure(14)
    semilogx(VL[:,D_GC], VL[:,Vx], '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel(r'$\mathrm{V_x}$ (km/s)', fontsize=24)
    axes = gca()
    axes.set_xlim(10**16, 10**20)
    axes.set_ylim(-600, 600)
    title('VL', fontsize=30)
    savefig('cosmo_plot_14.pdf', format='pdf')

    figure(15)
    semilogx(VL[:,D_GC], VL[:,Vy], '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel(r'$\mathrm{V_y}$ (km/s)', fontsize=24)
    axes = gca()
    axes.set_xlim(10**16, 10**20)
    axes.set_ylim(-600, 600)
    title('VL', fontsize=30)
    savefig('cosmo_plot_15.pdf', format='pdf')

    figure(16)
    semilogx(VL[:,D_GC], VL[:,Vz], '.')
    xlabel(r'$\mathrm{D_{GC}}$ (kpc)', fontsize=24)
    ylabel(r'$\mathrm{V_z}$ (km/s)', fontsize=24)
    axes = gca()
    axes.set_xlim(10**16, 10**20)
    axes.set_ylim(-600, 600)
    title('VL', fontsize=30)
    savefig('cosmo_plot_16.pdf', format='pdf')

if opts.make_hists:

    V_Max_AQ = 208.75
    V_Max_VL = 201.033

    figure(20)
    hist(AQ_Vnet, bins=50)
    axes = gca()
    axvline(x = V_Max_AQ, color = 'r', linewidth = 3)
    xlabel('V (km/s)', fontsize=24)
    ylabel('Number', fontsize=24)
    title('AQ', fontsize=30)

    figure(21)
    hist(VL_Vnet, bins=50)
    axes = gca()
    axvline(x = V_Max_VL, color = 'r', linewidth = 3)
    xlabel('V (km/s)', fontsize=24)
    ylabel('Number', fontsize=24)
    title('VL', fontsize=30)


if opts.hy_plots:

    print
    some_particles = particles
    some_halos = halos[0] 

    data = [[some_particles[:,X]], [some_particles[:,Y]]]
    labels = ['X coordinate', 'Y coordinate']
    titles = ['Particle Density']
    vis(data, labels, titles)

    figure(0)
    plot(some_particles[:,X], some_particles[:,Y], '.')
    plot(some_halos[:,X], some_halos[:,Y], 'x', color='r', markersize=20, markeredgewidth=5)
    xlabel('X coordinate', fontsize=24)
    ylabel('Y coordinate', fontsize=24)

    figure(1)

    plot(some_particles[:,Y], some_particles[:,Z], '.')
#    plot(labeled_particles[:,Y], labeled_particles[:,Z], '.')
    plot(some_halos[:,Y], some_halos[:,Z], 'x', color='r', markersize=20, markeredgewidth=5)
    xlabel('Y coordinate', fontsize=24)
    ylabel('Z coordinate', fontsize=24)

if opts.make_plots or opts.make_hists or opts.hy_plots: show()


