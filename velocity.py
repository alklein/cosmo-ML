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
#        print 'particle:',particle[ID],'halo:',halo[ID],'separation:',sep,'R200h:',halo[R200h]
        if separation(particle, halo) < 2*halo[R200h]: 
#            print 'binding particle',particle[ID],'to halo',halo[ID]
            return halo[ID]
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
particles = np.loadtxt('particles.txt')
halos = np.loadtxt('halos.txt')
#------------- Rescale Particle Positions --------------#
if opts.verbose: print '\t position conversion factor:',x_unit
particles[:,X] *= x_unit
particles[:,Y] *= x_unit
particles[:,Z] *= x_unit
#------------ Rescale Particle Velocities --------------#
if opts.verbose: print '\t velocity conversion factor:',vel_unit
particles[:,VX] *= vel_unit
particles[:,VY] *= vel_unit
particles[:,VZ] *= vel_unit
#--------------- Order and Cut by Mass -----------------#
""" TODO: cut by mass """
halos = halos[halos[:,M200h].argsort()][::-1]
print '\t highest mass of a halo: 10 ^',np.log10(halos[0][M200h])
print '\t lowest mass of a halo: 10 ^',np.log10(halos[-1][M200h])
print '\t range of halo Xs:',min(halos[:,X]),max(halos[:,X])
print '\t range of particle Xs:',min(particles[:,X]),max(particles[:,X])
print '\t range of halo VXs:',min(halos[:,VX]), max(halos[:,VX])
print '\t range of particle VXs:',min(particles[:,VX]), max(particles[:,VX])
print '\t range of halo virial radii:',min(halos[:,R200h]), max(halos[:,R200h])
#-------------- Assign Particles to Halos --------------#
if opts.verbose: print '\n ... Assigning particles to 2 heaviest halos ... \n'
halo_IDs = [get_halo_ID(particle, halos[:2]) for particle in particles]
print '\t number of assigned particles:',len([i for i in halo_IDs if i > -1])
print '\t 2*radius of most massive halo:',2*halos[0][R200h]
print '\t min sep of any particle to this halo:',min([separation(particle, halos[0]) for particle in particles])
print '\t position of 1st particle:',particles[0][X],particles[0][Y],particles[0][Z]
print '\t position of 1st halo:',halos[0][X],halos[0][Y],halos[0][Z]
print '\t separation:',separation(particles[0], halos[0])

#particles = np.column_stack((particles, halo_IDs))

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
    some_particles = np.array([p for p in particles if (p[X] > 0 and p[X] < 5 and p[Y] > 0 and p[Y] < 5 and p[Z] > 0 and p[Z] < 5)])
    some_halos = np.array([p for p in halos if (p[X] > 0 and p[X] < 5 and p[Y] > 0 and p[Y] < 5 and p[Z] > 0 and p[Z] < 5)])
    some_halos = some_halos[::-1] # does halo assignment change if the halos are considered in reverse order?
    print '\t first particle:'
    print some_particles[0]
    print '\t first halo:'
    print some_halos[0]
    print '\t # of particles in 5 x 5 square:',len(some_particles)
    print '\t # of halo centers in 5 x 5 square:',len(some_halos)
    print '\t range of halo R200h values:',min(some_halos[:,R200h]),max(some_halos[:,R200h])
    halo_IDs = [get_halo_ID(particle, some_halos) for particle in some_particles]
    print '\t number of assigned particles from 5 x 5 square:',len([i for i in halo_IDs if i > -1])
    print
    new_particles = []
    for i in range(len(some_particles)):
        cur_row = [val for val in some_particles[i]]
        cur_row.append(get_halo_ID(cur_row, some_halos))
        new_particles.append(cur_row)
    some_particles = np.array(new_particles)
    labeled_particles = some_particles[ some_particles[:,-1] > -1]
    print '\t first labeled particle:'
    print labeled_particles[0]
    print '\t first halo center:'
    first_halo = some_halos[0][ID]
    print first_halo
    print '\t its purported mass:',some_halos[0][M200h]
    sample = np.array([p for p in labeled_particles if p[-1] == first_halo])
    print '\t # of particles assigned to it:',len(sample)
    print '\t their X coords:'
    print sample[:,X]
    print '\t their Y coords:'
    print sample[:,Y]

    data = [[some_particles[:,Y]], [some_particles[:,Z]]]
    labels = ['Y coordinate', 'Z coordinate']
    titles = ['Particle Density']
#    vis(data, labels, titles)

    figure(0)
    plot(some_particles[:,X], some_particles[:,Y], '.')
#    plot(labeled_particles[:,X], labeled_particles[:,Y], '.')
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


