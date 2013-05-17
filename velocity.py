#!/usr/bin/env python                                                                                           

"""                                                                                                             

@file velocity.py
@brief Script to learn velocity distributions within DM halos
@author Andrea Klein       <alklein@alumni.stanford.edu>

"""

__author__   = "Andrea Klein"

import chunk_manager

from PyML import *
from constants import *
from PyML.classifiers.svm import SVR

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

def V_rel(p, halos):
    h = halos[p[H_ID]]
    return (p[VX] - h[VX], p[VY] - h[VY], p[VZ] - h[VZ])

def Vnet_rel(p, halos):
    (pvx, pvy, pvz) = V_rel(p, halos)
    return Vnet(pvx, pvy, pvz)

def Vs_rel(p, halos):
    h = halos[p[H_ID]]
    return (p[VX] - h[VX], p[VY] - h[VY], p[VZ] - h[VZ])

def Pos_rel(p, halos):
    h = halos[p[H_ID]]
    return (p[X] - h[X], p[Y] - h[Y], p[Z] - h[Z])

#---------------- Data Creation / Extraction ----------------#

# efficient for high-mass halos only
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
    P = np.array(P)
    if len(P) > 0: return np.column_stack((P[:,0], P[:,2:], P[:,1]))
    else: P = np.array(P)

def get_training_data(halos, min_halo_ID, max_halo_ID, filename_dict, num_halos, p_per_halo=1):
    data = []
    which_halos = np.random.randint(min_halo_ID, max_halo_ID + 1, num_halos)
    for h in which_halos:
        P = chunk_manager.fetch(h, filename_dict)
        Vs = [Vnet_rel(p, halos) for p in P]
        for i in range(p_per_halo): 
            p = P[np.random.randint(0, high=len(P))]
            p = [val for val in p]
            p.append(np.average(Vs))
            p.append(np.std(Vs))
            data.append(p)
    return np.array(data)

def save_training_data(halos, min_halo_ID, max_halo_ID, filename_dict, num_halos, size_name, how_many = 5):
    for i in range(how_many):
        train_data = get_training_data(halos, min_halo_ID, max_halo_ID, filename_dict, num_halos)
        filename = size_name + '_' + str(i) + '.txt'
        print 'saving file:',filename
        np.savetxt(filename, train_data)

#---------------------- Custom ML ----------------------#

def RMSE(Y1s, Y2s):
    return (sum([(Y1s[i] - Y2s[i])**2 for i in range(len(Y1s))])/len(Y1s))**.5

def cross_validate(s, data):
    errs = []
    for i in range(len(data)):
        trainset = data[i]
        testsets = data[:i] + data[i+1:]
        s.train(trainset)
        err_row = []
        for testset in testsets:
            r = s.test(testset)
            #print 'r:',r
            orig_Ys = r.getGivenLabels()
            SVR_Ys = r.Y
            print 'orig Ys:',orig_Ys
            print 'SVR Ys:',SVR_Ys
            err = RMSE(SVR_Ys, orig_Ys)
            err_row.append(err)
        errs.append(err_row)
    return errs

def test_errs(s, trainset, testset):
    s.train(trainset)
    r = s.test(testset)
    orig_Ys = r.getGivenLabels()
    SVR_Ys = r.Y
    print 'orig Ys:',orig_Ys
    print 'SVR Ys:',SVR_Ys
    return r

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

if opts.verbose: print '\n ... Loading Aquarius A and Via Lactea II ...'
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

if opts.verbose: print '\n ... Loading Hy Trac\'s Simulations ... '
#--------------- Select, Order, and Cut Halos -----------------#
halo_file = 'halos.txt'
halos = np.loadtxt(halo_file)
halos = halos[halos[:,M200a].argsort()][::-1] # already done?
halos = halos[halos[:,N200a] > 1000]
#----------------------- Manage Chunks ------------------------#
ID_bounds = chunk_manager.get_ID_bounds(halo_file)
filename_dict = chunk_manager.get_filename_dict(ID_bounds)

#-------------------------------------------------------#
#------------------------- SVR -------------------------#
#-------------------------------------------------------#

data = np.loadtxt('enhanced_halos.txt')
data = data[:1000] # cut to med data set size
As = ((math.pi/8.)**.5)*data[:,-2]
avg_A = np.average(As)
rmse = sum([(As[i] - avg_A)**2 for i in range(len(As))])/len(As)
print '\tRMSE using avg as estimate:',rmse

if opts.verbose: print '\n ... Training Classifier ... \n'
#------------------- Select Training Data ---------------------#
train_file = 'med_train_max.data'
test_file = 'med_test_max.data'
#----------------------- Prepare Data -------------------------#
traindata = SparseDataSet(train_file, numericLabels = True)
testdata = SparseDataSet(test_file, numericLabels = True)
traindata.normalize(1)
testdata.normalize(1)

"""
size = 'med'
D0 = SparseDataSet(size + '_max_0.data', numericLabels = True)
D1 = SparseDataSet(size + '_max_1.data', numericLabels = True)
D2 = SparseDataSet(size + '_max_2.data', numericLabels = True)
D3 = SparseDataSet(size + '_max_3.data', numericLabels = True)
D4 = SparseDataSet(size + '_max_4.data', numericLabels = True)

D0.normalize(1)
D1.normalize(1)
D2.normalize(1)
D3.normalize(1)
D4.normalize(1)

K0 = SparseDataSet('large_max_0.data', numericLabels = True)
K0.normalize(1)
"""
#Cs = [2e-5, 2e-4, 2e-3, 2e-2, 2e-1, 2e0, 2e1, 2e2, 2e3, 2e4, 2e5, 2e7, 2e9] #, 2e11, 2e13, 2e15]
Cs = [2e-5, 2e-4, 2e-3, 2e-2, 2e-1, 2e0, 2e1, 2e2, 2e3, 2e4, 2e5, 2*10**6.5, 2*10**7, 2*10**7.5, 2*10**8, 2*10**8.5] #, 2*10**9.]
#----------------------- Linear Kernel ------------------------#
lin_errs = []
for c in Cs:
    s = SVR(C = c)
    s.train(traindata)
    errs = s.cv(testdata, 4)
    lin_errs.append(errs)
    #errs = test_errs(s, K0, D0)
    #errs = cross_validate(s, [K0, D1])
    #for line in errs: print errs
    #s.train(traindata)
    #r = s.test(testdata)
    #given = r.getGivenLabels()
    #Ys = r.Y
    #print 'test labels:',Ys
    #lin_errs.append(r)

print 'Cs:'
print Cs
print 'linear errors:'
print lin_errs
# TODO: divide data into 5 sets per dataset size; write own CV function; plot expected vs. actual Ys
exit(0) # TEMP
#---------------------- Gaussian Kernel -----------------------#
Gs = [2e-15, 2e-13, 2e-11, 2e-9, 2e-7, 2e-5, 2e-3, 2e-1] #, 2e1, 2e3]
gauss_errs = []
for c in Cs:
    err_row = []
    for g in Gs:
        traindata.attachKernel('gaussian', gamma = g)
        s = SVR(C = c)
        s.train(traindata)
        r = s.cv(testdata, 4)
        err_row.append(r)
    gauss_errs.append(err_row)
#---------------------- Polynomial Kernel ----------------------#
poly_errs = []
traindata.attachKernel('polynomial')
for c in Cs:
    s = SVR(C = c)
    s.train(traindata)
    r = s.cv(testdata, 4)
    poly_errs.append(err_row)
#-------------------------- Results ----------------------------#
print
print 'Cs:'
print Cs
print 'linear errors:'
print lin_errs
print 'poly errors:'
print poly_errs
print
print 'minimum error with linear kernel:',min(lin_errs)
print 'minimum error with gaussian kernel:',min([min(row) for row in gauss_errs])
print 'minimum error with deg-2 poly kernel:',min([min(row) for row in poly_errs])
exit(0)

#-------------------------------------------------------#
#---------------------- WORKSPACE ----------------------#
#-------------------------------------------------------#

save_training_data(halos, 1, 10000, filename_dict, 10, 'tiny')

light_halos = halos[halos[:,M200a] < 10**14]
print 'first light halo:',light_halos[0]
print 'ID:',light_halos[0][ID]
print 'mass: 10 ^',np.log10(light_halos[0][M200a])

i = 0
h = halos[i]
while np.log10(h[M200a] > 13): 
    i += 1
    h = halos[i]
print 'first halo w/ mass below 10^13:',halos[i][ID]

Ns = [N for N in halos[:,N200a] if N > 1000]
Ms = [N*4.5*10**9 for N in Ns]
Vs = [Vnet_rel(p, halos) for p in particles]
Rs = [((p[X] - halos[which][X])**2 + (p[Y] - halos[which][Y])**2 + (p[Z] - halos[which][Z])**2)**.5 for p in particles]
Xrel, Yrel, Zrel = [], [], []
VXrel, VYrel, VZrel = [], [], []
for p in particles:
    (xx, yy, zz) = Pos_rel(p, halos)
    (vx, vy, vz) = Vs_rel(p, halos)
    Xrel.append(xx)
    Yrel.append(yy)
    Zrel.append(zz)
    VXrel.append(vx)
    VYrel.append(vy)
    VZrel.append(vz)
Rdot_rel = [R_dot(Xrel[i], Yrel[i], Zrel[i], VXrel[i], VYrel[i], VZrel[i]) for i in range(len(particles))]

#print 'number of halos with > 1000 particles:',len(Ns)
#print 'number of particles in those halos:',sum(Ns)
#print 'max mass: 10 ^',np.log10(max(Ms))
#print 'min mass: 10 ^',np.log10(min(Ms))
#print
p0 = particles[0]
h0 = halos[which]
print 'halo under consideration:',h0[ID]
print 'its # particles:',h0[N200a]
print 'its mass: 10 ^',np.log10(h0[M200a])

#figure(0)
#semilogx(halos[:,M200a], halos[:,vcmax], '.')
#xlabel('Mass (Msun)', fontsize=24)
#ylabel('Vcmax (km/s)', fontsize=24)
"""
figure(1)
hist(np.log10(halos[:,N200a]), bins=50, log=True)
xlabel('Log(Number of Particles)', fontsize=24)
title('Distribution of Halo Particle Counts', fontsize=30)
axes = gca()
axvline(x = 3, color = 'r', linewidth = 3)

print
print 'len Rs:',len(Rs)
print 'len Vs:',len(Vs)
print 'range Rs:',min(Rs),max(Rs)

data = [[np.array(Rs)], [np.log10(Vs)]]
labels = ['Radial Position (mpc/h)', 'Log(Vnet) (km/s)']
titles = ['Velocity as a function of radius']
#vis(data, labels, titles)

figure(2)
hist(Vs, bins=50)
xlabel('Vnet (km/s)', fontsize=24)
title('Distribution of Net Particle Velocities', fontsize=30)
axes = gca()
axvline(x = halos[which][vcmax], color = 'r', linewidth = 3)

figure(3)
plot(Rs, np.log10(Vs), '.')
xlabel('Radial Position (mpc/h)', fontsize=24)
ylabel('Log(Vnet) (km/s)', fontsize=24)
title('Velocity as a function of radius', fontsize=30)
axes = gca()
axvline(x = halos[which][rcmax], color = 'r', linewidth = 3)
axvline(x = halos[which][R200a], color = 'g', linewidth = 3)

figure(4)
plot(Rs, Rdot_rel, '.')
xlabel('Radial Position (mpc/h)', fontsize=24)
ylabel(r'$\mathrm{\dot{R}}$ (km/s)', fontsize=24)
title('Radial velocity component as a function of radius', fontsize=24)

show()
"""
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
#    vis(data, labels, titles)

    figure(0)
    plot(some_particles[:,X], some_particles[:,Y], '.')
    plot([halos[0][X]], [halos[0][Y]], 'x', color='r', markersize=20, markeredgewidth=5)
#    plot(some_halos[:,X], some_halos[:,Y], 'x', color='r', markersize=20, markeredgewidth=5)
    xlabel('X coordinate', fontsize=24)
    ylabel('Y coordinate', fontsize=24)

#    figure(1)

#    plot(some_particles[:,Y], some_particles[:,Z], '.')
#    plot(labeled_particles[:,Y], labeled_particles[:,Z], '.')
#    plot(some_halos[:,Y], some_halos[:,Z], 'x', color='r', markersize=20, markeredgewidth=5)
#    xlabel('Y coordinate', fontsize=24)
#    ylabel('Z coordinate', fontsize=24)

if opts.make_plots or opts.make_hists or opts.hy_plots: show()


