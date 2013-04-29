#!/usr/bin/python

from time import *
from constants import *

def get_filename(left_bound, right_bound):
    return str(int(left_bound)) + '_' + str(int(right_bound)) + '.txt'

def get_filename_dict(ID_bounds):
    filename_dict = {}
    for (L, R) in ID_bounds:
        for i in range(int(L), int(R)+1): #inclusive range
            filename_dict[i] = get_filename(L, R)
    return filename_dict

def get_ID_bounds(halo_file, max_chunk_sz = None, chunk_dir = 'pchunks'):
    halos = np.loadtxt(halo_file)
    if not max_chunk_sz: max_chunk_sz = halos[0][N200a]

    cur_sum = halos[0][N200a]
    left_bound = halos[0][ID]
    right_bound = halos[0][ID]
    
    ID_bounds = []
    for h in halos:
        cur_id, num_p = h[ID], h[N200a]
        if cur_sum + num_p > max_chunk_sz:
            ID_bounds.append((left_bound, right_bound))
            left_bound = cur_id
            right_bound = cur_id
            cur_sum = num_p
        else: 
            cur_sum += num_p
            right_bound += 1
    if (max_chunk_sz == halos[0][N200a]): return ID_bounds[1:]
    else: return ID_bounds

def make_chunks(ID_bounds):
    i = 0 
    (left_bound, right_bound) = ID_bounds[i]
    cur_chunk = []
    for line in open('particles.txt'):
        particle = [float(val) for val in line.split()]
        cur_ID = particle[1]
        if (cur_ID >= left_bound and cur_ID <= right_bound):
            cur_chunk.append(particle)
        else: 
            # time to save chunk
            os.chdir(chunk_dir)
            filename = get_filename(left_bound, right_bound)
            print '\nsaving array of length',len(cur_chunk),'to file:',filename
            print '\t Progress:',100*(i+1.)/len(ID_bounds),'%'
            assert len(cur_chunk) > 0
            np.savetxt(filename, cur_chunk)
            os.chdir('..')
            # start new chunk
            i += 1
            (left_bound, right_bound) = ID_bounds[i]
            cur_chunk = [] # reset to blank chunk
            if (cur_ID >= left_bound and cur_ID <= right_bound):
                cur_chunk.append(particle)

def fetch(halo_ID, filename_dict, chunk_dir = 'pchunks'):
    P = np.loadtxt(chunk_dir + '/' + str(filename_dict[halo_ID]))
    P = np.column_stack((P[:,0], P[:,2:], P[:,1]))
    return P[P[:,H_ID] == halo_ID]

def test(halo_file, filename_dict):
    fetch_times = []
    length_errors = []
    for h in np.loadtxt(halo_file):
        starttime = clock()
        P = fetch(h[ID], filename_dict)
        fetch_times.append(clock() - starttime)
        if (len(P) != h[N200a]): length_errors.append(len(P) - int(h[N200a]))
        for p in P: assert p[H_ID] == h[ID]
        print '\tHalo',h[ID],'verified. Current avg fetch time:',np.average(fetch_times)
    return fetch_times, length_errors

if (__name__ == '__main__'):
    halo_file = 'halos.txt'
    chunk_dir = 'pchunks'
    print '\nfinding ID bounds...'
    ID_bounds = get_ID_bounds(halo_file)
    print '\nnumber of chunks:',len(ID_bounds)
    
    print '\ngenerating dictionary...'
    filename_dict = get_filename_dict(ID_bounds)

#    print '\nmaking chunk files...'
#    make_chunks(ID_bounds)

    print '\ntesting chunk files... \n'
    fetch_times, length_errors = test(halo_file, filename_dict)
    print '\nunexpected length found',len(length_errors),'times'

    figure(0)
    hist(fetch_times, bins=50)
    xlabel('Time to fetch chunk (s)', fontsize=24)
    ylabel('Number', fontsize=24)
    title('Particle Retrieval Efficiency Distribution', fontsize=30)
    show()



