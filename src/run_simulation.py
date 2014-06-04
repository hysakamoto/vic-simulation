""" Run a Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

from vic_func import vic_sim
import time


def write_params(sim_name, mat_params, sim_params):
    import os
    try: 
        os.makedirs(sim_name)
    except OSError:
        print "Could not create a directory"
        return -1

    try:
        f = open(sim_name+'/params.txt', 'w')
    except IOError:
        print "Could not open file!"
        return -1

    with f:
        f.write('mat_params = ' + str(mat_params) + '\n')
        f.write('sim_params = ' + str(sim_params) + '\n')
        f.close()

    return 0


def runSim( base_name, mat_params, sim_params ):

    sim_name = base_name + str(sim_params['m_num']) \
               + '-' + str(sim_params['max_it']) + '/'

    # only the rank-0 process checks the starting condition
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank==0:
        print '**** Running simulation '+sim_name + ' ****'
        exist = write_params(sim_name, mat_params, sim_params)
    else:
        exist = None
    exist = comm.bcast(exist, root=0)

    # if the parameter file is already created, don't run the simulation
    if exist == -1:
        return 0

    ## starting time
    time_start = time.time()

    dt = sim_params['T_total']/float(sim_params['max_it'])
    vic_sim( sim_name, 
             sim_params['m_num'], 
             sim_params['p_order'], 
             dt, 
             sim_params['T_total'], 
             sim_params['max_it'], 
             sim_params['omega'],
             mat_params['Ee'], 
             mat_params['nu'], 
             mat_params['gamma'], 
             mat_params['tau'], 
             mat_params['perm'], 
             mat_params['top_trac'], 
             mat_params['body_force'] )

    ## ending time
    time_end = time.time()

    return time_end - time_start



