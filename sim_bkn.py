from run_simulation import runsim

## OMEGA = 1.0 :: BACKWARD EULER METHOD

base_name = 'bkn/'


## NOTE: Race condition may occur under parallel environment, so just set it up before the simulation.
## Generate the manufactured solutions, bcs, ics
# from mpi4py import MPI
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# if rank == 0:
#     print "Generating manufactured solutions from MATLAB output..."
#     execfile("convert_equations.py")
# comm.barrier()



# material parameters
mat_params = {'gamma'      : 0.0, 
              'tau'        : 1.0,
              'Ee'         : 100.0,
              'nu'         : 0.25,
              'perm'       : 1.0,
              'top_trac'   : (0.0,0.0,0.0),
              'body_force' : (0.0,0.0,0.0)
}



# simulation parameters
sim_params = {'omega'   : 1.0,   # forward:0, backward: 1, C-N: 0.5 
              'T_total' : 10.0,
              'max_it'  : 64,
              'm_num'   : 32,
              'p_order' : 1
          }


## convergence parameters
# mesh_refinement = [1,2,4,8,16,32]
mesh_refinement = [16,32]
max_iterations = [1,2,4,8,16,32,64]

max_mer = 32
max_mit = 64


# mesh convergence
sim_params['m_num'] = max_mer
sim_params['max_it'] = max_mit
for mer in mesh_refinement:
    sim_params['m_num'] = mer
    runsim(base_name, mat_params, sim_params)

# time convergence
sim_params['m_num'] = max_mer
sim_params['max_it'] = max_mit
for mit in max_iterations:
    sim_params['max_it'] = mit
    runsim(base_name, mat_params, sim_params)

# set the parameters back
sim_params['m_num'] = max_mer
sim_params['max_it'] = max_mit

