## Test Simulation

import sys
sys.path.append("src") 

from  dolfin import *

from run_simulation import runSim

base_name = 'test/'

set_log_level(PROGRESS)
# PETScOptions.set("pc_type", "ml") 
# PETScOptions.set("log_summary", "LOG")  

## Generate the manufactured solutions, bcs, ics
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    print "Generating manufactured solutions from MATLAB output..."
    execfile("src/convert_equations.py")
comm.barrier()

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


# mesh convergence
sim_params['m_num'] = 8
sim_params['max_it'] = 10



## Remove directory before start
if rank == 0:
    from subprocess import call
    call(["rm", "-r", "test"])
comm.barrier()



runSim(base_name, mat_params, sim_params)

# u_1.vector()[:] =  u_1.vector()[:]+omega*up.sub(0,deepcopy=True).vector()*dt

