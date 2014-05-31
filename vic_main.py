""" Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

from dolfin import *
from vic_func import *
import pdb
import matplotlib.pyplot as plt
import time

import vic_func
reload(vic_func)

## starting time
time_start = time.time()

## Generate the manufactured solutions, bcs, ics

# from mpi4py import MPI
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# if rank == 0:
#     print "Generating manufactured solutions from MATLAB output..."
#     execfile("conveq.py")
# comm.barrier()
# comm.Disconnect()
print "Generating manufactured solutions from MATLAB output..."
execfile("conveq.py")

# import sys
# sys.exit()

# Begin simulation
m_num   = 8
p_order = 1
T_total = 20.0
dt      = 0.5    # time step
max_it  = 20
omega   = 0.5    # forward Euler: 0, backward Euler: 1, Crank-Nicholson: 0.5 
gamma   = 0.0    # 0 for no viscoelasticity
tau     = 1.0
Ee      = 100.0
nu      = 0.25
perm    = 1.0
top_trac    = (0.0,0.0,0.0)
body_force = (0.0,0.0,0.0)


T_total = 10.0
dt      = 1    # time step
max_it  = int(T_total/dt)


errors_u = []
errors_p = []
u_maxmax = []
errors_u2 = []
errors_p2 = []


######################### TEST SIMULATION #########################

# sim_name = 'result/'
# max_m_nums = 16
# n_err_comp = max_m_nums/m_num+1
# u_max, Eus, Eps, Eus2, Eps2\
#     = vic_func.vic_sim( sim_name, 
#                         m_num, p_order, dt, T_total, max_it, omega,
#                         Ee, nu, gamma, tau, perm,
#                         top_trac, body_force,
#                         n_err_comp)

######################### MESH SIMULATION #########################

sim_basename = 'mesh/'
T_total = 64
max_it = 16
dt = T_total/float(max_it)

m_nums = [4,8,16]
for i in range(len(m_nums)):
    m_num = m_nums[i]
    sim_name = sim_basename + str(m_num)
    u_max, Eus, Eps, Eus2, Eps2  \
        = vic_func.vic_sim( sim_name, 
                            m_num, p_order, dt, T_total, max_it, omega,
                            Ee, nu, gamma, tau, perm,
                            top_trac, body_force, 1 )


######################### TIME SIMULATION #########################

# sim_basename = 'time/'
# m_num = 16

# T_total = 10.0
# max_its = [2,4,8,16,32,64]
# max_it_num = max(max_its)
# for i in range(len(max_its)):
#     max_it = max_its[i]
#     n_err_comp = max_it_num/max_it+1 # number of error computation per time step
#     n_err_comp = 2
#     dt = T_total/float(max_it)
#     sim_name = sim_basename + str(max_it)

#     u_max, Eus, Eps, Eus2, Eps2\
#         = vic_func.vic_sim( sim_name, 
#                             m_num, p_order, dt, T_total, max_it, omega,
#                             Ee, nu, gamma, tau, perm,
#                             top_trac, body_force,
#                             n_err_comp)




## ending time
time_end = time.time()



