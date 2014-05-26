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
print "Generating manufactured solutions from MATLAB output..."
execfile("conveq.py")


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
nu      = 0.3
perm    = 1.0
top_trac    = (0.0,0.0,0.0)
body_force = (0.0,0.0,0.0)

max_its = [2,4,8,16,32]
dts = [T_total/mit for mit in max_its]

T_total = 10.0
dt      = 1    # time step
max_it  = int(T_total/dt)


errors_u = []
errors_p = []
u_maxmax = []
errors_u2 = []
errors_p2 = []


######################### TEST SIMULATION #########################

sim_name = 'result/'
max_m_nums = 16
n_err_comp = max_m_nums/m_num+1
u_max, Eus, Eps, Eus2, Eps2\
    = vic_func.vic_sim( sim_name, 
                        m_num, p_order, dt, T_total, max_it, omega,
                        Ee, nu, gamma, tau, perm,
                        top_trac, body_force,
                        n_err_comp)

# m_nums = [1,2,4,8,16]
# sim_basename = 'mesh/'
# for i in range(len(m_nums)):
#     m_num = m_nums[i]
#     sim_name = sim_basename + str(m_num)
#     u_max, Eus, Eps, Eus2, Eps2  \
#         = vic_func.vic_sim( sim_name, 
#                             m_num, p_order, dt, T_total, max_it, omega,
#                             Ee, nu, gamma, tau, perm,
#                             top_trac, body_force )
#     errors_u.append(sum(Eus)*dt)
#     errors_p.append(sum(Eps)*dt)
#     u_maxmax.append(max(u_max))
#     errors_u2.append(sum(Eus2)*dt)
#     errors_p2.append(sum(Eps2)*dt)

# T_total = 10.0
# max_its = [2,4,8,16,32,64]
# max_it_num = max(max_its)
# sim_basename = 'time/'
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

#     errors_u.append(sum(Eus))
#     errors_p.append(sum(Eps))


#     ## Output convergence result
#     with open(sim_basename+'conv.txt', 'w') as f:
#         f.write(str(errors_u))
#         f.write('\n')
#         f.write(str(errors_p))
#         f.write('\n')
#         f.write(str(errors_u2))
#         f.write('\n')
#         f.write(str(errors_p2))

#     # u_maxmax.append(max(u_max))


# print '\nu-convergence rate: '+ str((np.log(errors_u[0])-np.log(errors_u[-1]))/(np.log(m_nums[0])-np.log(m_nums[-1])))
# print '\np-convergence rate: '+ str((np.log(errors_p[0])-np.log(errors_p[-1]))/(np.log(m_nums[0])-np.log(m_nums[-1])))

# print('\nu-convergence rate: '+ str((np.log(errors_u2[0])-np.log(errors_u2[-1]))/(np.log(m_nums[0])-np.log(m_nums[-1]))))
# print('\np-convergence rate: '+ str((np.log(errors_p2[0])-np.log(errors_p2[-1]))/(np.log(m_nums[0])-np.log(m_nums[-1]))))


## ending time
time_end = time.time()

# ## Output convergence result
# with open(sim_basename+'conv.txt', 'w') as f:
#     f.write(str(errors_u))
#     f.write('\n')
#     f.write(str(errors_p))
#     f.write('\n')
#     f.write(str(errors_u2))
#     f.write('\n')
#     f.write(str(errors_p2))

#     f.write('\nelapsed time: ' + str(time_end-time_start))


# ## plotting
# plt.loglog(max_its, errors_p2,'-o'); 
# plt.show()



