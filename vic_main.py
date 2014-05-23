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

# sim_name = 'result/'
# u_max, Eus, Eps = vic_func.vic_sim( sim_name, 
#                                     m_num, p_order, dt, T_total, max_it, omega,
#                                     Ee, nu, gamma, tau, perm,
#                                     top_trac, body_force )

m_nums = [1,2,4,8,16]
m_nums = [1,2]
sim_basename = 'mesh/'
for i in range(len(m_nums)):
    m_num = m_nums[i]
    sim_name = sim_basename + str(m_num)
    u_max, Eus, Eps, Eus2, Eps2  \
        = vic_func.vic_sim( sim_name, 
                            m_num, p_order, dt, T_total, max_it, omega,
                            Ee, nu, gamma, tau, perm,
                            top_trac, body_force )
    errors_u.append(sum(Eus)*dt)
    errors_p.append(sum(Eps)*dt)
    u_maxmax.append(max(u_max))
    errors_u2.append(sum(Eus2)*dt)
    errors_p2.append(sum(Eps2)*dt)

# T_total = 10.0
# max_its = [1,2,4,8,16,32,64]
# m_num = 8
# sim_basename = 'time_cn'
# for i in range(len(max_its)):
#     max_it = max_its[i]
#     dt = T_total/float(max_it)
#     sim_name = sim_basename + str(max_it)
#     u_max, Eus, Eps = vic_func.vic_sim( sim_name,
#                                         m_num, p_order, dt, T_total, max_it, omega,
#                                         Ee, nu, gamma, tau, perm,
#                                         top_trac, body_force )
#     errors_u.append(sum(Eus)*dt)
#     errors_p.append(sum(Eps)*dt)
#     u_maxmax.append(max(u_max))

print errors_u
print errors_p
# print u_maxmax

# plt.loglog(max_its, errors_u,'-o'); 
# plt.show()

print 'u-convergence rate: '+ str((np.log(errors_u[0])-np.log(errors_u[-1]))/(np.log(m_nums[0])-np.log(m_nums[-1])))
print 'p-convergence rate: '+ str((np.log(errors_p[0])-np.log(errors_p[-1]))/(np.log(m_nums[0])-np.log(m_nums[-1])))

## ending time
time_end = time.time()

## Output convergence result
with open(sim_basename+'conv.txt', 'w') as f:
    f.write(str(errors_u))
    f.write(str(errors_p))

    f.write('\nu-convergence rate: '+ str((np.log(errors_u[0])-np.log(errors_u[-1]))/(np.log(m_nums[0])-np.log(m_nums[-1]))))
    f.write('\np-convergence rate: '+ str((np.log(errors_p[0])-np.log(errors_p[-1]))/(np.log(m_nums[0])-np.log(m_nums[-1]))))

    f.write('\nelapsed time: ' + str(time_start-time_end))


print('\nu-convergence rate: '+ str((np.log(errors_u2[0])-np.log(errors_u2[-1]))/(np.log(m_nums[0])-np.log(m_nums[-1]))))
print('\np-convergence rate: '+ str((np.log(errors_p2[0])-np.log(errors_p2[-1]))/(np.log(m_nums[0])-np.log(m_nums[-1]))))
