""" Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

from dolfin import *
from vic_func import *
import pdb
import matplotlib.pyplot as plt

import vic_func
reload(vic_func)

# Begin simulation
m_num   = 4
p_order = 1
T_total = 20.0
dt      = 1.0    # time step
max_it  = 20
omega   = 0.5    # forward Euler: 0, backward Euler: 1, Crank-Nicholson: 0.5 
gamma   = 0.0    # 0 for no viscoelasticity
tau     = 1.0
Ee      = 100.0
nu      = 0.45
perm    = 1.0
top_trac    = (0.0,0.0,0.0)
body_force = (0.0,0.0,0.0)

max_its = [2,4,8,16,32]
dts = [T_total/mit for mit in max_its]
errors_u = []
errors_p = []
u_maxmax = []

# for i in range(len(dts)):
#     dt = dts[i]
#     max_it = max_its[i]

#     u_max, Eus, Eps = vic_func.vic_sim( m_num, p_order, dt, T_total, max_it, omega,
#                                         Ee, nu, gamma, tau, perm,
#                                         top_trac, body_force )
#     errors_u.append(sum(Eus)*dt)
#     errors_p.append(sum(Eps)*dt)
#     u_maxmax.append(max(u_max))

T_total = 10.0
dt      = 1.0    # time step
max_it  = int(T_total/dt)

m_nums = [1,2,4,8,16]
m_nums = [4]
for i in range(len(m_nums)):
    m_num = m_nums[i]
    u_max, Eus, Eps = vic_func.vic_sim( m_num, p_order, dt, T_total, max_it, omega,
                                        Ee, nu, gamma, tau, perm,
                                        top_trac, body_force )
    errors_u.append(sum(Eus)*dt)
    errors_p.append(sum(Eps)*dt)
    u_maxmax.append(max(u_max))


# plt.plot(E1s, '-o')
# plt.show()


## analytical solutions

#  X*(-1/((t - 20)^2/400 - 2))^(1/2) - X
#  Y*(-1/((t - 20)^2/400 - 2))^(1/2) - Y
#                -Z*((t - 20)^2/400 - 1)
#  p = -(((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2


