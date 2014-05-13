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
dt      = 4.0    # time step
omega   = 0.5    # forward Euler: 0, backward Euler: 1, Crank-Nicholson: 0.5 
gamma   = 0.0    # 0 for no viscoelasticity
tau     = 1.0
Ee      = 100.0
nu      = 0.45
perm    = 1.0
top_trac    = (0.0,0.0,0.0)
body_force = (0.0,0.0,0.0)

u_max = vic_func.vic_sim( m_num, p_order, dt, T_total, omega,
                          Ee, nu, gamma, tau, perm,
                          top_trac, body_force )

plt.plot(u_max, '-o')
plt.show()



 # ((t/200 - 1/10)/pow((pow((t - 20),2)/400 - 2),2) - (t/200 - 1/10)/(pow((-1/(pow((t - 20),2)/400 - 2)),(3/2))*pow((pow((t - 20),2)/400 - 2),2)))*(Z - Z*(pow((t - 20),2)/400 - 1))

