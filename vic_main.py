""" Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

from dolfin import *
from vic_func import *
import pdb
import matplotlib.pyplot as plt

import vic_func
reload(vic_func)

# Begin simulation
m_num = 4
p_order = 1
T_total =40.0
dt     = 4.0       # time step
omega = 0.5
gamma = 0.0
tau = 1.0
Ee = 4.0
nu = 0.3
perm = 0.01

u_max = vic_func.vic_sim( m_num, p_order, dt, T_total, omega, Ee, nu, gamma, tau, perm )

plt.plot(u_max, '-o')
plt.show()
