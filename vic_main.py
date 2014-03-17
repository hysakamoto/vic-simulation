""" Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

from dolfin import *
from vic_func import *
import pdb
import matplotlib.pyplot as plt

import vic_func
reload(vic_func)

# Begin simulation
m_num = 10
p_order = 1
T_total =4.0
dt     = 2.0       # time step
omega = 0.5
gamma = 20.0
tau = 1.0
Ee = 10.0
nu = 0.45
perm = 1.0

u_max = vic_func.vic_sim( m_num, p_order, dt, T_total, omega, Ee, nu, gamma, tau, perm )

plt.plot(u_max, '-o')
plt.show()
