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
T_total = 10.0
dt     = 0.5       # time step
omega = 0.5

u_max = vic_func.vic_sim( m_num, p_order, dt, T_total, omega )

plt.plot(u_max)
plt.show()
