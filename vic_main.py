""" Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

from dolfin import *
from vic_func import *
import pdb


# Begin simulation
m_num = 10
p_order = 1
dt     = 0.2       # time step
omega = 0.5

vic_sim( m_num, p_order, dt, omega )
