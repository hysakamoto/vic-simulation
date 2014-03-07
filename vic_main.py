""" Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

from dolfin import *
from vic_func import *
import pdb

reload(vic_func)

# Begin simulation
m_num = 10
p_order = 1
dt     = 0.2       # time step

vic_sim( m_num, p_order, dt )
