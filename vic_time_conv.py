""" Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

from dolfin import *
from vic_func import *
import pdb
import numpy as np

import vic_func
reload(vic_func)

# Begin simulation
m_num = 10
p_order = 1
dt     = 1.0       # time step
omega  = 0.5       # time integration weight

# up = []
# for i in range(6):
#     up.append(vic_sim( m_num, p_order, 1.0/2**i, omega ))



## Create mesh and define function space
mesh = UnitCubeMesh(m_num, m_num, m_num)
##  Function Spaces
Pu = VectorFunctionSpace(mesh, "Lagrange", p_order)  # space for displacements
Pp = FunctionSpace(mesh, "Lagrange", p_order)        # space for pressure


# ux,uy,uz,p = [],[],[],[]
# for i in range(6):
#     ux.append ( up[i].vector().array()[Pu.sub(0).dofmap().dofs()])
#     uy.append (up[i].vector().array()[Pu.sub(1).dofmap().dofs()])
#     uz.append (up[i].vector().array()[Pu.sub(2).dofmap().dofs()])
#     p.append  (up[i].vector().array()[Pp.dofmap().dofs()])
#     np.save('sol_vector_%d'%i,up[i].vector().array())

upva = []
for i in range(6):
    upva.append(np.load('sol_vector_%d.npy'%i))


ux,uy,uz,p = [],[],[],[]
for i in range(6):
    ux.append ( upva[i][Pu.sub(0).dofmap().dofs()])
    uy.append (upva[i][Pu.sub(1).dofmap().dofs()])
    uz.append (upva[i][Pu.sub(2).dofmap().dofs()])
    p.append  (upva[i][Pp.dofmap().dofs()])


ux = np.array(ux)
uy = np.array(uy)
uz = np.array(uz)
p = np.array(p)

er_asse = []
dts = []
for i in range(6-1):
    e_Pp = Function(Pp)
    e_Pp.vector()[:] = p[i][:]-p[5][:]
    error_x = e_Pp**2*dx
    er_asse.append (sqrt(assemble(error_x, mesh=Pu.mesh())))
    dts.append(1.0/2**i)

import matplotlib.pyplot as plt

plt.plot(dts,er_asse)
plt.xscale('log')
plt.yscale('log')
plt.show()

