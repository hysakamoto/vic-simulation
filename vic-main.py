""" Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

# Begin simulation

from dolfin import *
import pdb

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 2
parameters["num_threads"] = 2
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}


dt = 1.0
dt_const = Constant(dt)

m_num = 10 # number of mesh in each dim
# Create mesh and define function space
mesh = UnitCubeMesh(m_num, m_num, m_num)
V = VectorFunctionSpace(mesh, "Lagrange", 1)

# ## boundary mesh
# bmesh = BoundaryMesh(mesh, "exterior")
# bbtree = BoundaryBoxTree()
# bbtree.build(bmesh)
# distane_to_boundary = MeshFunction("double", mesh, 0)

## Boundary Conditions

# Create mesh function over cell facets
exterior_facet_domains = FacetFunction("size_t", mesh)

# Mark Neumann boundaries
class TopBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14   # tolerance for coordinate comparisons
        return on_boundary and abs(x[2] - 1.0) < tol

Gamma_T    = TopBoundary()
exterior_facet_domains.set_all(1)
Gamma_T.mark(exterior_facet_domains, 0)
ds_neumann = ds[exterior_facet_domains]

# Mark boundary subdomians
# top    = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)
bottom = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)

# Define Dirichlet boundary (x = 0 or x = 1)
# d_top    = Expression(("0.0", "0.0", "0.5"))
d_bottom = Expression(("0.0", "0.0", "0.0"))
# bc_top    = DirichletBC(V, d_top, top)
bc_bottom = DirichletBC(V, d_bottom, bottom)
# bcs = [bc_top, bc_bottom]
bcs = [bc_bottom]


# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration

## Contact
zeta = 0.02
pen_contact=1000
# penetration = Expression(("0.0","0.0",\
#                         "-pen_contact*((x[2])-(1.0-t*zeta))"), \
#                         t=0.0, zeta=zeta, pen_contact=pen_contact)
penetration = Expression(("0.0","0.0",\
                        "((x[2])-(1.0-t*zeta))"), \
                        t=0.0, zeta=zeta)
# penetration = Expression("1.0")
# pen_contact = 100
# contact_trac = penetration*pen_contact


Body  = Constant((0.0,  0.0, 0.0))  # Body force per unit volume
# Trac  = ((0.0,  0.0, contact_trac))  # Traction force on the boundary
# Trac = penetration


# Define function G such that G \cdot n = g
class NormalTraction(Expression):
    def __init__(self, mesh):
        self.mesh = mesh
    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        if ((ufc_cell.local_facet)!=-1):
            n = cell.normal(ufc_cell.local_facet)
            g = 1.0 # sin(5*x[0])   
            values[0] = g*n[0]
            values[1] = g*n[1]
            values[2] = g*n[2]
        else:
            values[0] = 0.0
            values[1] = 0.0
            values[2] = 0.0
    def value_shape(self):
        return (3,)


Trac = NormalTraction(mesh)






# Kinematics
d = u.geometric_dimension() # dimension
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient - seems like underyling coordinates is initial
F = variable(F)
C = F.T*F                   # Right Cauchy-Green tensor

# virtual Kinematics
ddotF = grad(v)
ddotE = 1.0/2.0 * (ddotF.T*F + F.T*ddotF)

# Invariants of deformation tensors
J  = det(F)
Ic = tr(C)
IIIc = det(C)

# Elasticity parameters
Ee, nu = 10.0, 0.4
mu, lmbda = Constant(Ee/(2*(1 + nu))), Constant(Ee*nu/((1 + nu)*(1 - 2*nu)))

# Strain energy density (compressible neo-Hookean model)
psi = (mu/2.0)*(Ic - 3) - mu*ln(J) + (lmbda/2.0)*(ln(J))**2
# psi = lmbda/2*(tr(E)**2) + mu*tr(E*E)

# # Strain energy density (nearly incompressible neo-hookean model)
# kappa = mu*10
# C_hat = (IIIc)**(-1.0/3.0)*C
# psi = mu/2.0*(tr(C_hat)-3) + 1.0/2.0*kappa*ln(J)**2.0
# # psi = mu/2*(tr(C_hat)-3) + kappa/2*(J-1)**2

# PK1 stress tensor
P = diff(psi,F)
# PK2 stress tensor
S = inv(F)*P

## Viscoelasticity!!!!

## Initial conditions
u_1 = Function(V)
u_1.interpolate(Constant((0.0, 0.0, 0.0)))

## viscoelastic parameters
gamma = Constant(20.0)
tau = Constant(10.0)

## Kinematics
F_1    = I + grad(u_1)             # Deformation gradient
F_1    = variable(F_1)             # Make F a variable for tensor differentiations
C_1    = F_1.T*F_1                   # Right Cauchy-Green tensor
invF_1 = inv(F_1)
# Invariants of deformation tensors
J_1    = det(F_1)
Ic_1   = tr(C_1)
IIIc_1 = det(C_1)
# Strain energy density (nearly incompressible neo-hookean model)
C_hat_1 = (IIIc_1)**(-1.0/3.0)*C_1
psi_1 = (mu/2.0)*(Ic_1 - 3) - mu*ln(J_1) + (lmbda/2.0)*(ln(J_1))**2
# psi_1   = mu/2.0*(tr(C_hat_1)-3) + 1.0/2.0*kappa*ln(J_1)**2.0
# PK1 stress tensor
P_1 = diff(psi_1,F_1)
# PK2 stress tensor
S_1 = inv(F_1)*P_1

H = exp(-dt_const/tau)*S_1 + (1-exp(-dt_const/tau))*(S-S_1/(dt_const/tau))    
Sc = S+gamma*H

# Compute residual
# R = tr(Sc*ddotE.T)*dx - dot(Body,v)*dx - (dot(Trac,v)+(-pen_contact*u[2])*v[2])*ds_neumann(0)
N=FacetNormal(mesh) # initial unit normal
n = F*N/dot(F*N,F*N) # current unit normal
R = tr(Sc*ddotE.T)*dx - dot(Body,v)*dx - dot(n,v)*ds_neumann(0)
# R = inner(S, ddotE)*dx - inner(B, v)*dx - inner(T, v)*ds

# Compute Jacobian of F
Jac = derivative(R, u, du)

set_log_level(DEBUG)
problem = NonlinearVariationalProblem(R, u, bcs=bcs, J=Jac)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["linear_solver"] = "bicgstab"
solver.parameters["newton_solver"]["preconditioner"] = "ilu"

file = File("results/displacement.pvd");

t = 0.0
tn = 0
T_total = 10.0

# Save solution in VTK format
file << (u, t);

### Run Simulation
while t<T_total:
    print 'time = ', t    
    penetration.t = t

    # solve
    solver.solve()

    # update
    t    += dt
    tn   += 1
    assign(u_1,u)

    # Save solution in VTK format
    file << (u, t);

    # plot(p_1, title = "pressure", axes=True, interactive = True)


