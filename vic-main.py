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

# Create mesh and define function space
m_num = 6 # number of mesh in each dim
mesh = UnitCubeMesh(m_num, m_num, m_num)
d = u.geometric_dimension() # dimension

# Define function spaces
porder = 2
Pu = VectorFunctionSpace(mesh, "Lagrange", porder) # space for displacements
Pp = FunctionSpace(mesh, "Lagrange", porder-1)     # space for pressure
V = MixedFunctionSpace([Pu,Pp])                    # mixed space

# Define functions
dup = TrialFunction(V)    # Incremental displacement-pressure
wq  = TestFunction(V)     # Test function
up  = Function(V)         # Displacement-pressure from previous iteration
u, p = split(up)          # Function in each subspace to write the functional
w, q = split(wq)          # Test Function split

# Loads
B  = Constant((0.0,  0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0,  0.0, 0.0))  # Traction force on the boundary
g_bar = Constant(0.0)            # Normal flux

# Mark boundary subdomians
top    = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)
bottom = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)

# Define Dirichlet boundary (x = 0 or x = 1)
d_top    = Expression(("0.0", "0.0", "0.5"))
d_bottom = Expression(("0.0", "0.0", "0.0"))
bc_top    = DirichletBC(V.sub(0), d_top, top)
bc_bottom = DirichletBC(V.sub(0), d_bottom, bottom)
bcs = [bc_top, bc_bottom]
## can predcribe dirichlet pressure condition

# Kinematics
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient - seems like underyling coordinates is initial
F = variable(F)             # Make F a variable for tensor differentiations
C = F.T*F                   # Right Cauchy-Green tensor
invF = inv(F)

# virtual Kinematics
ddotF = grad(w)
ddotE = 1.0/2.0 * (ddotF.T*F + F.T*ddotF)

# Invariants of deformation tensors
J  = det(F)
Ic = tr(C)
IIIc = det(C)

# Elasticity parameters
Ee, nu = 10.0, 0.45
mu, lmbda = Constant(Ee/(2*(1 + nu))), Constant(Ee*nu/((1 + nu)*(1 - 2*nu)))

# Permeability
perm = 0.1
K_perm = perm*I

# Strain energy density (compressible neo-Hookean model)
# psi = (mu/2.0)*(Ic - 3) - mu*ln(J) + (lmbda/2.0)*(ln(J))**2

# Strain energy density (nearly incompressible neo-hookean model)
kappa = mu*10
C_hat = (IIIc)**(-1.0/3.0)*C
psi = mu/2.0*(tr(C_hat)-3) + 1.0/2.0*kappa*ln(J)**2.0

# PK1 stress tensor
P = diff(psi,F)
# PK2 stress tensor
S = inv(F)*P

# Total potential energy
Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

# Compute residual (no velocity for now)
R = (inner(S, ddotE) + dot(J*(K_perm*invF.T*grad(p)), invF.T*grad(q)))*dx \
    - p*J*inner(ddotF, invF.T)*dx  \
    - (inner(T,w))*ds + (inner(g_bar,q))*ds 
# R = (inner(S, ddotE) + dot(J*(K_perm*invF.T*grad(p)), invF.T*grad(q)))*dx \
#     - p*J*inner(ddotF, invF.T)*dx + (q*J*inner(grad(v),invF.T))*dx \
#     - (inner(T,w))*ds + (inner(g_bar,q))*ds 

# Compute Jacobian of F
Jac = derivative(R, up, dup)

set_log_level(DEBUG)
problem = NonlinearVariationalProblem(R, up, bcs=bcs, J=Jac)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["linear_solver"] = "bicgstab"
solver.parameters["newton_solver"]["preconditioner"] = "ilu"
solver.solve()

# Save solution in VTK format
file = File("displacement.pvd");
file << up.sub(0);

# Plot and hold solution
plot(u, mode = "displacement", axes=True, interactive = True)
plot(p, title="pressure", axes=True, interactive = True)
