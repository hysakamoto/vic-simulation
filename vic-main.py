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
mesh = UnitCubeMesh(16, 16, 16)
V = VectorFunctionSpace(mesh, "Lagrange", 1)

# Mark boundary subdomians
top =   CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)
bottom = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)

# Define Dirichlet boundary (x = 0 or x = 1)
d_top = Expression(("0.0", "0.0", "0.1"))
d_bottom = Expression(("0.0", "0.0", "0.0"))

bct = DirichletBC(V, d_top, top)
bcb = DirichletBC(V, d_bottom, bottom)
# bcs = [bct, bcb]
bcs = [bcb]

# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration
B  = Constant((0.0,  0.0, 0.1))  # Body force per unit volume
T  = Constant((0.0,  0.0, 0.0))  # Traction force on the boundary

# Kinematics
d = u.geometric_dimension() # dimension
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient - seems like underyling coordinates is initial
F = variable(F)
C = transpose(F)*F                   # Right Cauchy-Green tensor
# C = variable(C)             # make it a variable for differntiation
E = (C-I)/2
# E = variable(E)
invC = inv(C)               # inverse of C

# virtual Kinematics
ddotF = grad(v)
ddotE = 1.0/2.0 * (ddotF.T*F + F.T*ddotF)

# Invariants of deformation tensors
J  = det(F)
Ic = tr(C)
IIIc = det(C)
# IIIc = J**2

# Elasticity parameters
Ee, nu = 10.0, 0.45
mu, lmbda = Constant(Ee/(2*(1 + nu))), Constant(Ee*nu/((1 + nu)*(1 - 2*nu)))

# Stored strain energy density (compressible neo-Hookean model)
# psi = (mu/2.0)*(Ic - 3) - mu*ln(J) + (lmbda/2.0)*(ln(J))**2
# psi = lmbda/2*(tr(E)**2) + mu*tr(E*E) #  (material model)

# nearly incompressible neo-hookean sed function
# kappa = mu*100
kappa = Constant(100)
C_hat = (IIIc)**(-1.0/3.0)*C
psi = mu/2*(tr(C_hat)-3) + 1/2*kappa*(J-1)**2
# psi = mu/2*(tr(C_hat)-3) + 1/2*kappa*ln(J)**2

# PK1 stress tensor
# P = diff(psi,F)
# PK2 stress tensor
# S = inv(F)*P

# PK2 stress tensor
# S = 2*diff(psi, C)
# S = diff(psi,E)
# S = mu*(I-invC)+lmbda*ln(J)*invC

## nearly incompressible neo-hookean PK2
# S = mu * IIIc**(-1/3)*(I-1/3*Ic*invC) + kappa*(J-1)*J*invC

# Total potential energy
Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

# Compute residual
# R = tr(S*ddotE.T)*dx - dot(B,v)*dx - dot(T,v)*ds
# R = inner(S, ddotE)*dx - inner(B, v)*dx - inner(T, v)*ds
R = derivative(Pi, u, v)

# Compute Jacobian of F
Jac = derivative(R, u, du)

set_log_level(DEBUG)
problem = NonlinearVariationalProblem(R, u, bcs=bcs, J=Jac)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["linear_solver"] = "bicgstab"
solver.parameters["newton_solver"]["preconditioner"] = "ilu"
solver.solve()

# Save solution in VTK format
file = File("displacement.pvd");
file << u;

# Plot and hold solution
plot(u, mode = "displacement", interactive = True)
