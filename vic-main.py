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

m_num = 10 # number of mesh in each dim
# Create mesh and define function space
mesh = UnitCubeMesh(m_num, m_num, m_num)
V = VectorFunctionSpace(mesh, "Lagrange", 1)

# Mark boundary subdomians
top    = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)
bottom = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)

# Define Dirichlet boundary (x = 0 or x = 1)
d_top    = Expression(("0.0", "0.0", "0.5"))
d_bottom = Expression(("0.0", "0.0", "0.0"))
bc_top    = DirichletBC(V, d_top, top)
bc_bottom = DirichletBC(V, d_bottom, bottom)
bcs = [bc_top, bc_bottom]

# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration
B  = Constant((0.0,  0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0,  0.0, 0.0))  # Traction force on the boundary

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
Ee, nu = 10.0, 0.45
mu, lmbda = Constant(Ee/(2*(1 + nu))), Constant(Ee*nu/((1 + nu)*(1 - 2*nu)))

# Strain energy density (compressible neo-Hookean model)
# psi = (mu/2.0)*(Ic - 3) - mu*ln(J) + (lmbda/2.0)*(ln(J))**2

# Strain energy density (nearly incompressible neo-hookean model)
kappa = Constant(mu*10)
C_hat = (IIIc)**(-1.0/3.0)*C
psi = mu/2.0*(tr(C_hat)-3) + 1.0/2.0*kappa*ln(J)**2.0

# Total potential energy
Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

# Compute residual
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
plot(u, mode = "displacement", axes=True, interactive = True)
