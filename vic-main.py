""" Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

# Begin simulation

from dolfin import *
import pdb

# DEBUG
set_log_level(DEBUG)

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
dim = mesh.topology().dim() 

##  Function Spaces
porder = 2
Pu = VectorFunctionSpace(mesh, "Lagrange", porder-1) # space for displacements
Pp = FunctionSpace(mesh, "Lagrange", porder-1)     # space for pressure
V = MixedFunctionSpace([Pu,Pp])                    # mixed space

# Functions
dup = TrialFunction(V)    # Incremental displacement-pressure
wq  = TestFunction(V)     # Test function
up  = Function(V)         # Displacement-pressure from previous iteration
u, p = split(up)          # Function in each subspace to write the functional
w, q = split(wq)          # Test Function split

# Loads
B  = Constant((0.0,  0.0, 0.0))  # Body force per unit volume
Trac = Expression(("0.0","0.0","1.0")) # Traction force on the boundary
g_bar = Constant(0.0)            # Normal flux


## Boundary Condition

# Create mesh function over cell facets
exterior_facet_domains = FacetFunction("size_t", mesh)

# Mark Neumann boundaries
class TopBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14   # tolerance for coordinate comparisons
        return on_boundary and abs(x[2] - 1.0) < tol

Gamma_T = TopBoundary()
exterior_facet_domains.set_all(1)
Gamma_T.mark(exterior_facet_domains, 0)
ds_neumann = ds[exterior_facet_domains]

# Mark Dirichlet boundaries
top    = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)
bottom = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)

# Assign Dirichlet boundaries (x = 0 or x = 1)
d_top    = Expression(("0.0", "0.0", "0.0"))
d_bottom = Expression(("0.0", "0.0", "0.0"))
bc_top    = DirichletBC(V.sub(0), d_top, top)
bc_bottom = DirichletBC(V.sub(0), d_bottom, bottom)

p_top    = Expression("0.0")
p_bottom = Expression("1.0")
bc_ptop = DirichletBC(V.sub(1), p_top, top)
bc_pbottom = DirichletBC(V.sub(1), p_bottom, bottom)

bcs = [bc_bottom, bc_ptop]

## Initial conditions
up_1 = Function(V)
u_1, p_1 = split(up_1)
u_1 = interpolate(Constant((0.0, 0.0, 0.0)), Pu)
p_1 = interpolate(Constant(0.0), Pp)

## Kinematics
I = Identity(dim)             # Identity tensor
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

## Potential Energy
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

## time steps
dt = 0.2 # time step
T_toal = 1.0
t = dt
tn=0

# Save initial conditions in VTK format
dfile = File("displacement_%d.pvd"%tn);
dfile << up.sub(0);
pfile = File("pressure_%d.pvd"%tn);
pfile << up.sub(1);

omega = 1.0

## Run Simulation
while t<T_toal:
    print 'time = ', t    

    # Compute residual
    v = (u-up_1.sub(0))/(dt*omega)
    R = (inner(S, ddotE) + dot(J*(K_perm*invF.T*grad(p)), invF.T*grad(q)))*dx \
        - p*J*inner(ddotF, invF.T)*dx \
        + (q*J*inner(grad(v),invF.T))*dx \
        + (inner(B,w))*dx - (inner(Trac,w))*ds_neumann(0) \
        + (inner(g_bar,q))*ds_neumann(1)

    # Compute Jacobian of F
    Jac = derivative(R, up, dup)
    
    # Set up the problem
    problem = NonlinearVariationalProblem(R, up, bcs=bcs, J=Jac )
    solver = NonlinearVariationalSolver(problem)
    solver.parameters["newton_solver"]["linear_solver"] = "bicgstab"
    solver.parameters["newton_solver"]["preconditioner"] = "ilu"

    # solve
    solver.solve()

    # update
    t += dt
    tn+=1
    # v_1 = (up.sub(0)-up_1)/(dt*omega) - (1-omega)/omega*v_1
    up_1.assign(up)

    # Save solution in VTK format
    dfile = File("displacement_%d.pvd"%tn);
    dfile << up.sub(0);
    pfile = File("pressure_%d.pvd"%tn);
    pfile << up.sub(1);

    # plot(p_1, title="pressure", axes=True, interactive = True)
    

# Plot and hold solution
# plot(u, mode = "displacement", title="displacement", axes=True, interactive = True)
plot(p, title="pressure", axes=True, interactive = True)
