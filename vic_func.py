""" Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

# Begin simulation

from dolfin import *
import numpy as np
import pdb
from newton_solve import newton_solver

# DEBUG
# set_log_level(DEBUG)

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 2
parameters["num_threads"] = 2
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}


def vic_sim( m_num, p_order, dt, T_total, omega, Ee, nu, gamma, tau, perm ):
    """ 
    run the VIC simulation
    """

    ## avoid recompilation
    dt_const = Constant(dt)

    ## Create mesh and define function space
    mesh = UnitCubeMesh(m_num, m_num, m_num)
    dim  = mesh.topology().dim() 

    ##  Function Spaces
    Pu = VectorFunctionSpace(mesh, "Lagrange", p_order)  # space for displacements
    Pp = FunctionSpace(mesh, "Lagrange", p_order)        # space for pressure
    V  = MixedFunctionSpace([Pu,Pp])                    # mixed space

    ## Functions
    dup  = TrialFunction(V)    # Incremental displacement-pressure
    wq   = TestFunction(V)     # Test function
    up   = Function(V)         # Displacement-pressure from previous iteration
    u, p = split(up)           # Function in each subspace to write the functional
    w, q = split(wq)           # Test Function split

    ## Loads
    B     = Constant((0.0,  0.0, 0.0))  # Body force per unit volume
    Trac  = Constant((0.0,  0.0, 2.0)) # Traction force on the boundary
    g_bar = Constant(0.0)            # Normal flux

    ## Boundary Conditions

    # Create mesh function over cell facets
    # exterior_facet_domains = FacetFunction("size_t", mesh)
    exterior_facet_domains = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    
    # Mark Neumann boundaries
    class TopBoundary(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14   # tolerance for coordinate comparisons
            return on_boundary and abs(x[2] - 1.0) < tol

    Gamma_T    = TopBoundary()
    exterior_facet_domains.set_all(0)
    Gamma_T.mark(exterior_facet_domains, 1)
    ds_neumann = ds[exterior_facet_domains]

    # Mark Dirichlet boundaries
    def side(x, on_boundary):
        tol = 1e-14
        return on_boundary \
            and (abs(x[0])<tol \
                 or (abs(x[0]-1)<tol) \
                 or (abs(x[1])<tol) \
                 or (abs(x[1]-1)<tol))

    def top(x, on_boundary):
        tol = 1e-14
        return on_boundary and (abs(x[2]-1.0)<tol)

    def bottom(x, on_boundary):
        tol = 1e-14
        return on_boundary and (abs(x[2])<tol)

    
    # top    = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)
    # bottom = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)

    # Assign Dirichlet boundaries (x = 0 or x = 1)
    d_top     = Expression(("0.0", "0.0", "0.0"))
    d_bottom  = Expression(("0.0", "0.0", "0.0"))

    bc_top    = DirichletBC(V.sub(0), d_top, top)
    bc_bottom = DirichletBC(V.sub(0), d_bottom, bottom)

    p_top = Expression("0.0")
    p_bottom   = Expression("1.0")
    p_side = Expression("0.0")

    bc_ptop    = DirichletBC(V.sub(1), p_top, top)
    bc_pbottom = DirichletBC(V.sub(1), p_bottom, bottom)
    bc_pside = DirichletBC(V.sub(1), p_side, side)

    bcs = [bc_bottom, bc_pbottom]

    ## Initial conditions
    up_1   = Function(V)         # Displacement-pressure from previous iteration
    u_1, p_1 = split(up_1)       # Function in each subspace to write the functional
    assign (up_1.sub(0), interpolate(Constant((0.0, 0.0, 0.0)),Pu))
    # assign (up_1.sub(1), interpolate(Expression('1.0-x[2]'),Pp))
    assign (up_1.sub(1), interpolate(Expression('0'),Pp))

    ## Kinematics
    I    = Identity(dim)           # Identity tensor
    F    = I + grad(u)             # Deformation gradient
    F    = variable(F)             # Make F a variable for tensor differentiations
    C    = F.T*F                   # Right Cauchy-Green tensor
    invF = inv(F)

    # virtual Kinematics
    ddotF = grad(w)
    ddotE = 1.0/2.0 * (ddotF.T*F + F.T*ddotF)

    # Invariants of deformation tensors
    J    = det(F)
    Ic   = tr(C)
    IIIc = det(C)

    # Elasticity parameters
    mu, lmbda = Constant(Ee/(2*(1 + nu))), Constant(Ee*nu/((1 + nu)*(1 - 2*nu)))

    # Permeability
    K_perm = Constant(np.eye(3)*perm) # avoid recompilation

    ## Potential Energy
    # Strain energy density (compressible neo-Hookean model)
    # psi = (mu/2.0)*(Ic - 3) - mu*ln(J) + (lmbda/2.0)*(ln(J))**2

    # Strain energy density (nearly incompressible neo-hookean model)
    kappa = mu*10
    C_hat = (IIIc)**(-1.0/3.0)*C
    psi   = mu/2.0*(tr(C_hat)-3) + 1.0/2.0*kappa*ln(J)**2.0

    # PK1 stress tensor
    P = diff(psi,F)
    # PK2 stress tensor
    S = inv(F)*P

    ## time steps
    t      = 0.0
    tn     = 0         # time step number

    ## Viscoelasticity!!!!

    ## H-value
    HS = TensorFunctionSpace(mesh, "Lagrange", p_order )
    H_1 = Function(HS)
    H_1.interpolate(Constant(np.zeros([3,3])))
    
    ## viscoelastic parameters
    gamma_const = Constant(gamma)
    tau_const = Constant(tau)

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
    psi_1   = mu/2.0*(tr(C_hat_1)-3) + 1.0/2.0*kappa*ln(J_1)**2.0
    # PK1 stress tensor
    P_1 = diff(psi_1,F_1)
    # PK2 stress tensor
    S_1 = invF_1*P_1

    # S_1 = Function(HS)
    # S_1.interpolate(Constant(np.zeros([3,3])))

    H = exp(-dt_const/tau_const)*H_1 + (1-exp(-dt_const/tau_const))*(S-S_1)/(dt_const/tau_const)
    Sc = S+gamma_const*H

    ## Poroelasticity!!!!
    omega_const = Constant(omega)
    v_1 = Function(Pu)
    v_1.interpolate(Constant((0.0, 0.0, 0.0)))
    v = (u-u_1)/(dt_const*omega_const) - (1-omega_const)/omega_const*v_1

    # Compute residual
    R = (inner(Sc, ddotE) + dot(J*(K_perm*invF.T*grad(p)), invF.T*grad(q)))*dx \
        - p*J*inner(ddotF, invF.T)*dx                                         \
        + (q*J*inner(grad(v),invF.T))*dx                                      \
        + (inner(B,w))*dx - (inner(Trac,w))*ds_neumann(1)                     \
        + (inner(g_bar,q))*ds_neumann(0)

    # Compute Jacobian of R
    Jac = derivative(R, up, dup)

    # Set up the problem
    problem = NonlinearVariationalProblem(R, up, bcs=bcs, J=Jac)
    solver = NonlinearVariationalSolver(problem)

    solver.parameters["nonlinear_solver"] = "newton"
    solver.parameters["newton_solver"]["linear_solver"] = "bicgstab"
    solver.parameters["newton_solver"]["preconditioner"] = "ilu"

    # solver.parameters["nonlinear_solver"] = "snes"
    # solver.parameters["snes_solver"]["line_search"] = "bt"
    # solver.parameters["snes_solver"]["linear_solver"] = "gmres"
    # solver.parameters["snes_solver"]["preconditioner"] = "ilu"        
    # solver.parameters["snes_solver"]["method"] = "tr"


    ## Save initial conditions in VTK format
    # assign(up, up_1)
    dfile = File("results/displacement.pvd");
    pfile = File("results/pressure.pvd");
    dfile << (up.sub(0),0.0);
    pfile << (up.sub(1),0.0);

    ### Run Simulation
    u_max = [0.0]
    while t<T_total:
        print 'time = ', t    

        # newton_solver(up, dup, up_1, R, Jac, bcs, 1e-10, 100, V)

        # solve
        solver.solve()

        # update
        t    += dt
        tn   += 1
        v_1  = project(v,Pu)
        assign(up_1.sub(0),up.sub(0))
        H_1 = project(H, HS)
        S_1 = project(S, HS)

        # Save solution in VTK format
        dfile << (up.sub(0), t);
        pfile << (up.sub(1), t);

        # plot(p_1, title = "pressure", axes=True, interactive = True)

        # pdb.set_trace()
        u_tent, p_tent = up.split(deepcopy=True) 
        # u_max.append( np.max(u_tent.vector().array()))
        u_max.append( np.max(p_tent.vector().array()))


    return u_max

    # Plot and hold solution
    # plot(u, mode = "displacement", title="displacement", axes=True, interactive = True)
    # plot(p, title  = "pressure", axes=True, interactive = True)
