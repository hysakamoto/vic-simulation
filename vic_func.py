""" Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

# Begin simulation

from dolfin import *
import numpy as np
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


def vic_sim( m_num, p_order, dt, T_total, omega ):
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
    Trac  = Constant((0.0,  0.0, 5.0)) # Traction force on the boundary
    g_bar = Constant(0.0)            # Normal flux

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

    # Mark Dirichlet boundaries
    top    = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)
    bottom = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)

    # Assign Dirichlet boundaries (x = 0 or x = 1)
    d_top     = Expression(("0.0", "0.0", "0.0"))
    d_bottom  = Expression(("0.0", "0.0", "0.0"))
    bc_top    = DirichletBC(V.sub(0), d_top, top)
    bc_bottom = DirichletBC(V.sub(0), d_bottom, bottom)

    p_top = Expression("0.0")
    p_bottom   = Expression("1.0")
    bc_ptop    = DirichletBC(V.sub(1), p_top, top)
    bc_pbottom = DirichletBC(V.sub(1), p_bottom, bottom)

    bcs = [bc_bottom, bc_ptop]

    ## Initial conditions
    u_1 = Function(Pu)
    u_1.interpolate(Constant((0.0, 0.0, 0.0)))
    p_1 = interpolate(Constant(0.0), Pp)

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
    Ee, nu    = 10.0, 0.45
    mu, lmbda = Constant(Ee/(2*(1 + nu))), Constant(Ee*nu/((1 + nu)*(1 - 2*nu)))

    # Permeability
    perm   = Constant(1.0)
    K_perm = perm*I

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
    psi_1   = mu/2.0*(tr(C_hat_1)-3) + 1.0/2.0*kappa*ln(J_1)**2.0
    # PK1 stress tensor
    P_1 = diff(psi_1,F_1)
    # PK2 stress tensor
    S_1 = invF_1*P_1

    H = exp(-dt_const/tau)*S_1 + (1-exp(-dt_const/tau))*(S-S_1)/(dt_const/tau)
    Sc = S+gamma*H


    ## Poroelasticity!!!!
    v_1 = Function(Pu)
    v_1.interpolate(Constant((0.0, 0.0, 0.0)))
    v = (u-u_1)/(dt_const*omega) - (1-omega)/omega*v_1

    # Compute residual
    R = (inner(Sc, ddotE) + dot(J*(K_perm*invF.T*grad(p)), invF.T*grad(q)))*dx \
        - p*J*inner(ddotF, invF.T)*dx                                         \
        + (q*J*inner(grad(v),invF.T))*dx                                      \
        + (inner(B,w))*dx - (inner(Trac,w))*ds_neumann(0)                     \
        + (inner(g_bar,q))*ds_neumann(1)

    # Compute Jacobian of R
    Jac = derivative(R, up, dup)

    # Set up the problem
    problem = NonlinearVariationalProblem(R, up, bcs=bcs, J=Jac )
    solver = NonlinearVariationalSolver(problem)
    solver.parameters["newton_solver"]["linear_solver"] = "bicgstab"
    solver.parameters["newton_solver"]["preconditioner"] = "ilu"

    ## Save initial conditions in VTK format
    dfile = File("results/displacement.pvd");
    pfile = File("results/pressure.pvd");
    dfile << (up.sub(0),0.0);
    pfile << (up.sub(1),0.0);

    ### Run Simulation
    u_max = []
    while t<T_total:
        print 'time = ', t    

        # solve
        solver.solve()

        # update
        t    += dt
        tn   += 1
        v_1  = project(v,Pu)
        assign(u_1,up.sub(0))

        # Save solution in VTK format
        dfile << (up.sub(0), t);
        pfile << (up.sub(1), t);

        # plot(p_1, title = "pressure", axes=True, interactive = True)

        # pdb.set_trace()
        u_max.append( np.max(u_1.vector().array()))


    return u_max

    # Plot and hold solution
    # plot(u, mode = "displacement", title="displacement", axes=True, interactive = True)
    # plot(p, title  = "pressure", axes=True, interactive = True)
