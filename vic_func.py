""" Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

# Begin simulation
import pdb

from dolfin import *
import numpy as np
import newton_solve
from vic_bcs import *

# DEBUG
# set_log_level(DEBUG) #PROGRESS

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 2
parameters["num_threads"] = 2
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}


def vic_sim( m_num, p_order, dt, T_total, max_it, \
             omega, Ee, nu, gamma, tau, perm, \
             top_trac, body_force ):
    """ 
    run the VIC simulation
    """


    # Exact solutions
    u_e = Expression(("((x[0]*pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(1.0/2.0)))-x[0])","((x[1]*pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(1.0/2.0)))-x[1])","((-x[2])*((pow((t-20.0),2.0)/400.0)-1.0))"), t=dt)
    p_e = Expression("((-(((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))))*pow((x[2]-(x[2]*((pow((t-20.0),2.0)/400.0)-1.0))),2.0)))/2.0)", t=dt)


    ## avoid recompilation
    dt_const = Constant(dt)
    
    # Elasticity parameters
    mu, lmbda = Constant(Ee/(2*(1 + nu))), Constant(Ee*nu/((1 + nu)*(1 - 2*nu)))
    
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

    ## Loads: Body Force
    B = Expression(("0.0", "0.0", 
                   "(((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))))*(x[2]-(x[2]*((pow((t-20.0),2.0)/400.0)-1.0))))"), 
                   t=dt)

    ## Boundary Conditions

    # Create mesh function over cell facets
    exterior_facet_domains = MeshFunction("size_t", mesh, mesh.topology().dim()-1)

    # Define Neumann boundaries
    bd_tol = 1E-14   # tolerance for coordinate comparisons
    ds_neumann = neumann_boundaries(bd_tol, exterior_facet_domains)

    # Get the neumann boundary conditions
    gbar_top, gbar_bottom, gbar_right, gbar_left, gbar_back, gbar_front, \
        tbar_top, tbar_bottom, tbar_right, tbar_left, tbar_back, tbar_front \
        = neumann_expressions(dt, mu)

    # Define Dirichlet boundaries
    bc_utop, bc_ubottom, bc_uright, bc_uleft, bc_uback, bc_ufront, \
        bc_ptop, bc_pbottom, bc_pright, bc_pleft, bc_pback, bc_pfront \
        = dirichlet_boundaries(bd_tol, V, dt)

    # bcs = [bc_ubottom, bc_uleft, bc_ufront, bc_utop, bc_uright, bc_uback, 
           # bc_pbottom, bc_pleft, bc_pright]
    bcs = [bc_ubottom, bc_uleft, bc_ufront, bc_pbottom, bc_pleft, bc_pfront]
    # bcs = [bc_ubottom, 
           # bc_pbottom, bc_ptop, bc_pleft, bc_pright, bc_pfront, bc_pback]

    ## Initial conditions
    up_1   = Function(V)         # Displacement-pressure from previous iteration
    u_1, p_1 = split(up_1)       # Function in each subspace to write the functional
    assign (up_1.sub(0), interpolate(Constant((0.0, 0.0, 0.0)),Pu))
    assign (up_1.sub(1), interpolate(Constant(0),Pp))

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
    psi = (mu/2.0)*(Ic - 3) - mu*ln(J) + (lmbda/2.0)*(ln(J))**2

    # Strain energy density (nearly incompressible neo-hookean model)
    kappa = mu*10
    C_hat = (IIIc)**(-1.0/3.0)*C
    # psi   = mu/2.0*(tr(C_hat)-3) + 1.0/2.0*kappa*ln(J)**2.0

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

    H = exp(-dt_const/tau_const)*H_1 + (1-exp(-dt_const/tau_const))*(S-S_1)/(dt_const/tau_const)
    Sc = S+gamma_const*H

    ## Poroelasticity!!!!
    omega_const = Constant(omega)
    v_1 = Function(Pu)
    v_1.interpolate(Constant((0.0, 0.0, 0.0)))
    v = ((u-u_1)/dt_const - (1-omega_const)*v_1)/omega_const

    # Compute residual
    R = (inner(Sc, ddotE) + dot(J*(K_perm*invF.T*grad(p)), invF.T*grad(q)))*dx \
        - p*J*inner(ddotF, invF.T)*dx                                         \
        + (q*J*inner(grad(v),invF.T))*dx                                      \
        + (inner(B,w))*dx \
        - (inner(tbar_top,w))*ds_neumann(0) \
        - (inner(tbar_right,w))*ds_neumann(2) \
        - (inner(tbar_back,w))*ds_neumann(4) \
        + (inner(gbar_top,q))*ds_neumann(0)\
        + (inner(gbar_right,q))*ds_neumann(2) \
        + (inner(gbar_back,q))*ds_neumann(4) \


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
    assign(up, up_1)
    dfile = File("results/displacement.pvd");
    pfile = File("results/pressure.pvd");
    dfile << (up.sub(0),0.0);
    pfile << (up.sub(1),0.0);

    # dup   = Function(V)

    ### Run Simulation
    u_max = []
    Eus = []
    Eps = []
    while tn<max_it:
        print 'time = ', (t+dt)

        # solve
        solver.solve()

        # update
        t    += dt
        tn   += 1

        # update body force
        B.t = t+dt
        # update boundary conditions
        gbar_top.t = t+dt
        gbar_bottom.t = t+dt
        gbar_right.t = t+dt
        gbar_left.t = t+dt
        gbar_back.t = t+dt
        gbar_front.t = t+dt
        tbar_top.t = t+dt
        tbar_bottom.t = t+dt
        tbar_right.t = t+dt
        tbar_left.t = t+dt
        tbar_back.t = t+dt
        tbar_front.t = t+dt 

        # Define Dirichlet boundaries
        bc_utop, bc_ubottom, bc_uright, bc_uleft, bc_uback, bc_ufront, \
            bc_ptop, bc_pbottom, bc_pright, bc_pleft, bc_pback, bc_pfront \
            = dirichlet_boundaries(bd_tol, V, t+dt)
        # bcs = [bc_ubottom,
               # bc_pbottom, bc_ptop, bc_pleft, bc_pright, bc_pfront, bc_pback]
        bcs = [bc_ubottom, bc_uleft, bc_ufront, bc_pbottom, bc_pleft, bc_pfront]

        # define problem
        problem = NonlinearVariationalProblem(R, up, bcs=bcs, J=Jac)
        solver = NonlinearVariationalSolver(problem)

        v_1  = project(v,Pu)
        up_1.assign(up)
        assign(up_1.sub(0),up.sub(0))
        H_1 = project(H, HS)
        S_1 = project(S, HS)

        # Save solution in VTK format
        dfile << (up.sub(0), t);
        pfile << (up.sub(1), t);

        # plot(p_1, title = "pressure", axes=True, interactive = True)

        # pdb.set_trace()
        u_tent, p_tent = up.split(deepcopy=True) 

        
        ### Error against exact solutions
        error_u = (u_tent-u_e)**2*dx
        Eu = sqrt(assemble(error_u))

        error_p = (p_tent-p_e)**2*dx
        Ep = sqrt(assemble(error_p))

        # Explicit interpolation of u_e onto the same space as u:
        # u_e_V = interpolate(u_e, V)
        # error_u = (u - u_e_V)**2*dx
        # E2 = sqrt(assemble(error_u))

        Eus.append(Eu)
        Eps.append(Ep)
        u_max.append( np.max(u_tent.vector().array()))

        # update body force
        u_e.t = t+dt
        p_e.t = t+dt

    return u_max, Eus, Eps

    # Plot and hold solution
    # plot(u, mode = "displacement", title="displacement", axes=True, interactive = True)
    # plot(p, title  = "pressure", axes=True, interactive = True)
