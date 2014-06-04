""" Simulation of VIC"""

# Simulation of poro-vosco-hyperelastic solid in initial configuration

import pdb
import numpy as np
import scipy.io

from dolfin import *

import newton_solve
from vic_bcs import *
from manufactured_solutions import getManuSolutions
from cui_toolbar import toolbar

# DEBUG
# set_log_level(DEBUG) #PROGRESS
# set_log_level(PROGRESS) #

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
# parameters["form_compiler"]["quadrature_degree"] = 2
# parameters["num_threads"] = 1
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}


def vic_sim( sim_name, \
             m_num, p_order, dt, T_total, max_it, \
             omega, Ee, nu, gamma, tau, perm, \
             top_trac, body_force \
             ):
    """ 
    run the VIC simulation
    """

    ## avoid recompilation
    dt_const = Constant(dt)
    
    # Elasticity parameters
    mu, lmbda = Constant(Ee/(2*(1 + nu))), Constant(Ee*nu/((1 + nu)*(1 - 2*nu)))

    # Exact solutions
    [u_e, p_e, v_e, source, body_force, tbars, gbars, u_initial, p_initial, v_initial] \
        = getManuSolutions(dt, perm, mu, lmbda)
    
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
    up   = Function(V)         # Velocity-pressure from previous iteration
    v, p = split(up)           # Function in each subspace to write the functional
    w, q = split(wq)           # Test Function split

    ## Boundary Conditions

    # Create mesh function over cell facets
    exterior_facet_domains = MeshFunction("size_t", mesh, mesh.topology().dim()-1)

    # Define Neumann boundaries
    bd_tol = 1E-14   # tolerance for coordinate comparisons
    ds_neumann = neumann_boundaries(bd_tol, exterior_facet_domains)

    # Get the neumann boundary conditions
    gbar_top, gbar_bottom, gbar_right, gbar_left, gbar_back, gbar_front \
        = gbars
    tbar_top, tbar_bottom, tbar_right, tbar_left, tbar_back, tbar_front \
        = tbars

    # Define Dirichlet boundaries
    bc_utop, bc_ubottom, bc_uright, bc_uleft, bc_uback, bc_ufront, \
        bc_ptop, bc_pbottom, bc_pright, bc_pleft, bc_pback, bc_pfront \
        = dirichlet_boundaries(bd_tol, V, dt, v_e, p_e)

    # bcs = [bc_utop, bc_ubottom, bc_uright, bc_uleft, bc_uback, bc_ufront, \
    #     bc_ptop, bc_pbottom, bc_pright, bc_pleft, bc_pback, bc_pfront]

    bcs = [bc_ubottom,\
           bc_pbottom, bc_pleft, bc_pfront]

    ## Initial conditions
    u_1 = Function(Pu)
    assign (u_1, interpolate(u_initial,Pu))

    ## Poroelasticity!!!!
    omega_const = Constant(omega)
    v_1 = Function(Pu)
    assign (v_1, interpolate(v_initial,Pu))
    u = (omega_const*v + (1.0-omega_const)*v_1)*dt_const + u_1

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
    # kappa = mu*10
    # C_hat = (IIIc)**(-1.0/3.0)*C
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
    # HS = TensorFunctionSpace(mesh, "Lagrange", p_order )
    # H_1 = Function(HS)
    # H_1.interpolate(Constant(np.zeros([3,3])))
    
    ## viscoelastic parameters
    # gamma_const = Constant(gamma)
    # tau_const = Constant(tau)

    # ## Kinematics
    # F_1    = I + grad(u_1)             # Deformation gradient
    # F_1    = variable(F_1)             # Make F a variable for tensor differentiations
    # C_1    = F_1.T*F_1                   # Right Cauchy-Green tensor
    # invF_1 = inv(F_1)
    # # Invariants of deformation tensors
    # J_1    = det(F_1)
    # Ic_1   = tr(C_1)
    # IIIc_1 = det(C_1)
    # # Strain energy density (nearly incompressible neo-hookean model)
    # C_hat_1 = (IIIc_1)**(-1.0/3.0)*C_1
    # # psi_1   = mu/2.0*(tr(C_hat_1)-3) + 1.0/2.0*kappa*ln(J_1)**2.0
    # # compressible neo-Hookean
    # psi_1 = (mu/2.0)*(Ic_1 - 3.0) - mu*ln(J_1) + (lmbda/2.0)*(ln(J_1))**2.0
    # # PK1 stress tensor
    # P_1 = diff(psi_1,F_1)
    # # PK2 stress tensor
    # S_1 = invF_1*P_1

    # H = exp(-dt_const/tau_const)*H_1 + (1-exp(-dt_const/tau_const))*(S-S_1)/(dt_const/tau_const)
    # Sc = S+gamma_const*H
    Sc = S

    # Compute residual
    R = (inner(Sc, ddotE) + dot(J*(K_perm*invF.T*grad(p)), invF.T*grad(q)))*dx \
        - p*J*inner(ddotF, invF.T)*dx                                         \
        + (q*J*inner(grad(v),invF.T))*dx                                      \
        - (source*q)*dx \
        + (dot(body_force,w))*dx \
        - (dot(tbar_top,w))*ds_neumann(0) \
        - (dot(tbar_right,w))*ds_neumann(2) \
        - (dot(tbar_back,w))*ds_neumann(4) \
        + (gbar_top*q)*ds_neumann(0)\
        + (gbar_right*q)*ds_neumann(2) \
        + (gbar_back*q)*ds_neumann(4) \
        - (dot(tbar_left,w))*ds_neumann(3) \
        - (dot(tbar_front,w))*ds_neumann(5) \

    # Compute Jacobian of R
    Jac = derivative(R, up, dup)

    # Set up the problem
    problem = NonlinearVariationalProblem(R, up, bcs=bcs, J=Jac)
    solver = NonlinearVariationalSolver(problem)

    # solver.parameters["nonlinear_solver"] = "newton"
    # solver.parameters["newton_solver"]["linear_solver"] = "bicgstab"
    # solver.parameters["newton_solver"]["preconditioner"] = "sor"

    solver.parameters["nonlinear_solver"] = "snes"
    solver.parameters["snes_solver"]["line_search"] = "bt"
    solver.parameters["snes_solver"]["linear_solver"] = "gmres"
    solver.parameters["snes_solver"]["preconditioner"] = "sor" # ilu
    solver.parameters["snes_solver"]["method"] = "tr"

    ## Save initial conditions in VTK format
    assign (up.sub(0), interpolate(v_initial,Pu))
    assign (up.sub(1), interpolate(p_initial,Pp))

    dfile = File(sim_name + "/displacement.pvd");
    vfile = File(sim_name + "/velocity.pvd");
    pfile = File(sim_name + "/pressure.pvd");
    dpfile = File(sim_name + "/up.pvd");
    u_1.rename('u', 'u_solution')
    dfile << (u_1,0.0);
    vfile << (up.sub(0), 0.0)
    pfile << (up.sub(1),0.0);
    # Save solutions in xml format
    File(sim_name+ '/up_%d.xml' %0) << up
    File(sim_name+ '/u_%d.xml' %0) << u_1;

    # Save mesh
    # File(sim_name+'/mesh.xdmf') << mesh

    ### Run Simulation
    while tn<max_it:
        print 'time = ', (t+dt)

        ## solve
        solver.solve()
        
        ## Update the _1 values
        # Calculate the new value of v_1 point-wise
        # u_tent = up.sub(0,deepcopy=True).vector().get_local()
        # u_1_tent = up_1.sub(0,deepcopy=True).vector().get_local()
        # v_1_tent = v_1.vector().get_local()
        # v_1.vector().set_local(((u_tent-u_1_tent)/dt - (1.0-omega)*v_1_tent)/omega)

        # Project the new value of v_1
        # u_1_local = u_1.vector().get_local()
        # v_local = up.sub(0,deepcopy=True).vector().get_local()
        # v_1_local = v_1.vector().get_local()
        # u_1.vector().set_local( u_1_local + (omega*v_local + (1.0-omega)*v_1_local)*dt )

        u_1.vector()[:] =  u_1.vector() + (omega*up.sub(0,deepcopy=True).vector() + (1.0-omega)*v_1.vector()) *dt
        # u_1.vector().axpy(1.0, (omega*up.sub(0,deepcopy=True).vector() + (1.0-omega)*v_1.vector()) *dt)

        # u_1 = project(u,Pu)

        assign(v_1, up.sub(0))

        ## update time
        t    += dt
        tn   += 1

        ## Update parameters
        body_force.t = t+dt
        source.t = t+dt
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
        # update exact solutions
        u_e.t = t+dt
        v_e.t = t+dt
        p_e.t = t+dt

        # Define Dirichlet boundaries
        bc_utop, bc_ubottom, bc_uright, bc_uleft, bc_uback, bc_ufront, \
            bc_ptop, bc_pbottom, bc_pright, bc_pleft, bc_pback, bc_pfront \
            = dirichlet_boundaries(bd_tol, V, t+dt, v_e, p_e)

        bcs = [bc_ubottom,\
               bc_pbottom, bc_pleft, bc_pfront]
        # bcs = [bc_utop, bc_ubottom, bc_uright, bc_uleft, bc_uback, bc_ufront, \
        #        bc_ptop, bc_pbottom, bc_pright, bc_pleft, bc_pback, bc_pfront]

        # define problem
        problem = NonlinearVariationalProblem(R, up, bcs=bcs, J=Jac)
        solver = NonlinearVariationalSolver(problem)

        # Save solution in VTK format
        u_1.rename('u', 'u_solution') 
        dfile << (u_1, t);
        vfile << (up.sub(0), t);
        pfile << (up.sub(1), t);
        # Save solutions in xml format
        File(sim_name+ '/up_%d.xml' %tn) << up
        File(sim_name+ '/u_%d.xml' %tn) << u_1

        # plot(p_1, title = "pressure", axes=True, interactive = True)

    return 0


def save_mat(up, sim_name):
    # Save solutions in MATLAB format using scipy.io
    uvec = [up.sub(0).sub(0,deepcopy=True).vector().array(), 
            up.sub(0).sub(1,deepcopy=True).vector().array(), 
            up.sub(0).sub(2,deepcopy=True).vector().array()]
    pvec = up.sub(1,deepcopy=True).vector().array()
    scipy.io.savemat(sim_name + '/u_%d' %tn, { 'u': uvec })
    scipy.io.savemat(sim_name + '/p_%d' %tn, { 'p': pvec })
