## generation of the body force and boundary traction from the 
## exact solutions.

from dolfin import *

# mesh = UnitCubeMesh(m_num, m_num, m_num)
# dim  = mesh.topology().dim() 

# ##  Function Spaces
# Pu = VectorFunctionSpace(mesh, "Lagrange", p_order)  # space for displacements
# Pp = FunctionSpace(mesh, "Lagrange", p_order)        # space for pressure
# V  = MixedFunctionSpace([Pu,Pp])                    # mixed space


def gen_conditions( u_e_str, p_e_str, v_e_str, Pu, Pp ):

    ## exact solution on function space
    u_e = Expression(u_e_str, t=dt, element=Pu.ufl_element())

    p_e = Expression(p_e_str, t=dt, element=Pp.ufl_element())

    v_e = Expression(v_e_str, t=dt, element=Pu.ufl_element())


    ## Kinematics
    I    = Identity(dim)           # Identity tensor
    F    = I + grad(u_e)             # Deformation gradient
    F    = variable(F)             # Make F a variable for tensor differentiations
    C    = F.T*F                   # Right Cauchy-Green tensor
    invF = inv(F)

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

    # source term
    src = J*inner(grad(v_e),invF.T) - div( J*invF*K_perm*invF.T*grad(p_e))

    # body force
    bf = div(F*S-p_e*J*invF.T)

    flx = -(invF*J*K_perm*invF.T*grad(p_e))
    SE = (-p_e*J*invF.T + F*S)

    # unit normal vectors
    n0 = as_vector([1.0, 0.0, 0.0])

    # boundary conditions
    g0_bar = inner(flx,n0)
    t0_bar = SE*n0

    return src, bf, g0_bar, t0_bar
