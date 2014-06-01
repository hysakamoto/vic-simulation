## Kinematics

from dolfin import *

def identity_tensor(dim):
    return Identity(dim)

def deformation_gradient(u):
    return variable (identity_tensor(3) + grad(u))

def right_CG_tensor(u):
    F = deformation_gradient(u)
    return variable(F.T*F)

def inverse_deformation_gradient(u):
    F = deformation_gradient(u)
    return variable(inv(F))

def volume_change(u):
    F = deformation_gradient(u)
    return variable(det(F))

def invariants(A):
    I1 = tr(A)
    I2 = 0.5 * (tr(A)**2 - tr(A*A))
    I3 = det(A)    
    return variable(I1), variable(I2), variable(I3)

def isochronic_right_CG_tensor(u):
    C = right_CG_tensor(u)
    Ic1, Ic2, Ic3 = invariants(C)
    return variable(Ic3**(-1.0/3.0)*C)

def neo_hookean_incompressible_energy(u, mu, kappa):
    C_hat = isochronic_right_CG_tensor(u)
    J = volume_change(u)
    return variable(mu/2.0*(tr(C_hat)-3) + 1.0/2.0*kappa*ln(J)**2.0)

def neo_hookean_incompressible_PK2_stress(u, mu, kappa):
    psi = neo_hookean_incompressible_energy(u, mu, kappa)
    F = deformation_gradient(u)
    return variable(inv(F)*diff(psi,F))

def neo_hookean_compressible_energy(u, mu, lmbda):
    C = right_CG_tensor(u)
    Ic1, Ic2, Ic3 = invariants(C)
    J = volume_change(u)
    # Strain energy density (compressible neo-Hookean model)
    return (mu/2.0)*(Ic1 - 3) - mu*ln(J) + (lmbda/2.0)*(ln(J))**2
