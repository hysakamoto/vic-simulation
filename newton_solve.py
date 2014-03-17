from dolfin import *
import numpy
import pdb
import pdb

def newton_solver( u, du, u_1, R, Jac, bcs, tol, maxiter, V ):

    # initial solve
    # A,b = assemble_system(Jac, R, bcs)
    # solve(A, u_1.vector(), b, 'lu')
    # u.vector()[:] = u_1.vector() - rlx*du.vector()
    # u_1.assign(u)
    

    # Mark Dirichlet boundaries
    top    = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)
    bottom = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)

    ## ZERO the Dirichlet boundaries
    d_top     = Expression(("0.0", "0.0", "0.0"))
    d_bottom  = Expression(("0.0", "0.0", "0.0"))
    bc_top    = DirichletBC(V.sub(0), d_top, top)
    bc_bottom = DirichletBC(V.sub(0), d_bottom, bottom)
    p_top = Expression("0.0")
    p_bottom   = Expression("0.0")
    bc_ptop    = DirichletBC(V.sub(1), p_top, top)
    bc_pbottom = DirichletBC(V.sub(1), p_bottom, bottom)
    bcs = [bc_bottom, bc_ptop, bc_pbottom]

    iter = 0
    eps = tol*10
    rlx = 1.0
    while eps > tol and iter < maxiter:
        # pdb.set_trace()
        iter += 1
        A, b = assemble_system(Jac, -R, bcs)
        solve(A, du.vector(), b)
        eps = numpy.linalg.norm(du.vector().array())
        print 'Norm:', eps
        u.vector()[:] += rlx*du.vector()
