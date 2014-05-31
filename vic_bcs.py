### Define boundary conditions
from dolfin import *
import numpy as np

def neumann_boundaries(tol, exterior_facet_domains):
    '''Mark Neumann Boundaries'''
    class neum_top(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[2] - 1.0) < tol
    class neum_bottom(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[2]) < tol
    class neum_right(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[0] - 1.0) < tol
    class neum_left(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[0]) < tol
    class neum_back(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[1] - 1.0) < tol
    class neum_front(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[1]) < tol

    # first mark everything with 0
    exterior_facet_domains.set_all(0)

    gamma_top = neum_top()
    gamma_bottom = neum_bottom()
    gamma_right = neum_right()
    gamma_left = neum_left()
    gamma_back = neum_back()
    gamma_front = neum_front()

    # mark boundaries with 0~5
    gamma_top.mark(exterior_facet_domains, 0)
    gamma_bottom.mark(exterior_facet_domains, 1)
    gamma_right.mark(exterior_facet_domains, 2)
    gamma_left.mark(exterior_facet_domains, 3)
    gamma_back.mark(exterior_facet_domains, 4)
    gamma_front.mark(exterior_facet_domains, 5)

    ds_neumann = ds[exterior_facet_domains]

    return ds_neumann

def neumann_expressions(t, mu):
    '''Return the tuple of the neumann boundary conditions'''

    ## normal flux neumann bc
    gbar_top= Expression("(((400.0*((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)))))*(x[2]-(x[2]*((pow((t-20.0),2.0)/400.0)-1.0))))*pow((1.0/pow((((-pow(t,2.0))+(40.0*t))+400.0),2.0)),(1.0/2.0)))", t=t)
    gbar_bottom= Expression("(((-400.0*((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)))))*(x[2]-(x[2]*((pow((t-20.0),2.0)/400.0)-1.0))))*pow((1.0/pow((((-pow(t,2.0))+(40.0*t))+400.0),2.0)),(1.0/2.0)))", t=t)
    gbar_right= Expression("0.0", t=t)
    gbar_left= Expression("0.0", t=t)
    gbar_back= Expression("0.0", t=t)
    gbar_front= Expression("0.0", t=t)

    ## traction neumann bc
    tbar_top= Expression(("0.0",
                          "0.0",
                          "((400.0*(((((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))))*pow((x[2]-(x[2]*((pow((t-20.0),2.0)/400.0)-1.0))),2.0))/2.0)+(mu*(pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)-1.0))))*pow((1.0/pow((((-pow(t,2.0))+(40.0*t))+400.0),2.0)),(1.0/2.0)))"), t=t, mu=mu)
    tbar_bottom= Expression(("0.0",
                             "0.0",
                             "((-400.0*(((((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))))*pow((x[2]-(x[2]*((pow((t-20.0),2.0)/400.0)-1.0))),2.0))/2.0)+(mu*(pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)-1.0))))*pow((1.0/pow((((-pow(t,2.0))+(40.0*t))+400.0),2.0)),(1.0/2.0)))"), t=t, mu=mu)
    tbar_right= Expression(("((((((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))))*pow((x[2]-(x[2]*((pow((t-20.0),2.0)/400.0)-1.0))),2.0))/2.0)-(mu*((1.0/((pow((t-20.0),2.0)/400.0)-2.0))+1.0)))*pow((((400.0/(((-pow(t,2.0))+(40.0*t))+400.0))*pow((((-pow(t,2.0))+(40.0*t))+400.0),2.0))/160000.0),(1.0/2.0)))",
                            "0.0",
                            "0.0"), t=t, mu=mu)
    tbar_left= Expression(("((-(((((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))))*pow((x[2]-(x[2]*((pow((t-20.0),2.0)/400.0)-1.0))),2.0))/2.0)-(mu*((1.0/((pow((t-20.0),2.0)/400.0)-2.0))+1.0))))*pow((((400.0/(((-pow(t,2.0))+(40.0*t))+400.0))*pow((((-pow(t,2.0))+(40.0*t))+400.0),2.0))/160000.0),(1.0/2.0)))",
                           "0.0",
                           "0.0"), t=t, mu=mu)
    tbar_back= Expression(("0.0",
                           "((((((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))))*pow((x[2]-(x[2]*((pow((t-20.0),2.0)/400.0)-1.0))),2.0))/2.0)-(mu*((1.0/((pow((t-20.0),2.0)/400.0)-2.0))+1.0)))*pow((((400.0/(((-pow(t,2.0))+(40.0*t))+400.0))*pow((((-pow(t,2.0))+(40.0*t))+400.0),2.0))/160000.0),(1.0/2.0)))",
                           "0.0"), t=t, mu=mu)
    tbar_front= Expression(("0.0",
                            "((-(((((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))))*pow((x[2]-(x[2]*((pow((t-20.0),2.0)/400.0)-1.0))),2.0))/2.0)-(mu*((1.0/((pow((t-20.0),2.0)/400.0)-2.0))+1.0))))*pow((((400.0/(((-pow(t,2.0))+(40.0*t))+400.0))*pow((((-pow(t,2.0))+(40.0*t))+400.0),2.0))/160000.0),(1.0/2.0)))",
                            "0.0"), t=t, mu=mu)

    return gbar_top, gbar_bottom, gbar_right, gbar_left, gbar_back, gbar_front, \
        tbar_top, tbar_bottom, tbar_right, tbar_left, tbar_back, tbar_front


def dirichlet_boundaries(tol, V, t, u_e, p_e):
    '''Define Dirichlet Boundaries'''

    # Define Dirichlet boundaries
    def diri_top(x, on_boundary):
        return on_boundary and (abs(x[2]-1.0)<tol)

    def diri_bottom(x, on_boundary):
        return on_boundary and (abs(x[2])<tol)

    def diri_right(x, on_boundary):
        return on_boundary and (abs(x[0]-1.0)<tol)

    def diri_left(x, on_boundary):
        return on_boundary and (abs(x[0])<tol)

    def diri_back(x, on_boundary):
        return on_boundary and (abs(x[1]-1.0)<tol)

    def diri_front(x, on_boundary):
        return on_boundary and (abs(x[1])<tol)

    # Assign Dirichlet boundaries - displacements
    bc_utop    = DirichletBC(V.sub(0), u_e, diri_top)
    bc_ubottom = DirichletBC(V.sub(0), u_e, diri_bottom)
    bc_uright  = DirichletBC(V.sub(0), u_e, diri_right)
    bc_uleft   = DirichletBC(V.sub(0), u_e, diri_left)
    bc_uback   = DirichletBC(V.sub(0), u_e, diri_back)
    bc_ufront  = DirichletBC(V.sub(0), u_e, diri_front)

    # Assign Dirichlet boundaries - pressure
    bc_ptop    = DirichletBC(V.sub(1), p_e, diri_top)
    bc_pbottom = DirichletBC(V.sub(1), p_e, diri_bottom)
    bc_pright = DirichletBC(V.sub(1), p_e, diri_right)
    bc_pleft = DirichletBC(V.sub(1), p_e, diri_left)
    bc_pback = DirichletBC(V.sub(1), p_e, diri_back)
    bc_pfront = DirichletBC(V.sub(1), p_e, diri_front)
    
    return bc_utop, bc_ubottom, \
        bc_uright, bc_uleft, \
        bc_uback, bc_ufront, \
        bc_ptop, bc_pbottom, bc_pright, bc_pleft, bc_pback, bc_pfront
