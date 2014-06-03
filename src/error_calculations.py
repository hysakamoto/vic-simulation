## Calculate L2 errors

import os
import pdb
from time import time as ttime

from dolfin import *
from manufactured_solutions import getManuSolutions
from cui_toolbar import toolbar

def errorCalc(base_name, max_mer, max_mit, sim_params, mat_params):

    sim_name = base_name + str(sim_params['m_num']) \
               + '-' + str(sim_params['max_it']) + '/'

    if not os.path.exists(sim_name):
        print 'The directory does not exist!'
        return 

    # Get manufactured solutions
    mu = mat_params['Ee']/(2*(1 + mat_params['nu']))
    lmbda = mat_params['Ee']*mat_params['nu'] \
            /((1 + mat_params['nu'])*(1 - 2*mat_params['nu']))
    [u_e, p_e, v_e, source, body_force, tbars, gbars, u_initial, p_initial, v_initial] \
        = getManuSolutions(0.0, mat_params['perm'], mu, lmbda)

    Eu = 0.0
    Ep = 0.0
    Eu2 = 0.0
    Ep2 = 0.0
    Eu3 = 0.0
    Ep3 = 0.0
    Eu4 = 0.0
    Ep4 = 0.0

    ##  Finest Function Spaces
    mesh_e = UnitCubeMesh(max_mer,max_mer,max_mer)
    Pu_e = VectorFunctionSpace(mesh_e, "Lagrange", sim_params['p_order'])
    Pp_e = FunctionSpace(mesh_e, "Lagrange", sim_params['p_order'])
    V_e  = MixedFunctionSpace([Pu_e,Pp_e])
    
    # number of error computation per time step
    n_err_comp = max_mit/sim_params['max_it'] 
    dt = sim_params['T_total']/float(sim_params['max_it'] )

    ##  Mesh and Function Spaces
    # mesh = Mesh(sim_name+'/mesh.xdmf')
    mesh = UnitCubeMesh(sim_params['m_num'], sim_params['m_num'], sim_params['m_num'])
    File(sim_name+'/mesh.xdmf') << mesh
    Pu = VectorFunctionSpace(mesh, "Lagrange", sim_params['p_order'])
    Pp = FunctionSpace(mesh, "Lagrange", sim_params['p_order'])
    V  = MixedFunctionSpace([Pu,Pp])

    t = 0.0
    tb = toolbar(30)
    for tn in range(1, sim_params['max_it']+1):

        # load solutions
        up_1 = Function(V, sim_name + '/up_%d.xml'%(tn-1))
        u_1 = up_1.sub(0)
        p_1 = up_1.sub(1)

        up_2 = Function(V, sim_name + '/up_%d.xml'%(tn))
        u_2 = up_2.sub(0)
        p_2 = up_2.sub(1)

        t1 = Constant(t)
        t2 = Constant(t+dt)

        t_1 = t
        t_2 = t+dt

        ddt = dt/(n_err_comp)

        # error computations within a timestep
        for j in range(n_err_comp):
            t += ddt
            tb.update(t/sim_params['T_total'])

            # set time (use a half-point for a better integration)
            u_e.t = t-ddt/2.0
            p_e.t = t-ddt/2.0
            t_const = Constant(t-ddt/2.0)

            uh = (u_2-u_1)/(t2-t1)*(t_const-t1) + u_1
            ph = (p_2-p_1)/(t2-t1)*(t_const-t1) + p_1

            ### L2 norms
            error_u = (uh-u_e)**2*dx
            error_p = (ph-p_e)**2*dx
            Eu += (assemble(error_u))*ddt
            Ep += (assemble(error_p))*ddt

                        
            ### H1 seminorms
            # set the Expression function space
            # u_e2 = Expression(u_e.cppcode, t=u_e.t, element=Pu_e.ufl_element())
            # p_e2 = Expression(p_e.cppcode, t=p_e.t, k=mat_params['perm'], element=Pp_e.ufl_element())
            # point_error_u = uh-u_e2
            # point_error_p = ph-p_e2
            # H1_error_u = \
            #     inner(grad(point_error_u[0]), grad(point_error_u[0]))*dx\
            #     + inner(grad(point_error_u[1]), grad(point_error_u[1]))*dx\
            #     + inner(grad(point_error_u[2]), grad(point_error_u[2]))*dx
            # H1_error_p = inner(grad(point_error_p), grad(point_error_p))*dx
            # Eu2 += assemble(H1_error_u, mesh=mesh_e)*ddt
            # Ep2 += assemble(H1_error_p, mesh=mesh_e)*ddt

            # ####### Using errornorm

            # ### L2 norm by dolfin function
            # Eu2 += errornorm(u_e, u_h, norm_type='L2', degree_rise=0)**2
            # Ep2 += errornorm(p_e, p_h, norm_type='L2', degree_rise=0)**2

            ### H1 norm by dolfin function with projection
            u_h = project(uh, Pu)
            p_h = project(ph, Pp)
            Eu3 += errornorm(u_e, u_h, norm_type='H1', degree_rise=0)**2*ddt
            Ep3 += errornorm(p_e, p_h, norm_type='H1', degree_rise=0)**2*ddt

            ### H1 norm by direct vector manipulation (NOT WORKING)
            # uph = Function(V)
            # uph.vector()[:] = (up_2.vector().array()-up_1.vector().array())/(t_2-t_1)*(t-t_1) + up_1.vector().array()
            # pdb.set_trace()
            # Eu4 = errornorm(u_e, uph.sub(0), norm_type='H1', degree_rise=0)**2*ddt
            # Ep4 = errornorm(p_e, uph.sub(1), norm_type='H1', degree_rise=0)**2*ddt
        
    Eu = sqrt(Eu)
    Ep = sqrt(Ep)

    Eu2 = sqrt(Eu2)
    Ep2 = sqrt(Ep2)

    Eu3 = sqrt(Eu3)
    Ep3 = sqrt(Ep3)

    Eu4 = sqrt(Eu4)
    Ep4 = sqrt(Ep4)

    return [Eu, Eu2, Eu3, Eu4], [Ep, Ep2, Ep3, Ep4]

