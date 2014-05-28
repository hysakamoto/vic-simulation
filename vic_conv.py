from dolfin import *
from manufactured_solutions import manufactured_solutions

sim_name = 'result'
p_order = 1

T_total = 10.0
dt      = 1    # time step
max_it  = int(T_total/dt)

# Material parameters
Ee, nu = 100.0, 0.25
mu, lmbda = Constant(Ee/(2*(1 + nu))), Constant(Ee*nu/((1 + nu)*(1 - 2*nu)))
perm    = 1.0

# get manufactured solutions
[u_e, p_e, v_e, source, body_force, tbars, gbars, u_initial, p_initial] \
    = manufactured_solutions(0.0, perm, mu, lmbda)



##### TIME-STEP CONVERGENCE #####
Eus1 = []
Eps1 = []
L2u1 = []
L2p1 = []


##  Finest Function Spaces
mesh_e = UnitCubeMesh(16,16,16)
Pu_e = VectorFunctionSpace(mesh_e, "Lagrange", p_order)  # space for displacements
Pp_e = FunctionSpace(mesh_e, "Lagrange", p_order)        # space for pressure
V_e  = MixedFunctionSpace([Pu_e,Pp_e])                    # mixed space


T_total = 10.0
max_its = [2,4,8,16,32,64]
max_it_num = max(max_its)
sim_basename = 'time/'
for i in range(len(max_its)):
    max_it = max_its[i]
    n_err_comp = max_it_num/max_it # number of error computation per time step
    dt = T_total/float(max_it)
    sim_name = sim_basename + str(max_it)

    # load mesh
    # mesh = Mesh(sim_name+'/mesh.xdmf')
    m_num = 16
    mesh = UnitCubeMesh(m_num, m_num, m_num)
    File(sim_name+'/mesh.xdmf') << mesh


    ##  Function Spaces
    Pu = VectorFunctionSpace(mesh, "Lagrange", p_order)  # space for displacements
    Pp = FunctionSpace(mesh, "Lagrange", p_order)        # space for pressure
    V  = MixedFunctionSpace([Pu,Pp])                    # mixed space


    ddt = dt/(n_err_comp)

    Eu = 0.0
    Ep = 0.0
    L2u = 0.0
    L2p = 0.0
    t = 0.0
    for tn in range(1,max_it+1):

        # load solutions
        up_1 = Function(V, sim_name + '/up_%d.xml'%(tn-1))
        u_1 = up_1.sub(0)
        p_1 = up_1.sub(1)
        
        up_2 = Function(V, sim_name + '/up_%d.xml'%(tn))
        u_2 = up_2.sub(0)
        p_2 = up_2.sub(1)

        t1 = Constant(t)
        t2 = Constant(t+dt)

        for j in range(n_err_comp):

            t += ddt
        
            # set time
            u_e.t = t
            p_e.t = t

            t_const = Constant(t)

            uh = (u_2-u_1)/(t2-t1)*(t_const-t1) + u_1
            ph = (p_2-p_1)/(t2-t1)*(t_const-t1) + p_1

            ### Error against exact solutions
            error_u = (uh-u_e)**2*dx
            Eu += sqrt(assemble(error_u))*ddt

            error_p = (ph-p_e)**2*dx
            Ep += sqrt(assemble(error_p))*ddt

            # L2 norm of these
            u_e_intp = interpolate(u_e, Pu_e)
            p_e_intp = interpolate(p_e, Pp_e)

            uesq = u_e_intp**2*dx
            pesq = p_e_intp**2*dx
            L2u += sqrt(assemble(uesq))*ddt
            L2p += sqrt(assemble(pesq))*ddt

            
    Eu = Eu
    Ep = Ep

    Eus1.append(Eu)
    Eps1.append(Ep)

    L2u1.append(L2u)
    L2p1.append(L2p)


##### MESH CONVERGENCE #####

##  Finest Function Spaces
mesh_e = UnitCubeMesh(16,16,16)
Pu_e = VectorFunctionSpace(mesh_e, "Lagrange", p_order)  # space for displacements
Pp_e = FunctionSpace(mesh_e, "Lagrange", p_order)        # space for pressure
V_e  = MixedFunctionSpace([Pu_e,Pp_e])                    # mixed space

### Run Simulation
u_max = []
Eus2 = []
Eps2 = []
L2u2 = []
L2p2 = []
m_nums = [1,2,4,8,16]
sim_basename = 'mesh/'

T_total = 10.0
max_it = 16
dt = T_total/float(max_it)

for i in range(len(m_nums)):
    m_num = m_nums[i]
    sim_name = sim_basename + str(m_num)

    # load mesh
    # mesh = Mesh(sim_name+'/mesh.xdmf')
    mesh = UnitCubeMesh(m_num, m_num, m_num)
    File(sim_name+'/mesh.xdmf') << mesh

    ##  Function Spaces
    Pu = VectorFunctionSpace(mesh, "Lagrange", p_order)  # space for displacements
    Pp = FunctionSpace(mesh, "Lagrange", p_order)        # space for pressure
    V  = MixedFunctionSpace([Pu,Pp])                    # mixed space

    Ep = 0.0
    Eu = 0.0
    L2u = 0.0
    L2p = 0.0
    t = dt
    for tn in range(1, max_it+1):

        # load solutions
        up = Function(V, sim_name + '/up_%d.xml'%tn)
        u_tent, p_tent = up.split(deepcopy=True) 

        u_e.t = t
        p_e.t = t

        # Explicit interpolation of u_e onto the same space as u:
        u_intp = interpolate(u_tent, Pu_e)
        u_e_intp = interpolate(u_e, Pu_e)
        error_u_intp = (u_intp - u_e_intp)**2*dx
        Eu += sqrt(assemble(error_u_intp))*dt

        p_intp = interpolate(p_tent, Pp_e)
        p_e_intp = interpolate(p_e, Pp_e)
        error_p_intp = (p_intp - p_e_intp)**2*dx
        Ep += sqrt(assemble(error_p_intp))*dt

        # L2 norm of these
        uesq = u_e_intp**2*dx
        pesq = p_e_intp**2*dx
        L2u += sqrt(assemble(uesq))*dt
        L2p += sqrt(assemble(pesq))*dt

        t += dt


    Eus2.append(Eu)
    Eps2.append(Ep)
    L2u2.append(L2u)
    L2p2.append(L2p)



# print Eus2
# print Eps2


##### PLOT RESULTS ######

relEu1 = [Eus1[i]/L2u1[i] for i in range(len(Eus1))]
relEp1 = [Eps1[i]/L2p1[i] for i in range(len(Eus1))]

plt.loglog(max_its,relEu1)
plt.loglog(max_its,relEp1)
plt.legend(('displacement', 'pressure'))
plt.show()


relEu2 = [Eus2[i]/L2u2[i] for i in range(len(Eus2))]
relEp2 = [Eps2[i]/L2p2[i] for i in range(len(Eus2))]

plt.loglog(m_nums,relEu2)
plt.loglog(m_nums,relEp2)
plt.legend(('displacement', 'pressure'))
plt.show()

