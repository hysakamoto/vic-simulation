from dolfin import *
from manufactured_solutions import manufactured_solutions

base_base_name = 'cn'
p_order = 1

print "Generating manufactured solutions from MATLAB output..."
execfile("conveq.py")

T_total = 10.0
dt      = 1    # time step
max_it  = int(T_total/dt)

# Material parameters
Ee, nu = 100.0, 0.25
mu, lmbda = Constant(Ee/(2*(1 + nu))), Constant(Ee*nu/((1 + nu)*(1 - 2*nu)))
perm    = 1.0

# get manufactured solutions
[u_e, p_e, v_e, source, body_force, tbars, gbars, u_initial, p_initial, v_initial] \
    = manufactured_solutions(0.0, perm, mu, lmbda)

##### TIME-STEP CONVERGENCE #####
print "time-step convergence analysis"
Eus1 = []
Eps1 = []

##  Finest Function Spaces
mesh_e = UnitCubeMesh(32,32,32)
Pu_e = VectorFunctionSpace(mesh_e, "Lagrange", p_order)  # space for displacements
Pp_e = FunctionSpace(mesh_e, "Lagrange", p_order)        # space for pressure
V_e  = MixedFunctionSpace([Pu_e,Pp_e])                    # mixed space


T_total = 10.0
max_its = [1, 2,4,8,16,32,64]
max_it_num = max(max_its)
sim_basename = base_base_name+ '/' +'time/'
for i in range(len(max_its)):
    print i

    max_it = max_its[i]
    n_err_comp = max_it_num/max_it # number of error computation per time step
    dt = T_total/float(max_it)
    sim_name = sim_basename + str(max_it)

    # load mesh
    # mesh = Mesh(sim_name+'/mesh.xdmf')
    m_num = 32
    mesh = UnitCubeMesh(m_num, m_num, m_num)
    File(sim_name+'/mesh.xdmf') << mesh


    ##  Function Spaces
    Pu = VectorFunctionSpace(mesh, "Lagrange", p_order)  # space for displacements
    Pp = FunctionSpace(mesh, "Lagrange", p_order)        # space for pressure
    V  = MixedFunctionSpace([Pu,Pp])                    # mixed space

    Eu = 0.0
    Ep = 0.0
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

        # n_err_comp = 2
        ddt = dt/(n_err_comp)

        for j in range(n_err_comp):

            t += ddt
        
            # set time
            u_e.t = t-ddt/2.0
            p_e.t = t-ddt/2.0

            t_const = Constant(t-ddt/2.0)

            uh = (u_2-u_1)/(t2-t1)*(t_const-t1) + u_1
            ph = (p_2-p_1)/(t2-t1)*(t_const-t1) + p_1

            ### Error against exact solutions
            error_u = (uh-u_e)**2*dx
            Eu += (assemble(error_u))*ddt

            error_p = (ph-p_e)**2*dx
            Ep += (assemble(error_p))*ddt

    Eus1.append(sqrt(Eu))
    Eps1.append(sqrt(Ep))


print Eus1
print Eps1


##### MESH CONVERGENCE #####
print "mesh convergence analysis"

##  Finest Function Spaces
mesh_e = UnitCubeMesh(32,32,32)
Pu_e = VectorFunctionSpace(mesh_e, "Lagrange", p_order)  # space for displacements
Pp_e = FunctionSpace(mesh_e, "Lagrange", p_order)        # space for pressure
V_e  = MixedFunctionSpace([Pu_e,Pp_e])                    # mixed space

### Run Simulation
u_max = []
Eus2 = []
Eps2 = []
m_nums = [1,2,4,8,16,32]
sim_basename = base_base_name+'/' + 'mesh/'

# L2 norm of solutions
L2u = 0.0
L2p = 0.0

T_total = 10.0
max_it = 64
dt = T_total/float(max_it)

for i in range(len(m_nums)):
    print i
    
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
        Eu += (assemble(error_u_intp))*dt

        p_intp = interpolate(p_tent, Pp_e)
        p_e_intp = interpolate(p_e, Pp_e)
        error_p_intp = (p_intp - p_e_intp)**2*dx
        Ep += (assemble(error_p_intp))*dt

        # L2 norm of these
        if i==0:
            uesq = u_e_intp**2*dx
            pesq = p_e_intp**2*dx
            L2u += (assemble(uesq))*dt
            L2p += (assemble(pesq))*dt

        t += dt

    Eus2.append(sqrt(Eu))
    Eps2.append(sqrt(Ep))

print Eus2
print Eps2

# ##### PLOT RESULTS ######

# time
relEu1 = [Eus1[i]/L2u for i in range(len(Eus1))]
relEp1 = [Eps1[i]/L2p for i in range(len(Eus1))]

plt.loglog(max_its,relEu1, '-o')
plt.loglog(max_its,relEp1, '-o')
plt.legend(('displacement', 'pressure'))
plt.savefig(base_base_name + '/'+'timestep-up-conv.jpg')
plt.show()

# mesh
relEu2 = [Eus2[i]/L2u for i in range(len(Eus2))]
relEp2 = [Eps2[i]/L2p for i in range(len(Eus2))]

plt.loglog(m_nums,relEu2,'-o')
plt.loglog(m_nums,relEp2,'-o')
plt.legend(('displacement', 'pressure'))
plt.savefig(base_base_name + '/' + 'mesh-up-conv.jpg')
plt.show()

with open(base_base_name + '/convergence.txt', 'w') as f:
    f.write(str(Eus1)+'\n')
    f.write(str(Eps1)+'\n')
    f.write(str(Eus2)+'\n')
    f.write(str(Eps2)+'\n')


print (np.log(Eus1[0])-np.log(Eus1[3]))/(np.log(max_its[0])-np.log(max_its[3]))
print (np.log(Eps1[0])-np.log(Eps1[4]))/(np.log(max_its[0])-np.log(max_its[4]))


print (np.log(Eus2[0])-np.log(Eus2[5]))/(np.log(m_nums[0])-np.log(m_nums[5]))
print (np.log(Eps2[0])-np.log(Eps2[4]))/(np.log(m_nums[0])-np.log(m_nums[4]))





# A = [0.016984484000945266, 0.003620986790733252, 0.0009605201341572476, 0.001231537706287947, 0.0013764656561359874, 0.001415489361134992, 0.0014254100429035416]
# B = [0.5022122432696464, 0.11371038131558306, 0.030826510023513484, 0.009615154207164462, 0.004518699535452111, 0.003383338520246388, 0.0031324720533711745]
# C = [0.08819447320545222, 0.04710652215528491, 0.016930969720366555, 0.004920670105373137, 0.0014468970724907848]
# D = [0.1920775348159219, 0.0969587566721823, 0.033902293135667554, 0.009841815798289924, 0.0029976417772129976]
