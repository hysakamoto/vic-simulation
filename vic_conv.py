from dolfin import *
from error_calculations import errorCalc

# mesh-time

base_name = 'bkn/'

# material parameters
mat_params = {'gamma'      : 0.0, 
              'tau'        : 1.0,
              'Ee'         : 100.0,
              'nu'         : 0.25,
              'perm'       : 1.0,
              'top_trac'   : (0.0,0.0,0.0),
              'body_force' : (0.0,0.0,0.0)
}

# simulation parameters
sim_params = {'omega'   : 1.0,   # forward:0, backward: 1, C-N: 0.5 
              'T_total' : 10.0,
              'max_it'  : 64,
              'm_num'   : 32,
              'p_order' : 1
          }


## convergence parameters
mesh_refinement = [1,2,4,8,16,32]
max_iterations = [1,2,4,8,16,32,64]

max_mer = 32
max_mit = 64

# errors
Eus_mesh = []
Eps_mesh = []
Eus_time = []
Eps_time = []

# mesh convergence
print 'mesh convergence analysis'
sim_params['m_num'] = max_mer
sim_params['max_it'] = max_mit
for mer in mesh_refinement:
    print mer

    sim_params['m_num'] = mer
    Eu, Ep = errorCalc(base_name, max_mer, max_mit, sim_params, mat_params)

    Eus_mesh.append(Eu)
    Eps_mesh.append(Ep)

# time convergence
print 'time step convergence analysis'
sim_params['m_num'] = max_mer
sim_params['max_it'] = max_mit
for mit in max_iterations:
    print mit

    sim_params['max_it'] = mit
    Eu, Ep = errorCalc(base_name, max_mer, max_mit, sim_params, mat_params)

    Eus_time.append(Eu)
    Eps_time.append(Ep)

# set the parameters back
sim_params['m_num'] = max_mer
sim_params['max_it'] = max_mit


##### PLOT RESULTS ######

# mesh
plt.loglog(mesh_refinement, Eus_mesh, '-o')
plt.loglog(mesh_refinement, Eps_mesh, '-o')
plt.legend(('displacement', 'pressure'))
plt.title('mesh refinement vs L2 error')
plt.savefig(base_name + '/' + 'L2_errors_mesh.jpg')
plt.show()

# time
plt.loglog(max_iterations, Eus_time, '-o')
plt.loglog(max_iterations, Eps_time, '-o')
plt.legend(('displacement', 'pressure'))
plt.title('time steps vs L2 error')
plt.savefig(base_name + '/' + 'L2_errors_time.jpg')
plt.show()


with open(base_name + '/L2_errros.py', 'w') as f:
    f.write('Eus_mesh = ' + str(Eus_mesh)+'\n')
    f.write('Eps_mesh = ' + str(Eps_mesh)+'\n')
    f.write('Eus_time = ' + str(Eus_time)+'\n')
    f.write('Eps_time = ' + str(Eps_time)+'\n')


print 'Mesh Displacement L2 error convergence rate = ', \
    (np.log(Eus_mesh[0])-np.log(Eus_mesh[5]))/\
    (np.log(mesh_refinement[0])-np.log(mesh_refinement[5]))
print 'Mesh Pressure L2 error convergence rate = ', \
    (np.log(Eps_mesh[0])-np.log(Eps_mesh[5]))/\
    (np.log(mesh_refinement[0])-np.log(mesh_refinement[5]))

print 'Time Displacement L2 error convergence rate = ', \
    (np.log(Eus_time[0])-np.log(Eus_time[6]))/\
    (np.log(max_iterations[0])-np.log(max_iterations[6]))
print 'Time Pressure L2 error convergence rate = ', \
    (np.log(Eps_time[0])-np.log(Eps_time[6]))/\
    (np.log(max_iterations[0])-np.log(max_iterations[6]))


# A = [0.016984484000945266, 0.003620986790733252, 0.0009605201341572476, 0.001231537706287947, 0.0013764656561359874, 0.001415489361134992, 0.0014254100429035416]
# B = [0.5022122432696464, 0.11371038131558306, 0.030826510023513484, 0.009615154207164462, 0.004518699535452111, 0.003383338520246388, 0.0031324720533711745]
# C = [0.08819447320545222, 0.04710652215528491, 0.016930969720366555, 0.004920670105373137, 0.0014468970724907848]
# D = [0.1920775348159219, 0.0969587566721823, 0.033902293135667554, 0.009841815798289924, 0.0029976417772129976]
