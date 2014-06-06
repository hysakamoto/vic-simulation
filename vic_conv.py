## Error calculations of the vic simulations.

import sys
sys.path.append("src") 

import matplotlib.pyplot as plt
import numpy as np

from dolfin import *
from error_calculations import errorCalc

# mesh-time

set_log_level(ERROR)
base_name = 'crn/'

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
sim_params = {'omega'   : 0.5,   # forward:0, backward: 1, C-N: 0.5 
              'T_total' : 10.0,
              'max_it'  : 64,
              'm_num'   : 32,
              'p_order' : 1
          }


## convergence parameters
# mesh_refinement = [1,2,4,8,16,32]
mesh_refinement = [1,2,4,8,16,32]
max_iterations = [4,8,16,32,64]

max_mer = 32
max_mit = 64

# errors
Eus_mesh = []
Eps_mesh = []
Eus_time = []
Eps_time = []

Eus2_mesh = []
Eps2_mesh = []
Eus2_time = []
Eps2_time = []

Eus3_mesh = []
Eps3_mesh = []
Eus3_time = []
Eps3_time = []

Eus4_mesh = []
Eps4_mesh = []
Eus4_time = []
Eps4_time = []


# mesh convergence
print 'mesh convergence analysis'
sim_params['m_num'] = max_mer
sim_params['max_it'] = max_mit
for mer in mesh_refinement:
    print 'mesh #: ', mer

    sim_params['m_num'] = mer
    Eus, Eps\
        = errorCalc(base_name, max_mer, max_mit, sim_params, mat_params)

    print Eus
    print Eps

    Eus_mesh.append(Eus[0])
    Eps_mesh.append(Eps[0])
    Eus2_mesh.append(Eus[1])
    Eps2_mesh.append(Eps[1])
    Eus3_mesh.append(Eus[2])
    Eps3_mesh.append(Eps[2])
    Eus4_mesh.append(Eus[3])
    Eps4_mesh.append(Eps[3])

# time convergence
print 'time step convergence analysis'
sim_params['m_num'] = max_mer
sim_params['max_it'] = max_mit
for mit in max_iterations:
    print 'iteration #: ', mit

    sim_params['max_it'] = mit
    Eus, Eps\
        = errorCalc(base_name, max_mer, max_mit, sim_params, mat_params)

    print Eus
    print Eps

    Eus_time.append(Eus[0])
    Eps_time.append(Eps[0])
    Eus2_time.append(Eus[1])
    Eps2_time.append(Eps[1])
    Eus3_time.append(Eus[2])
    Eps3_time.append(Eps[2])
    Eus4_time.append(Eus[3])
    Eps4_time.append(Eps[3])

# set the parameters back
sim_params['m_num'] = max_mer
sim_params['max_it'] = max_mit


##### PLOT RESULTS ######

# mesh-L2
plt.loglog(mesh_refinement, Eus2_mesh, '-o')
plt.loglog(mesh_refinement, Eps2_mesh, '-o')
plt.legend(('displacement', 'pressure'))
plt.title('mesh refinement vs L2 error')
plt.savefig(base_name + '/' + 'L2_errors_mesh.jpg')
plt.show()

# mesh-H1
plt.loglog(mesh_refinement, Eus3_mesh, '-o')
plt.loglog(mesh_refinement, Eps3_mesh, '-o')
plt.legend(('displacement', 'pressure'))
plt.title('mesh refinement vs H1 error')
plt.savefig(base_name + '/' + 'H1_errors_mesh.jpg')
plt.show()

# time-L2
plt.loglog(max_iterations, Eus2_time, '-o')
plt.loglog(max_iterations, Eps2_time, '-o')
plt.legend(('displacement', 'pressure'))
plt.title('time steps vs L2 error')
plt.savefig(base_name + '/' + 'L2_errors_time.jpg')
plt.show()

# time-H1
plt.loglog(max_iterations, Eus3_time, '-o')
plt.loglog(max_iterations, Eps3_time, '-o')
plt.legend(('displacement', 'pressure'))
plt.title('time steps vs H1 error')
plt.savefig(base_name + '/' + 'H1_errors_time.jpg')
plt.show()


## Convergence Rates

EuL2_mesh_rate = (np.log(Eus2_mesh[0])-np.log(Eus2_mesh[-1]))/\
                 (np.log(mesh_refinement[0])-np.log(mesh_refinement[-1]))
EpL2_mesh_rate = (np.log(Eps2_mesh[0])-np.log(Eps2_mesh[-1]))/\
                 (np.log(mesh_refinement[0])-np.log(mesh_refinement[-1]))
EuL2_time_rate = (np.log(Eus2_time[0])-np.log(Eus2_time[-1]))/\
                 (np.log(max_iterations[0])-np.log(max_iterations[-1]))
EpL2_time_rate = (np.log(Eps2_time[0])-np.log(Eps2_time[-1]))/\
                 (np.log(max_iterations[0])-np.log(max_iterations[-1]))
EuH1_mesh_rate = (np.log(Eus3_mesh[0])-np.log(Eus3_mesh[-2]))/\
                 (np.log(mesh_refinement[0])-np.log(mesh_refinement[-2]))
EpH1_mesh_rate = (np.log(Eps3_mesh[0])-np.log(Eps3_mesh[-2]))/\
                 (np.log(mesh_refinement[0])-np.log(mesh_refinement[-2]))
EuH1_time_rate = (np.log(Eus3_time[0])-np.log(Eus3_time[-1]))/\
                 (np.log(max_iterations[0])-np.log(max_iterations[-1]))
EpH1_time_rate = (np.log(Eps3_time[0])-np.log(Eps3_time[-1]))/\
                 (np.log(max_iterations[0])-np.log(max_iterations[-1]))

print 'Mesh Displacement L2 error convergence rate = ', \
    EuL2_mesh_rate
print 'Mesh Pressure L2 error convergence rate = ', \
    EpL2_mesh_rate
print 'Time Displacement L2 error convergence rate = ', \
    EuL2_time_rate
print 'Time Pressure L2 error convergence rate = ', \
    EpL2_time_rate
print 'Mesh Displacement H1 error convergence rate = ', \
    EuH1_mesh_rate
print 'Mesh Pressure H1 error convergence rate = ', \
    EpH1_mesh_rate
print 'Time Displacement H1 error convergence rate = ', \
    EuH1_time_rate
print 'Time Pressure H1 error convergence rate = ', \
    EpH1_time_rate


## Output to file
with open(base_name + '/errors.py', 'w') as f:
    f.write('EuL2_mesh = ' + str(Eus2_mesh)+'\n')
    f.write('EpL2_mesh = ' + str(Eps2_mesh)+'\n')
    f.write('EuH1_mesh = ' + str(Eus3_mesh)+'\n')
    f.write('EpH1_mesh = ' + str(Eps3_mesh)+'\n')
    f.write('EuL2_time = ' + str(Eus2_time)+'\n')
    f.write('EpL2_time = ' + str(Eps2_time)+'\n')
    f.write('EuH1_time = ' + str(Eus3_time)+'\n')
    f.write('EpH1_time = ' + str(Eps3_time)+'\n')

    f.write('\n')

    f.write('EuL2_mesh_rate = ' + str(EuL2_mesh_rate)+'\n')
    f.write('EpL2_mesh_rate = ' + str(EpL2_mesh_rate)+'\n')
    f.write('EuL2_time_rate = ' + str(EuL2_time_rate)+'\n')
    f.write('EpL2_time_rate = ' + str(EpL2_time_rate)+'\n')
    f.write('EuH1_mesh_rate = ' + str(EuH1_mesh_rate)+'\n')
    f.write('EpH1_mesh_rate = ' + str(EpH1_mesh_rate)+'\n')
    f.write('EuH1_time_rate = ' + str(EuH1_time_rate)+'\n')
    f.write('EpH1_time_rate = ' + str(EpH1_time_rate)+'\n')


