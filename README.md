## Observations 

* Mixed variable assignment does not get reflected unless explicitly
  define residual and jacobian at each time step.
* Backward Euler and Crank-Nicolson schemes implemented. 
* The computation is very large on a fine mesh --> 32\*32\*32 mesh and
  64 time steps simulation takes ~30min for 16*8 nodes on Stampede.

[FEniCS Q&A](http://fenicsproject.org/qa/)

## Convergence Analysis against Manufactured Solution

* Poroelastic material with compressible neo-Hookean solid and
  isotropic permeability.
* The analytical solutions are manufactured using MATLAB and converted
  to Python (Dolfin) readable script.
* The analytical solutions are derived in two ways: 1) from current
  configuration and 2) from initial configuration. The option *2*
  seems working but not *1* --> why?

### Computation of the error

1. Solve for **u_h** and p_h at each times step.
2. Get eacact values for **u** and p.
3. Compute the L2 error in space.
4. Compute the L2 error of the error in *3*.
5. Numerical integration of the errors are done on the finest mesh and
   on the smallest time step.

### TODO

1. Convergence analysis in time and apce.
2. Implementation of H1 error calucation and convergence analysis.
3. Converence rate in Backward Euler scheme and Crank-Nicholson scheme.

### Current and Initial Normal Tractions

1. Holzapfel: p111~113, (3.1)~(3.10) -- explains the conversion
   from/to the current/initial normal tractions.
       * The direction of the current and initial tractions are the
         same.
	   * The initial traction is defined by the Cauchy stress tensor
         and initial normal vector.
	   * The current traction is defined by the first Piola-Kirchhoff
         stress and current normal vector.
	   * The second PK stress needs to be converted to the first PK
         stress in order to use it for the traction boundary
         conditions.
2. Wood: p129~131, (4.30), (3.68) -- conversion from/to the
   current/initial normal tractions using area change term.
