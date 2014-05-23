## Observations 

* Mixed variable assignment does not get reflected unless explicitly
  define residual and jacobian at each time step.
* Backward Euler implemented. 

[FEniCS Q&A](http://fenicsproject.org/qa/)


## Convergence Analysis against Manufactured Solution

* Uniaxial stretch of a box.
* Stretch is time-dependent and defined by a quadratic function.
* Poroelastic material with compressible neo-Hookean solid and
  isotropic permeability.
* Deformation is defined so that there is no (local) volume change.
* Initial displacement and pressure = 0.

### Computation of the error

1. Solve for **u_h** and p_h at each times step.
2. Get eacact values for **u** and p.
3. Compute the (relative) error.
4. Integrate it over the (initial) domain.
5. Repeat it for every time step.
6. Integrate the error over time.
7. **The integration over space need to be done on the fine mesh!!!**


### TODO

1. Check the nonlinear boundary conditions. - read Wood
2. Implement the nonlinear boundary conditions.
3. Convergence analysis.
4. Time- and mesh- convergence.

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
