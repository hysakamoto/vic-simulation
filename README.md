## Observations ##

* Mixed variable assignment does not get reflected unless explicitly
  define residual and jacobian at each time step.
* Backward Euler implemented. 

[FEniCS Q&A](http://fenicsproject.org/qa/)


## Convergence Analysis against Manufactured Solution ##

* Uniaxial stretch of a box.
* Stretch is time-dependent and defined by a quadratic function.
* Poroelastic material with compressible neo-Hookean solid and
  isotropic permeability.
* Deformation is defined so that there is no (local) volume change.
* Initial displacement and pressure = 0.
