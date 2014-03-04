## Observations ##

* Minimization of total potential energy works well.
..* Compressible neo-Hookean
..* Nearly incompressible neo-Hookean
* Explicity deriving the PK2 strees tensor works well too.
..* inner(A,B) where A and B are 2-rank tensors is the tensor
contraction.
..* tr(A*B.T) is the tensor contraction too.
* Defining a strin energy density function and symbolically deriving
  PK2 tensor does work some time.
..* It works for material model for both C and E derivative cases.
..* It does/doesn't work for different cases.
..* Maybe the definition of different tensors

**Conclusion**
* The variable (tensor) needs to be defined as "variable" first before
  all the dependencies are constructed.
* Thus, Right C-G tensor C cannot be used as a differential for F
  because C is constructed from F.



