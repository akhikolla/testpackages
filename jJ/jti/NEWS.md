## jti v0.6.0 (2020-12-16)

 * There was a bug in the creation of the junction tree when calling Kruskals algorithm.
 
 * It is now possible to specify variables of interest in advance, such that we are 
 guaranteed to be able to query the joint pmf of these variables.
 * Some refactoring making compilation much faster. When potentials is assigned to
 a clique we no longer start by creating a unity table and then multiply. This was killing
 the advantage of the sparsity.
 
## jti v0.5.2 (2020-11-24)

 * A new way of entering evidence that is much more robust
 
## jti v0.5.1 (2020-11-13)

 * A more optimal triangulation method has been implemented, which in general leads to faster run time of the **jt** function.
 * Some mis-spelled words fixed

## jti v0.5.0 (2020-11-09)

 * Initial version
