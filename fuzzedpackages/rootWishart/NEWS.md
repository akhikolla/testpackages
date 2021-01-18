## rootWishart 0.4.1
 - Fix NOTE from R-devel
 - Add testing suite

## rootWishart 0.4.0
 - Register native routines
 - Parameter `mprec` has been removed. The type of precision is now controlled by the parameter `type`.
   - `type = 'double'` gives double precision
   - `type = 'multi'` gives multiprecision
 - The default behaviour when `type` is unspecified is to decide adaptively based on the input parameters which type of precision to choose.

## rootWishart 0.3.0
 - First release
