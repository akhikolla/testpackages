[![Travis-CI Build Status](https://travis-ci.org/nlmixrdevelopment/lbfgsb3c.svg?branch=master)](https://travis-ci.org/nlmixrdevelopment/lbfgsb3c)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/nlmixrdevelopment/lbfgsb3c?branch=master&svg=true)](https://ci.appveyor.com/project/nlmixrdevelopment/lbfgsb3c)

# libfgsb3c interface from C
This is the fork of the libfgsb3 from cran with the following differences:
- The return type has changed is is very similar to what `optim` returns
- Allows a direct C/C++ interface through a R registered function,
  similar to C interface to `optim` with 2 additional arguments.
- Allows adjustment of tolerances for minimization success.
- Added `xtolAtol` and `xtolRtol` minimization success criterion.
- Added `maxit` termination
