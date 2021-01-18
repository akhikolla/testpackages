# 2019-12-12 sokol@insa-toulouse.fr

# build internal lib for sundials
SUNTOP=/usr/local/src/cvodes-5.0.0
MYTOP=$HOME/dev/R/rcpp-pkgs/r2sundials
mkdir -p $MYTOP/src/lib/

# copy sources
cp -a $SUNTOP/src/{cvodes,nvector,sundials,sunlinsol,sunmatrix,sunnonlinsol} $MYTOP/src/lib/
# remove fmod dirs
find $MYTOP/src/lib/ -type d -name fmod -exec rm -fr {} \;
find $MYTOP/src/lib/ -name CMakeLists.txt -exec rm -fr {} \;
# remove pthreads etc
rm -rf $MYTOP/src/lib/nvector/{cuda,manyvector,mpiplusx,openmp,openmpdev,parallel,parhyp,petsc,raja,trilinos}
rm -rf $MYTOP/src/lib/sunnonlinsol/petscsnes
rm -rf $MYTOP/src/lib/sunlinsol/{klu,superludist,superlumt}
rm -rf $MYTOP/src/lib/sunmatrix/slunrloc

( cd $MYTOP/inst/include/nvector/ && rm -rf $(ls -1 | grep -v serial) )

# copy includes
cp -a $SUNTOP/include/{cvodes,nvector,sundials,sunlinsol,sunmatrix,sunnonlinsol} $MYTOP/inst/include
cp -a $SUNTOP/build/include/sundials/sundials_config.h $MYTOP/inst/include/sundials/
