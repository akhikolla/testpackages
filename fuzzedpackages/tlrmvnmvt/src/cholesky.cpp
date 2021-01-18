#include <RcppEigen.h>
#include <R_ext/Lapack.h>

using namespace std;
using namespace Eigen;

int cholesky(MatrixXd& A)
{
        int n = A.rows();
        int success_code;
        F77_CALL(dpotrf)("L", &n, A.data(), &n, &success_code);
        A.triangularView<StrictlyUpper>().setZero();
        return success_code;
}

