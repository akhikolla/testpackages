#include <RcppEigen.h>
#include <R_ext/Lapack.h>

using namespace std;
using namespace Eigen;

int inverse(MatrixXd& A)
{
        int n = A.rows();
        int success_code;
        F77_CALL(dtrtri)("L", "N", &n, A.data(), &n, &success_code);
        return success_code;
}

