#include "rcor.h"
#include "rcor_bm.h"
#include "rcor_permtest.h"
#include "rcor_exacttest.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


static const R_CallMethodDef callMethods[] = {
    {"rcor", (DL_FUNC) &rcor, 3},
    {"rcor_matrix_classical", (DL_FUNC) &rcor_matrix_classical, 2},
    {"rcor_matrix_linear", (DL_FUNC) &rcor_matrix_linear, 2},
    {"rcor_matrix_exp", (DL_FUNC) &rcor_matrix_exp, 2},
    {"rcor_matrix_gauss", (DL_FUNC) &rcor_matrix_gauss, 2},
    {"rcor_matrix_epstol", (DL_FUNC) &rcor_matrix_epstol, 2},
    {"rcor_matrices_classical", (DL_FUNC) &rcor_matrices_classical, 4},
    {"rcor_matrices_linear", (DL_FUNC) &rcor_matrices_linear, 4},
    {"rcor_matrices_exp", (DL_FUNC) &rcor_matrices_exp, 4},
    {"rcor_matrices_gauss", (DL_FUNC) &rcor_matrices_gauss, 4},
    {"rcor_matrices_epstol", (DL_FUNC) &rcor_matrices_epstol, 4},
    {"rcor_exacttest", (DL_FUNC) &rcor_exacttest, 7},
    {"rcor_permtest", (DL_FUNC) &rcor_permtest, 7},
    {"permNextWrapper", (DL_FUNC) &permNextWrapper, 2},
    {NULL, NULL, 0}
};

extern "C"
{
    void attribute_visible R_init_rococo(DllInfo *info) {
		/* Register routines, allocate resources. */
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
    }

    void R_unload_rococo(DllInfo *info) {
		/* Release resources. */
    }
}
