#include <RcppArmadillo.h>
#include "RProblem.h"

RCPP_MODULE(ManifoldOptim_module) {
	class_<RProblem>("RProblem")
	.constructor<Function,Function,Function>()
	.constructor<Function,Function>()
	.constructor<Function>()
	.method("objFun", &RProblem::objFun)
	.method("gradFun", &RProblem::gradFun)
	.method("hessEtaFun", &RProblem::hessEtaFun)
	;
}
