// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_ADMMnet_RCPPEXPORTS_H_GEN_
#define RCPP_ADMMnet_RCPPEXPORTS_H_GEN_

#include <RcppEigen.h>
#include <Rcpp.h>

namespace ADMMnet {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("ADMMnet", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("ADMMnet", "_ADMMnet_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in ADMMnet");
            }
        }
    }

    inline List scaleC(Eigen::MatrixXd X) {
        typedef SEXP(*Ptr_scaleC)(SEXP);
        static Ptr_scaleC p_scaleC = NULL;
        if (p_scaleC == NULL) {
            validateSignature("List(*scaleC)(Eigen::MatrixXd)");
            p_scaleC = (Ptr_scaleC)R_GetCCallable("ADMMnet", "_ADMMnet_scaleC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_scaleC(Shield<SEXP>(Rcpp::wrap(X)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List OmegaC(Eigen::MatrixXd& Omega, Eigen::VectorXi& sgn) {
        typedef SEXP(*Ptr_OmegaC)(SEXP,SEXP);
        static Ptr_OmegaC p_OmegaC = NULL;
        if (p_OmegaC == NULL) {
            validateSignature("List(*OmegaC)(Eigen::MatrixXd&,Eigen::VectorXi&)");
            p_OmegaC = (Ptr_OmegaC)R_GetCCallable("ADMMnet", "_ADMMnet_OmegaC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_OmegaC(Shield<SEXP>(Rcpp::wrap(Omega)), Shield<SEXP>(Rcpp::wrap(sgn)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List OmegaSC(Eigen::SparseMatrix<double>& OmegaS, Eigen::VectorXi& sgn) {
        typedef SEXP(*Ptr_OmegaSC)(SEXP,SEXP);
        static Ptr_OmegaSC p_OmegaSC = NULL;
        if (p_OmegaSC == NULL) {
            validateSignature("List(*OmegaSC)(Eigen::SparseMatrix<double>&,Eigen::VectorXi&)");
            p_OmegaSC = (Ptr_OmegaSC)R_GetCCallable("ADMMnet", "_ADMMnet_OmegaSC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_OmegaSC(Shield<SEXP>(Rcpp::wrap(OmegaS)), Shield<SEXP>(Rcpp::wrap(sgn)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline double maxLambdaLmC(Eigen::MatrixXd X, Eigen::VectorXd y, double alpha, Eigen::VectorXd wbeta, int N0) {
        typedef SEXP(*Ptr_maxLambdaLmC)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_maxLambdaLmC p_maxLambdaLmC = NULL;
        if (p_maxLambdaLmC == NULL) {
            validateSignature("double(*maxLambdaLmC)(Eigen::MatrixXd,Eigen::VectorXd,double,Eigen::VectorXd,int)");
            p_maxLambdaLmC = (Ptr_maxLambdaLmC)R_GetCCallable("ADMMnet", "_ADMMnet_maxLambdaLmC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_maxLambdaLmC(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(wbeta)), Shield<SEXP>(Rcpp::wrap(N0)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline Eigen::VectorXd cvTrimLmC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF) {
        typedef SEXP(*Ptr_cvTrimLmC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cvTrimLmC p_cvTrimLmC = NULL;
        if (p_cvTrimLmC == NULL) {
            validateSignature("Eigen::VectorXd(*cvTrimLmC)(Eigen::VectorXd,int,int,Eigen::VectorXi,Eigen::MatrixXd,Eigen::VectorXd,int)");
            p_cvTrimLmC = (Ptr_cvTrimLmC)R_GetCCallable("ADMMnet", "_ADMMnet_cvTrimLmC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cvTrimLmC(Shield<SEXP>(Rcpp::wrap(beta)), Shield<SEXP>(Rcpp::wrap(nn)), Shield<SEXP>(Rcpp::wrap(nn2)), Shield<SEXP>(Rcpp::wrap(loco)), Shield<SEXP>(Rcpp::wrap(XF)), Shield<SEXP>(Rcpp::wrap(yF)), Shield<SEXP>(Rcpp::wrap(NF)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Eigen::VectorXd >(rcpp_result_gen);
    }

    inline List EnetLmC(Eigen::MatrixXd X, Eigen::VectorXd y, double alpha, Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta, int isd, int p, int N0, double thresh, int maxit) {
        typedef SEXP(*Ptr_EnetLmC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_EnetLmC p_EnetLmC = NULL;
        if (p_EnetLmC == NULL) {
            validateSignature("List(*EnetLmC)(Eigen::MatrixXd,Eigen::VectorXd,double,Eigen::VectorXd,int,int,Eigen::VectorXd,int,int,int,double,int)");
            p_EnetLmC = (Ptr_EnetLmC)R_GetCCallable("ADMMnet", "_ADMMnet_EnetLmC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_EnetLmC(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(nlambda)), Shield<SEXP>(Rcpp::wrap(ilambda)), Shield<SEXP>(Rcpp::wrap(wbeta)), Shield<SEXP>(Rcpp::wrap(isd)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(N0)), Shield<SEXP>(Rcpp::wrap(thresh)), Shield<SEXP>(Rcpp::wrap(maxit)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List cvEnetLmC(Eigen::MatrixXd X, Eigen::VectorXd y, double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, int N, int p, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF) {
        typedef SEXP(*Ptr_cvEnetLmC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cvEnetLmC p_cvEnetLmC = NULL;
        if (p_cvEnetLmC == NULL) {
            validateSignature("List(*cvEnetLmC)(Eigen::MatrixXd,Eigen::VectorXd,double,Eigen::VectorXd,int,Eigen::VectorXd,int,int,double,int,Eigen::MatrixXd,Eigen::VectorXd,int)");
            p_cvEnetLmC = (Ptr_cvEnetLmC)R_GetCCallable("ADMMnet", "_ADMMnet_cvEnetLmC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cvEnetLmC(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(nlambda)), Shield<SEXP>(Rcpp::wrap(wbeta)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(thresh)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(XF)), Shield<SEXP>(Rcpp::wrap(yF)), Shield<SEXP>(Rcpp::wrap(NF)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List NetLmC(Eigen::MatrixXd& X, Eigen::VectorXd& y, double alpha, Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta, Eigen::SparseMatrix<double>& Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj, int isd, int p, int N0, double thresh, int maxit) {
        typedef SEXP(*Ptr_NetLmC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_NetLmC p_NetLmC = NULL;
        if (p_NetLmC == NULL) {
            validateSignature("List(*NetLmC)(Eigen::MatrixXd&,Eigen::VectorXd&,double,Eigen::VectorXd,int,int,Eigen::VectorXd,Eigen::SparseMatrix<double>&,Eigen::MatrixXd,Eigen::VectorXi,int,int,int,double,int)");
            p_NetLmC = (Ptr_NetLmC)R_GetCCallable("ADMMnet", "_ADMMnet_NetLmC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_NetLmC(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(nlambda)), Shield<SEXP>(Rcpp::wrap(ilambda)), Shield<SEXP>(Rcpp::wrap(wbeta)), Shield<SEXP>(Rcpp::wrap(Omega)), Shield<SEXP>(Rcpp::wrap(loc)), Shield<SEXP>(Rcpp::wrap(nadj)), Shield<SEXP>(Rcpp::wrap(isd)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(N0)), Shield<SEXP>(Rcpp::wrap(thresh)), Shield<SEXP>(Rcpp::wrap(maxit)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List cvNetLmC(Eigen::MatrixXd& X, Eigen::VectorXd& y, double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, Eigen::SparseMatrix<double>& Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj, int N, int p, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF) {
        typedef SEXP(*Ptr_cvNetLmC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cvNetLmC p_cvNetLmC = NULL;
        if (p_cvNetLmC == NULL) {
            validateSignature("List(*cvNetLmC)(Eigen::MatrixXd&,Eigen::VectorXd&,double,Eigen::VectorXd,int,Eigen::VectorXd,Eigen::SparseMatrix<double>&,Eigen::MatrixXd,Eigen::VectorXi,int,int,double,int,Eigen::MatrixXd,Eigen::VectorXd,int)");
            p_cvNetLmC = (Ptr_cvNetLmC)R_GetCCallable("ADMMnet", "_ADMMnet_cvNetLmC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cvNetLmC(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(nlambda)), Shield<SEXP>(Rcpp::wrap(wbeta)), Shield<SEXP>(Rcpp::wrap(Omega)), Shield<SEXP>(Rcpp::wrap(loc)), Shield<SEXP>(Rcpp::wrap(nadj)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(thresh)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(XF)), Shield<SEXP>(Rcpp::wrap(yF)), Shield<SEXP>(Rcpp::wrap(NF)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline double maxLambdaCoxC(Eigen::MatrixXd X, Eigen::VectorXd tevent, int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, double alpha, Eigen::VectorXd wbeta, int N0) {
        typedef SEXP(*Ptr_maxLambdaCoxC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_maxLambdaCoxC p_maxLambdaCoxC = NULL;
        if (p_maxLambdaCoxC == NULL) {
            validateSignature("double(*maxLambdaCoxC)(Eigen::MatrixXd,Eigen::VectorXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,double,Eigen::VectorXd,int)");
            p_maxLambdaCoxC = (Ptr_maxLambdaCoxC)R_GetCCallable("ADMMnet", "_ADMMnet_maxLambdaCoxC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_maxLambdaCoxC(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(tevent)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(nevent)), Shield<SEXP>(Rcpp::wrap(nevent1)), Shield<SEXP>(Rcpp::wrap(loc1)), Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(wbeta)), Shield<SEXP>(Rcpp::wrap(N0)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double pletaCm(Eigen::VectorXd& xb, Eigen::VectorXd& exb, Eigen::VectorXi& nevent, Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n, int& ifast, int& itwo) {
        typedef SEXP(*Ptr_pletaCm)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_pletaCm p_pletaCm = NULL;
        if (p_pletaCm == NULL) {
            validateSignature("double(*pletaCm)(Eigen::VectorXd&,Eigen::VectorXd&,Eigen::VectorXi&,Eigen::VectorXi&,Eigen::VectorXi&,int&,int&,int&)");
            p_pletaCm = (Ptr_pletaCm)R_GetCCallable("ADMMnet", "_ADMMnet_pletaCm");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_pletaCm(Shield<SEXP>(Rcpp::wrap(xb)), Shield<SEXP>(Rcpp::wrap(exb)), Shield<SEXP>(Rcpp::wrap(nevent)), Shield<SEXP>(Rcpp::wrap(nevent1)), Shield<SEXP>(Rcpp::wrap(loc1)), Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(ifast)), Shield<SEXP>(Rcpp::wrap(itwo)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline Eigen::VectorXd cvTrimCoxC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco, Eigen::MatrixXd XF, int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF, Eigen::MatrixXd X, int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, int ifast, int itwo) {
        typedef SEXP(*Ptr_cvTrimCoxC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cvTrimCoxC p_cvTrimCoxC = NULL;
        if (p_cvTrimCoxC == NULL) {
            validateSignature("Eigen::VectorXd(*cvTrimCoxC)(Eigen::VectorXd,int,int,Eigen::VectorXi,Eigen::MatrixXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,Eigen::MatrixXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,int,int)");
            p_cvTrimCoxC = (Ptr_cvTrimCoxC)R_GetCCallable("ADMMnet", "_ADMMnet_cvTrimCoxC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cvTrimCoxC(Shield<SEXP>(Rcpp::wrap(beta)), Shield<SEXP>(Rcpp::wrap(nn)), Shield<SEXP>(Rcpp::wrap(nn2)), Shield<SEXP>(Rcpp::wrap(loco)), Shield<SEXP>(Rcpp::wrap(XF)), Shield<SEXP>(Rcpp::wrap(NF)), Shield<SEXP>(Rcpp::wrap(neventF)), Shield<SEXP>(Rcpp::wrap(nevent1F)), Shield<SEXP>(Rcpp::wrap(loc1F)), Shield<SEXP>(Rcpp::wrap(nF)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(nevent)), Shield<SEXP>(Rcpp::wrap(nevent1)), Shield<SEXP>(Rcpp::wrap(loc1)), Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(ifast)), Shield<SEXP>(Rcpp::wrap(itwo)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Eigen::VectorXd >(rcpp_result_gen);
    }

    inline List EnetCoxC(Eigen::MatrixXd X, Eigen::VectorXd tevent, double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, int p, int N0, double thresh, int maxit, int ifast) {
        typedef SEXP(*Ptr_EnetCoxC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_EnetCoxC p_EnetCoxC = NULL;
        if (p_EnetCoxC == NULL) {
            validateSignature("List(*EnetCoxC)(Eigen::MatrixXd,Eigen::VectorXd,double,Eigen::VectorXd,int,Eigen::VectorXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,int,int,double,int,int)");
            p_EnetCoxC = (Ptr_EnetCoxC)R_GetCCallable("ADMMnet", "_ADMMnet_EnetCoxC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_EnetCoxC(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(tevent)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(nlambda)), Shield<SEXP>(Rcpp::wrap(wbeta)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(nevent)), Shield<SEXP>(Rcpp::wrap(nevent1)), Shield<SEXP>(Rcpp::wrap(loc1)), Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(N0)), Shield<SEXP>(Rcpp::wrap(thresh)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(ifast)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List cvEnetCoxC(Eigen::MatrixXd X, Eigen::VectorXd tevent, double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, int p, int N0, double thresh, int maxit, int ifast, Eigen::MatrixXd XF, int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF) {
        typedef SEXP(*Ptr_cvEnetCoxC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cvEnetCoxC p_cvEnetCoxC = NULL;
        if (p_cvEnetCoxC == NULL) {
            validateSignature("List(*cvEnetCoxC)(Eigen::MatrixXd,Eigen::VectorXd,double,Eigen::VectorXd,int,Eigen::VectorXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,int,int,double,int,int,Eigen::MatrixXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int)");
            p_cvEnetCoxC = (Ptr_cvEnetCoxC)R_GetCCallable("ADMMnet", "_ADMMnet_cvEnetCoxC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cvEnetCoxC(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(tevent)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(nlambda)), Shield<SEXP>(Rcpp::wrap(wbeta)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(nevent)), Shield<SEXP>(Rcpp::wrap(nevent1)), Shield<SEXP>(Rcpp::wrap(loc1)), Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(N0)), Shield<SEXP>(Rcpp::wrap(thresh)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(ifast)), Shield<SEXP>(Rcpp::wrap(XF)), Shield<SEXP>(Rcpp::wrap(NF)), Shield<SEXP>(Rcpp::wrap(neventF)), Shield<SEXP>(Rcpp::wrap(nevent1F)), Shield<SEXP>(Rcpp::wrap(loc1F)), Shield<SEXP>(Rcpp::wrap(nF)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List NetCoxC(Eigen::MatrixXd& X, Eigen::VectorXd tevent, double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, Eigen::SparseMatrix<double>& Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj, int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, int p, int N0, double thresh, int maxit, int ifast) {
        typedef SEXP(*Ptr_NetCoxC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_NetCoxC p_NetCoxC = NULL;
        if (p_NetCoxC == NULL) {
            validateSignature("List(*NetCoxC)(Eigen::MatrixXd&,Eigen::VectorXd,double,Eigen::VectorXd,int,Eigen::VectorXd,Eigen::SparseMatrix<double>&,Eigen::MatrixXd,Eigen::VectorXi,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,int,int,double,int,int)");
            p_NetCoxC = (Ptr_NetCoxC)R_GetCCallable("ADMMnet", "_ADMMnet_NetCoxC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_NetCoxC(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(tevent)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(nlambda)), Shield<SEXP>(Rcpp::wrap(wbeta)), Shield<SEXP>(Rcpp::wrap(Omega)), Shield<SEXP>(Rcpp::wrap(loc)), Shield<SEXP>(Rcpp::wrap(nadj)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(nevent)), Shield<SEXP>(Rcpp::wrap(nevent1)), Shield<SEXP>(Rcpp::wrap(loc1)), Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(N0)), Shield<SEXP>(Rcpp::wrap(thresh)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(ifast)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List cvNetCoxC(Eigen::MatrixXd& X, Eigen::VectorXd tevent, double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, Eigen::SparseMatrix<double>& Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj, int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, int p, int N0, double thresh, int maxit, int ifast, Eigen::MatrixXd XF, int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF) {
        typedef SEXP(*Ptr_cvNetCoxC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cvNetCoxC p_cvNetCoxC = NULL;
        if (p_cvNetCoxC == NULL) {
            validateSignature("List(*cvNetCoxC)(Eigen::MatrixXd&,Eigen::VectorXd,double,Eigen::VectorXd,int,Eigen::VectorXd,Eigen::SparseMatrix<double>&,Eigen::MatrixXd,Eigen::VectorXi,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,int,int,double,int,int,Eigen::MatrixXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int)");
            p_cvNetCoxC = (Ptr_cvNetCoxC)R_GetCCallable("ADMMnet", "_ADMMnet_cvNetCoxC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cvNetCoxC(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(tevent)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(nlambda)), Shield<SEXP>(Rcpp::wrap(wbeta)), Shield<SEXP>(Rcpp::wrap(Omega)), Shield<SEXP>(Rcpp::wrap(loc)), Shield<SEXP>(Rcpp::wrap(nadj)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(nevent)), Shield<SEXP>(Rcpp::wrap(nevent1)), Shield<SEXP>(Rcpp::wrap(loc1)), Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(N0)), Shield<SEXP>(Rcpp::wrap(thresh)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(ifast)), Shield<SEXP>(Rcpp::wrap(XF)), Shield<SEXP>(Rcpp::wrap(NF)), Shield<SEXP>(Rcpp::wrap(neventF)), Shield<SEXP>(Rcpp::wrap(nevent1F)), Shield<SEXP>(Rcpp::wrap(loc1F)), Shield<SEXP>(Rcpp::wrap(nF)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

}

#endif // RCPP_ADMMnet_RCPPEXPORTS_H_GEN_
