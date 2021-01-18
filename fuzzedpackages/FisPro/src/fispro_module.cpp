#include <Rcpp.h>
#include <string>
#include "fis.h"

//' @title fis class
//' @name fis
//' @docType class
//' @description The fis class is the main class (Rcpp class) to open and use a Fuzzy Inference System.
//'
//' @field name the name of the FIS.
//' @field input_size the number of inputs in the FIS.
//' @field output_size the number of outputs in the FIS.
//'
//' @section Methods:
//' \describe{
//' \item{constructor}{argument: \code{fis_file} a FIS configuration file (can be built with FisPro software \url{https://www.fispro.org}). \cr return: an object of fis.}
//' \item{infer}{infers a value for each output from the input values \cr argument: \code{values} numerical vector of input values. \cr return: all infered outputs.}
//' \item{infer_output}{infers a value for a single output from the input values \cr argument: \code{values} numerical vector of input values. \cr argument: \code{output_number} the number of the output to infer. \cr return: the infered output.}
//' }
//' 
//' @references{
//'   \href{https://www.fispro.org/documentation/en/inline-help/node39.html}{Fuzzy Logic Elementary Glossary}
//' }
//'
//' @examples
//' # build the FIS
//' fis_file <- system.file("extdata", "test.fis", package = "FisPro")
//' fis <- new(fis, fis_file)
//'
//' # infer all outputs
//' infered <- fis$infer(c(0.25, 0.75))
//' #infer first output
//' infered_output1 <- fis$infer_output(c(0.25, 0.75), 0)
//' #infer second output
//' infered_output2 <- fis$infer_output(c(0.25, 0.75), 1)
//'
//' @export fis
class fis_wrapper {

private:
	FIS *fis;

	void check_infer_values(Rcpp::NumericVector values) {
		if(values.size() != fis->GetNbIn())
			Rcpp::stop("values must be equal to input size");
	}

public:
	fis_wrapper() : fis(NULL) {
	    Rcpp::stop("fis default constructor not allowed");
	}

	fis_wrapper(const char *fis_file) : fis(new FIS(fis_file)) {}

	~fis_wrapper() {
		delete fis;
		fis = NULL;
	}

	const char *get_name() const {
		return fis->Name;
	}

	void set_name(const char *name) {
		fis->SetName(name);
	}

	int get_input_size() {
		return fis->GetNbIn();
	}

	int get_output_size() {
		return fis->GetNbOut();
	}

	Rcpp::NumericVector infer(Rcpp::NumericVector values) {
		check_infer_values(values);
		fis->Infer(values.begin());
		return Rcpp::NumericVector(fis->OutValue, fis->OutValue + fis->GetNbOut());
	}

	double infer_output(Rcpp::NumericVector values, int output_number) {
		check_infer_values(values);
		if((output_number < 0) || (output_number >= fis->GetNbOut()))
			Rcpp::stop("output_number must be in range [0, output size)");
		fis->Infer(values.begin(), output_number);
		return fis->OutValue[output_number];
	}
};

//' @title mf class
//' @docType class
//' @name mf
//' @description Base class of all MF classes. This class cannot be instantiate. Use derived classes \link{mf_triangular}, \link{mf_trapezoidal}, \link{mf_trapezoidal_inf} or \link{mf_trapezoidal_sup}
//'
//' @field name the name of the mf.
//'
//' @section Methods:
//' \describe{
//' \item{\code{degree(value)}}{compute the membership degree of a numerical value. \cr argument: \code{value} numerical value to compute the membership degree. \cr return: the membership degree.}
//' }
//' 
//' @references{
//'   \href{https://www.fispro.org/documentation/en/inline-help/node39.html}{Fuzzy Logic Elementary Glossary}
//' }
//'
//' @export mf
class mf_wrapper {

private:
	MF* mf;

protected:
	mf_wrapper(MF* mf) : mf(mf) {}

public:
	mf_wrapper() : mf(NULL) {
		Rcpp::stop("mf class is not instantiable, use derived classes");
	}

	virtual ~mf_wrapper() {
		delete mf;
		mf = NULL;
	}

	const char *get_name() const {
		return mf->Name;
	}

	void set_name(const char *name) {
		mf->SetName(name);
	}

	double get_degree(double value) const {
		return mf->GetDeg(value);
	}
};

//' @title mf_triangular class
//' @name mf_triangular
//' @docType class
//' @description Class to build triangular MF.
//'
//' @section Inherits:
//' mf_triangular class inherits all fields and methods of \link{mf} class.
//'
//' @section Methods:
//' \describe{
//' \item{\code{constructor(lower_support, kernel, upper_support)}}{argument: \code{lower_support} numerical lower value of support. \cr argument: \code{kernel} numerical value of kernel. \cr argument: \code{upper_support} numerical upper value of support. \cr return: an object of mf_triangular.}
//' }
//'
//' @examples
//' mf <- new(mf_triangular, 0, 1, 2)
//' mf$degree(0.5)
//'
//' @export mf_triangular
class mf_triangular_wrapper : public mf_wrapper {

private:
	MFTRI *mf_triangular;

public:
	mf_triangular_wrapper() : mf_triangular_wrapper(NULL) {
	    Rcpp::stop("mf_triangular default constructor not allowed");
	}
	mf_triangular_wrapper(double lower_support, double kernel, double upper_support) : mf_triangular_wrapper(new MFTRI(lower_support, kernel, upper_support)) {}

	mf_triangular_wrapper(MFTRI *mf_triangular) : mf_wrapper(mf_triangular), mf_triangular(mf_triangular) {}
};

//' @title mf_trapezoidal_inf class
//' @name mf_trapezoidal_inf
//' @docType class
//' @description Class to build trapezoidal inf MF.
//'
//' @section Inherits:
//' mf_trapezoidal_inf class inherits all fields and methods of \link{mf} class.
//'
//' @section Methods:
//' \describe{
//' \item{\code{constructor(upper_kernel, upper_support)}}{argument: \code{upper_kernel} numerical upper value of kernel. \cr argument: \code{upper_support} numerical upper value of support. \cr return: an object of mf_trapezoidal_inf.}
//' }
//'
//' @examples
//' mf <- new(mf_trapezoidal_inf, 0, 1)
//' mf$degree(0.5)
//'
//' @export mf_trapezoidal_inf
class mf_trapezoidal_inf_wrapper : public mf_wrapper {

private:
	MFTRAPINF *mf_trapezoidal_inf;

public:
	mf_trapezoidal_inf_wrapper() : mf_trapezoidal_inf_wrapper(NULL) {
	    Rcpp::stop("mf_trapezoidal_inf default constructor not allowed");
	}
	mf_trapezoidal_inf_wrapper(double upper_kernel, double upper_support) : mf_trapezoidal_inf_wrapper(new MFTRAPINF(upper_kernel, upper_kernel, upper_support)) {}

	mf_trapezoidal_inf_wrapper(MFTRAPINF *mf_trapezoidal_inf) : mf_wrapper(mf_trapezoidal_inf), mf_trapezoidal_inf(mf_trapezoidal_inf) {}
};

//' @title mf_trapezoidal_sup class
//' @name mf_trapezoidal_sup
//' @docType class
//' @description Class to build trapezoidal sup MF.
//'
//' @section Inherits:
//' mf_trapezoidal_sup class inherits all fields and methods of \link{mf} class.
//'
//' @section Methods:
//' \describe{
//' \item{\code{constructor(lower_support, lower_kernel)}}{argument: \code{lower_support} numerical lower value of support. \cr argument: \code{lower_kernel} numerical lower value of kernel. \cr return: an object of mf_trapezoidal_sup.}
//' }
//'
//' @examples
//' mf <- new(mf_trapezoidal_sup, 0, 1)
//' mf$degree(0.5)
//'
//' @export mf_trapezoidal_sup
class mf_trapezoidal_sup_wrapper : public mf_wrapper {

private:
	MFTRAPSUP *mf_trapezoidal_sup;

public:
	mf_trapezoidal_sup_wrapper() : mf_trapezoidal_sup_wrapper(NULL) {
	    Rcpp::stop("mf_trapezoidal_sup default constructor not allowed");
	}

	mf_trapezoidal_sup_wrapper(double lower_support, double lower_kernel) : mf_trapezoidal_sup_wrapper(new MFTRAPSUP(lower_support, lower_kernel, lower_kernel)) {}

	mf_trapezoidal_sup_wrapper(MFTRAPSUP *mf_trapezoidal_sup) : mf_wrapper(mf_trapezoidal_sup), mf_trapezoidal_sup(mf_trapezoidal_sup) {}
};

//' @title mf_trapezoidal class
//' @name mf_trapezoidal
//' @docType class
//' @description Class to build trapezoidal MF.
//'
//' @section Inherits:
//' mf_trapezoidal class inherits all fields and methods of \link{mf} class.
//'
//' @section Methods:
//' \describe{
//' \item{\code{constructor(lower_support, lower_kernel, upper_kernel, upper_support)}}{argument: \code{lower_support} numerical lower value of support. \cr argument: \code{lower_kernel} numerical lower value of kernel. \cr argument: \code{upper_kernel} numerical upper value of kernel. \cr argument: \code{upper_support} numerical upper value of support. \cr return: an object of mf_trapezoidal.}
//' }
//'
//' @examples
//' mf <- new(mf_trapezoidal, 0, 1, 2, 3)
//' mf$degree(0.5)
//'
//' @export mf_trapezoidal
class mf_trapezoidal_wrapper : public mf_wrapper {

private:
	MFTRAP *mf_trapezoidal;

public:
	mf_trapezoidal_wrapper() : mf_trapezoidal_wrapper(NULL) {
	    Rcpp::stop("mf_trapezoidal default constructor not allowed");
	}

	mf_trapezoidal_wrapper(double lower_support, double lower_kernel, double upper_kernel, double upper_support) : mf_trapezoidal_wrapper(new MFTRAP(lower_support, lower_kernel, upper_kernel, upper_support)) {}

	mf_trapezoidal_wrapper(MFTRAP *mf_trapezoidal) : mf_wrapper(mf_trapezoidal), mf_trapezoidal(mf_trapezoidal) {}
};

RCPP_MODULE(FisPro) {
	using namespace Rcpp;

	class_<fis_wrapper>("fis")
	.default_constructor()
	.constructor<const char *>()
	.property("name", &fis_wrapper::get_name, &fis_wrapper::set_name)
	.property("input_size", &fis_wrapper::get_input_size)
	.property("output_size", &fis_wrapper::get_output_size)
	.method("infer", &fis_wrapper::infer)
	.method("infer_output", &fis_wrapper::infer_output);

	class_<mf_wrapper>("mf")
	.default_constructor()
	.property("name", &mf_wrapper::get_name, &mf_wrapper::set_name)
	.method("degree", &mf_wrapper::get_degree);

	class_<mf_triangular_wrapper>("mf_triangular")
	.derives<mf_wrapper>("mf")
	.default_constructor()
	.constructor<double, double, double>();

	class_<mf_trapezoidal_inf_wrapper>("mf_trapezoidal_inf")
	.derives<mf_wrapper>("mf")
	.default_constructor()
	.constructor<double, double>();

	class_<mf_trapezoidal_sup_wrapper>("mf_trapezoidal_sup")
	.derives<mf_wrapper>("mf")
	.default_constructor()
	.constructor<double, double>();

	class_<mf_trapezoidal_wrapper>("mf_trapezoidal")
	.derives<mf_wrapper>("mf")
	.default_constructor()
	.constructor<double, double, double, double>();
}
