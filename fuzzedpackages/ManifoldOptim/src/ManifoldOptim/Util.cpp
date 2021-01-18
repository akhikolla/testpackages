#include "Util.h"

void CopyFrom(SmartSpace* s, const arma::mat& x)
{
	size_t n = x.n_rows;
	size_t p = x.n_cols;

	double* sptr = s->ObtainWriteEntireData();
	for (size_t i = 0; i < n; ++i) {
		for (integer j = 0; j < p; ++j) {
			sptr[COLM(i, j, n)] = x(i,j);
		}
	}
}

arma::vec ToArmaVec(const SmartSpace* s)
{
	size_t nn = s->Getlength();
	const double* sptr = s->ObtainReadData();

	arma::vec x(nn);
	for (size_t i = 0; i < nn; ++i) {
		x(i) = sptr[i];
	}

	return x;
}

arma::mat ToArmaMat(const SmartSpace* s)
{
	if (s->Getls() == 1) {
		throw ManifoldOptimException("Expect Element to have exactly two dimensions (it has one)");
	} else if (s->Getls() > 2 && s->Getsize()[2] > 1) {
		throw ManifoldOptimException("Expect Element to have exactly two dimensions (has a non-trival third dimension)");
	}

	integer n = s->Getsize()[0];
	integer p = s->Getsize()[1];
	const double* sptr = s->ObtainReadData();

	arma::mat x(n, p);
	for (integer i = 0; i < n; ++i) {
		for (integer j = 0; j < p; ++j) {
			x(i,j) = sptr[COLM(i, j, n)];
		}
	}

	return x;
}

arma::mat ToArmaMat(const ProductElement* s)
{
	integer nElem = s->GetNumofElement();
	if (nElem != 1) {
		throw ManifoldOptimException("Expect ProductElement to have exactly one element");
	}
	Rprintf("ProductElement had exactly one element!");
	return ToArmaMat(s->GetElement(0));
}

void CopyFrom(SmartSpace* s, const NumericMatrix& x)
{
	size_t n = x.nrow();
	size_t p = x.ncol();

	double* sptr = s->ObtainWriteEntireData();
	for (size_t i = 0; i < n; ++i) {
		for (integer j = 0; j < p; ++j) {
			sptr[COLM(i, j, n)] = x(i,j);
		}
	}
}

void CopyFrom(NumericMatrix& x, const SmartSpace* s)
{
	size_t n = x.nrow();
	size_t p = x.ncol();

	const double* sptr = s->ObtainReadData();
	for (size_t i = 0; i < n; ++i) {
		for (integer j = 0; j < p; ++j) {
			x(i,j) = sptr[COLM(i, j, n)];
		}
	}
}

void CopyFrom(SmartSpace* s, const NumericVector& x)
{
	size_t n = x.size();
	double* sptr = s->ObtainWriteEntireData();
	for (size_t i = 0; i < n; ++i) {
		sptr[i] = x(i);
	}
}

NumericVector ToNumericVector(const SmartSpace* s)
{
	integer n = s->Getlength();
	const double* sptr = s->ObtainReadData();

	NumericVector x(n);
	for (integer i = 0; i < n; ++i) {
		x(i) = sptr[i];
	}

	return x;
}

List ExtractElements(const SmartSpace* s)
{
	integer nn = s->Getlength();
	integer ndim = s->Getls();
	const integer* sizes = s->Getsize();
	const double* sptr = s->ObtainReadData();

	Rcpp::IntegerVector ds(ndim);
	for (integer i = 0; i < ndim; ++i) {
		ds(i) = sizes[i];
	}

	Rcpp::Dimension dim(ds);
	Rcpp::NumericVector arr(dim);
	
	for (integer i = 0; i < nn; ++i) {
		arr(i) = sptr[i];
	}

	return Rcpp::List::create(arr);
}

List ExtractElements(const ProductElement* s)
{
	integer nElem = s->GetNumofElement();
	List ret(nElem);

	for (int i = 0; i < nElem; i++) {
		const Element* el = s->GetElement(i);
		const List& singletonList = ExtractElements(el);
		ret[i] = singletonList[0];
	}

	return ret;
}

