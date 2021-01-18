//
//  PLS.h
//

#ifndef GenAlgPLS_PLS_h
#define GenAlgPLS_PLS_h

#include <memory>
#include "config.h"

#include <RcppArmadillo.h>

enum PLSMethod {
	SIMPLS = 0
};

class PLS {
protected:
	enum ViewState {
		UNKNOWN = 0,
		COLUMNS,
		ROWS
	};

public:
	PLS(const arma::mat &X, const arma::vec &Y) : X(X), Y(Y), currentViewState(UNKNOWN) {};
	virtual ~PLS() {};

	static std::unique_ptr<PLS> getInstance(PLSMethod method, const arma::mat &X, const arma::vec &Y);

	/**
	 * Reset the current view to be the original X and original Y matrix
	 * and then select the given columns from the X matrix
	 */
	virtual void viewSelectColumns(const arma::uvec &columns);

	/**
	 * Select the given rows from the current column-view for further
	 * processing
	 */
	virtual void viewSelectRows(const arma::uvec &rows);

	/**
	 * Select all rows from the current column-view for further
	 * processing
	 */
	virtual void viewSelectAllRows();


	/**
	 * Fit a PLS model to the data with the previously set view
	 * with up to ncomp components
	 * @param ncomp The maximum number of components to fit (0 ... fit as many as possible)
	 */
	virtual void fit(uint16_t ncomp = 0) =0;

	/**
	 * Get the coefficients of the last fit (i.e. coefficients
	 * that are obtained with ncomp specified in the last call
	 * to PLS::fit).
	 */
	virtual const arma::mat& getCoefficients() const = 0;

	/**
	 * Returns the intercept term for every number of components
	 * i.e. ncomp x nresp matrix
	 */
	virtual const arma::vec& getIntercepts() const = 0;

	/**
	 * Get dimensions of original matrices Y and X
	 */
	arma::uword getNumberOfPredictorVariables() const { return this->X.n_cols; }
	arma::uword getNumberOfResponseVariables() const { return this->Y.n_cols; }
	arma::uword getNumberOfObservations() const { return this->X.n_rows; }

	/**
	 * Returns the number components the last fit was performed with
	 */
	uint16_t getResultNComp() const { return this->resultNComp; }

	/**
	 * Predict the values with all or `ncomp` components
	 * @param ncomp The 1-based number of components
	 */
	arma::vec predict(const arma::mat &newX, uint16_t ncomp) const;
	arma::mat predict(const arma::mat &newX) const;

	const arma::mat & getXColumnView() const { return this->viewXCol; }
	const arma::vec & getY() const { return this->Y; }

	virtual std::unique_ptr<PLS> clone() const = 0;

protected:
	const arma::mat X;
	const arma::vec Y;

	uint16_t resultNComp;

	ViewState currentViewState;
	arma::mat viewXCol;
	arma::vec viewY;
	arma::mat viewX;
};

#include "PLSSimpls.h"

#endif
