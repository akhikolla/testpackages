#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::cube csmooth_sam(arma::cube myCube, int wind, int bins, int nslices){
	// wind = window size.
	// bins = how many clusters we had/have?
	// nslices = how many slices, starting from zero, should be smoothered.
	int max_drow = myCube.n_rows;
  int max_dcol = myCube.n_cols;
	
	arma::cube tmp_Cube(size(myCube.slices(0, nslices)));
	
	for(int dslice = 0; dslice <= nslices; ++dslice)
	{
	  arma::mat myCubeSlice = myCube.slice(dslice);
	  
		for(int drow = 0; drow < max_drow; ++drow)
		{
			for(int dcol = 0; dcol < max_dcol; ++dcol)
			{
				int first_row = (drow - wind > 0) ? drow - wind : 0;
				int first_col = (dcol - wind > 0) ? dcol - wind : 0;
				int last_row = (drow + wind < max_drow - 1) ? drow + wind : max_drow -1;
				int last_col = (dcol + wind < max_dcol - 1) ? dcol + wind : max_dcol -1;

				arma::vec electors = arma::vectorise(myCubeSlice.submat(first_row, first_col, last_row, last_col));

				// find the mode
				// we can do that using an histogram
				std::vector<int> histogram(bins+1, 0);
				for(arma::uword i = 0; i < electors.n_elem; ++i)
					++histogram[ electors[i] ];
				int vote = std::max_element( histogram.begin(), histogram.end() ) - histogram.begin();

				tmp_Cube(drow, dcol, dslice) = vote;
			}
		}
	}
	return tmp_Cube;
}

