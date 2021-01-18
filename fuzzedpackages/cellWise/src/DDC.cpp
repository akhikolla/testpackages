
#include "DDC.h"
#include "ANN.h"

arma::uvec DDC::vinter(const arma::uvec &first, const arma::uvec &second) {
  //returns uvec with elements from first which are also in second
  std::vector<arma::uword > output;
  std::set_intersection(first.begin(), first.end(), second.begin(), second.end(),
                        std::back_inserter(output));
  arma::uvec result = arma::conv_to< arma::uvec >::from(output);
  return result;
}

arma::uvec DDC::vdiff(const arma::uvec &first, const arma::uvec &second)
{ //returns uvec with elements from first which are not in second
  std::vector<arma::uword > output;
  std::set_difference(first.begin(), first.end(), second.begin(), second.end(),
                      std::back_inserter(output));
  
  arma::uvec result = arma::conv_to< arma::uvec >::from(output);
  return result;
}


arma::uvec DDC::col2cell(const arma::uvec &colNrs, const int n) {
  // Transforms column indices to cellwise indices.
  // Here colNrs is a vector with column numbers between 0 and d-1.
  arma::umat cindex(n, colNrs.size());
  // cindex.each_row() = ((colNrs - 1) * n).t();
  cindex.each_row() = ((colNrs) * n).t();
  cindex.each_col() += arma::regspace<arma::uvec>(0, (n - 1));
  return(arma::vectorise(cindex));
}

arma::uvec DDC::row2cell(const arma::uvec &rowNrs, const int n, const int d) {
  // Transforms row indices to cellwise indices.
  // Here rowNrs is a vector with row numbers between 0 and n-1.
  
  arma::umat cindex(d, rowNrs.size(),arma::fill::zeros);
  cindex.each_row() = (rowNrs).t();
  cindex.each_col() += arma::regspace<arma::uvec>(0, n, n * (d - 1));
  return(arma::vectorise(cindex));
}

double DDC::weightedMedian(arma::vec x, arma::vec weights) {
  double wmed = 0;
  // throw out m	issing/inf elements or elements with missing/inf weights:
  // 
  
  x = x(arma::find_finite(weights));
  weights = weights(arma::find_finite(weights));
  weights = weights(arma::find_finite(x));
  x = x(arma::find_finite(x));
  x = x(arma::find(weights > 0));
  weights = weights(arma::find(weights > 0));
  
  if (x.size() == 0) {
    return(wmed);
  }
  if (x.size() == 1) {
    return(x(0));
  }
  
  double wtotal = sum(weights);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* All weights equal?  Happens if +Inf were detected.                  */
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  arma::uvec uniquew = arma::find_unique(weights);
  if (uniquew.size() == 1) {
    wmed = arma::median(x);
    return (wmed);
  }
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* Sort x and calculated the cumulative sum of weights (normalize to   */
  /* one) according to the reordered vector.                             */
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* (a) Sort x */
  arma::uvec xorder = arma::sort_index(x);
  x = x(xorder);
  weights = weights(xorder);
  
  /* (b) Normalized cumulative weights */
  arma::vec wcum = arma::cumsum(weights);
  double tmp_d = 0;
  double tmp_d2 = 0;
  /* Index where cumulative weight passed 1/2 */
  unsigned int half = x.size() + 1; /* Default is last */
  
  /* Adjust */
  for (unsigned int ii = 0; ii < x.size(); ii++) {
    tmp_d = weights(xorder(ii)) / wtotal;
    tmp_d2 = tmp_d2 + tmp_d;
    wcum(ii) = tmp_d2 - (tmp_d / 2);
    if (wcum(ii) >= 0.5) {
      half = ii;
      /* Early stopping - no need to continue */
      break;
    }
  }
  
  /* Two special cases where more than half of the total weight is at
   a) the first, or b) the last value */
  if (half == 0 || half == x.size()) {
    wmed = x(half);
    return wmed;
  }
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* Alt 1: Linearly interpolated weighted median                        */
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  /* The width and the height of the "rectangle". */
  double dx = x(half) - x(half - 1);
  double Dy = wcum(half) - wcum(half - 1);
  
  /* The width and the height of the triangle which upper corner touches
   the level where the cumulative sum of weights *equals* half the
   total weight. */
  double dy = 0.5 - wcum(half);
  dx = (dy / Dy) * dx;
  /*    printf("dx=%g, dy=%g\n", dx, dy); */
  
  /* The corresponding x value */
  wmed = dx + x(half);
  
  return wmed;
  
}


arma::vec DDC::predictCol(const arma::vec &colj, const arma::mat &U, const int coln,
                          const arma::umat &ngbrs, const arma::mat &corrweight,
                          const arma::mat &robslopes, const int combinRule) {
  // Predicts the values in column colj using the set 'ngbrs' of
  // columns of U, by applying the combination rule 'combinRule' whose
  // inputs are the weights in 'corrweight' and the slopes in 'robslopes'.
  // Assumes that the first entry of colj is the number of the column.
  // Assumes the remainder of colj is a vector with same length as the
  // remainder of the columns of U, and that all of these columns are
  // already centered.
  
  arma::uvec contributors = arma::find(corrweight.row(coln) > 0);
  arma::vec estcol(colj.size(), arma::fill::zeros);
  
  if (contributors.size() > 1) {
    arma::uvec ngb1 = ngbrs.row(coln).t();
    ngb1 = ngb1(contributors);
    arma::vec slopes1 = robslopes.row(coln).t();
    slopes1 = slopes1(contributors);
    arma::vec corrwt1 = corrweight.row(coln).t();
    corrwt1 = corrwt1(contributors);
    arma::mat ZestAllh = U.cols(ngb1); // is n by k matrix
    // Predicts column j from each column h, using slope(Z_j ~Z_h).
    // This array has the estimates from k variables.
    ZestAllh.each_row() %= slopes1.t();
    
    switch(combinRule) {
    case 1:
      for (unsigned int i = 0; i < colj.size(); i++) {
        arma::vec dotprod = ZestAllh.row(i).t() % corrwt1;
        estcol(i) = sum(dotprod(arma::find_finite(dotprod))) / 
          sum(corrwt1(arma::find_finite(dotprod))); //weighted mean
      }
      break;
    case 2:
      for (unsigned int i = 0; i < colj.size(); i++) {
        arma::vec temprow = ZestAllh.row(i).t();
        estcol(i) = weightedMedian(temprow, corrwt1); //weighted median
      }
      break;
    case 3:
      for (unsigned int i = 0; i < colj.size(); i++) {
        arma::vec temprow = ZestAllh.row(i).t();
        estcol(i) = mean(temprow(find_finite(temprow))); //mean
      }
      break;
    case 4:
      for (unsigned int i = 0; i < colj.size(); i++) {
        arma::vec temprow = ZestAllh.row(i).t();
        estcol(i) = median(temprow(arma::find_finite(temprow))); //median
      }
      break;
    }
  }
  return(estcol);
}

double DDC::slopeMedWLS(const arma::vec &xcol, const arma::vec &colj,
                        double qRegr, double precScale) {
  // Computes the slope of a robust regression without intercept
  // of the column colj on the column xcol.
  // The computation starts by Median regression, followed by
  // weighted least squares(WLS) in which the weights are
  // determined by the quantile qRegr.
  // Assumes that colj is a vector with the same length as xcol
  // and that both columns are already centered.
  
  arma::vec ratio = colj / xcol;
  double slope = 0;
  arma::uvec ratiofin = find_finite(ratio);
  if (ratiofin.size() <= 3) {
    slope = 0;
  }
  else {
    double rawb = median(ratio(find_finite(ratio))); // raw slope
    if (!std::isfinite(rawb)) {
      slope = 0;
    }
    else {
      rawb = rawb > 2 ? 2 : rawb;
      rawb = rawb < -2 ? -2 : rawb;
      // Now compute weighted LS slope :
      arma::vec r = colj - rawb * xcol; // raw residuals
      double cutoff = qRegr * LocScaleEstimators::scale1StepM(r, LocScaleEstimators::rhoHuber25,
                                                              arma::datum::nan, precScale);
      // cutoff can be zero, which is okay.
      arma::uvec rowSel = find(arma::abs(r) <= cutoff); // selects the inliers
      rowSel = rowSel(find_finite(rowSel));
      if (rowSel.size() > 0) {
        arma::vec xw = xcol(rowSel);
        arma::vec yw = colj(rowSel);
        slope = dot(xw, yw) / pow(norm(xw), 2); // slope of colj ~xcol + 0
        if (!std::isfinite(slope)) { slope = 0; }
      }
      else {
        slope = 0;
      }
    }
  }
  return(slope);
}

arma::vec DDC::compSlopes(const arma::vec &colj, arma::uvec ngbrs, const arma::mat &U,
                          double qRegr, double precScale) {
  // For a given column colj this computes the slopes(obtained by
  // robSlope) of regressing colj on each of k given columns of
  // the matrix U.
  // Assumes that ngbrs are the indices of
  // those k columns.
  // Assumes the remainder of colj is a vector with same length as the
  // remainder of the columns of U, and that all of these columns are
  // already centered.
  const int k = ngbrs.size();
  arma::vec slopes(k, arma::fill::zeros); // initialize
  ngbrs = ngbrs(arma::find_finite(ngbrs)); // keep only the necessary neighbors
  ngbrs = ngbrs(arma::find(ngbrs < U.n_cols)); // keep only neighbors with inbound indices
  if (ngbrs.size() > 0) {
    arma::vec b(ngbrs.size(), arma::fill::zeros);
    for (unsigned int i = 0; i < ngbrs.size(); i++) {
      b(i) = slopeMedWLS(U.col(ngbrs(i)), colj, qRegr, precScale);
    }
    slopes.head(ngbrs.size()) = b;
  }
  return(slopes);
}

double DDC::corrGKWLS(arma::vec xcol, double qCorr, arma::vec colj, double precScale) {
  // Computes a robust correlation between the columns xcol and
  // colj using the formula of Gnanadesikan - Kettenring(GK),
  // followed by a Weighted Pearson correlation.
  // qCorr is a quantile used to determine the weights.
  // Assumes colj is a vector with same length as xcol
  // and that all normalizations have already happened.
  double corr = 0;
  arma::vec sumvec = xcol + colj;
  arma::uvec sumvecfin = find_finite(sumvec);
  if (sumvecfin.size() <= 3) {
    //we need at least 3 valids
    corr = 0;
  }
  else {
    arma::vec difvec = xcol - colj;
    corr = (std::pow(LocScaleEstimators::scale1StepM(sumvec, LocScaleEstimators::rhoHuber25,
                                                     arma::datum::nan, precScale), 2) - 
                                                       std::pow(LocScaleEstimators::scale1StepM(difvec, LocScaleEstimators::rhoHuber25,
                                                                                                arma::datum::nan, precScale), 2)) / 4;
    // This corr should not be NA since data columns with too many NAs
    // have already been taken out.But just in case someone increases
    // the allowed fraction fracNA too much, we add the precaution :
    if (!std::isfinite(corr)) { corr = 0; }
    else {
      corr = corr > 0.99 ? 0.99 : corr;
      corr = corr < -0.99 ? -0.99 : corr;
      
      
      // Compute reweighted Pearson correlation :
      arma::mat corrMatInv(2, 2, arma::fill::eye);
      
      corrMatInv(1) = -corr; // replace off - diagonal elements
      corrMatInv(2) = -corr;
      
      corrMatInv = corrMatInv / std::abs(1 - std::pow(corr, 2)); // divide by determinant
      arma::mat xtemp = join_rows(xcol, colj);
      arma::vec RDi2 = sum((xtemp * corrMatInv) % xtemp, 1);
      arma::uvec rowSel = find(RDi2 < qCorr);
      rowSel = rowSel(arma::find_finite(rowSel));
      if (rowSel.size()>0) {
        xcol = xcol(rowSel);
        colj = colj(rowSel);
        corr = arma::norm_dot(xcol, colj);
      }
      else {
        corr = 0;
      }
    }
  }
  return(corr);
}

arma::vec DDC::limitFilt(arma::vec v, double qCut) {
  // Detects outliers and sets them to NA.
  // Assumes that the data have already been standardized.
  v.for_each([qCut](arma::mat::elem_type &value) {
    value = std::abs(value) > qCut ? arma::datum::nan : value;
  });
  return v;
}

arma::vec DDC::rawEquiGYfilt(const arma::vec &v, double qCut) {
  // raw version, has no iteration
  // assumes qCut is quantile of chi square distribution
  arma::vec v2 = arma::pow(v, 2);
  arma::vec u = arma::sort(v2);
  arma::uvec indices = arma::find(u < qCut);
  arma::uvec i0 = arma::find(u < qCut); // initial inliers
  double n0 = 0;
  
  if (i0.size() > 0) {
    double dn = max(arma::erf(sqrt((u.tail((u.size() - i0.size())) / 2))) -
                    (arma::regspace(i0.size(), (v.size() - 1)) - 1) / v.size());
    dn = dn > 0 ? dn : 0;
    n0 = round(dn * v.size());
  }
  arma::vec vna = v;
  if (n0 > 0) {
    double cutoff = u[(v.size() - n0 + 1)];
    vna(arma::find(v2 >= cutoff)).fill(arma::datum::nan);
  }
  return vna;
}

arma::vec DDC::equiGYfilt(const arma::vec &v, double qCut, const int miter) {
  // Detects outliers and sets them to NA.This is a
  // permutation - equivariant version of the GY filter.
  // Assumes that the data have already been standardized.
  // assumes qCut is sqrt of a chi square quantile
  qCut = pow(qCut,2);
  
  int converge = 0;
  int	iter = 0;
  arma::uvec observed = find_finite(v);
  arma::vec	vobs = v(observed);
  arma::uvec id = arma::regspace<arma::uvec>(0,(vobs.size()-1));
  while ((converge == 0) && (iter < miter)) {
    iter = iter + 1;
    vobs = rawEquiGYfilt(vobs, qCut);
    id = id(arma::find_finite(vobs));
    if (vobs.is_finite() ){
      converge = 1;
    }
    vobs = vobs(arma::find_finite(vobs));
  }
  arma::vec vobsout(observed.size());
  vobsout.fill(arma::datum::nan);
  vobsout(id) = vobs;
  arma:: vec vout(v.size());
  vout.fill(arma::datum::nan);
  vout(observed) = vobsout;
  return(vout);
}

DDC::kbestcorr DDC::kBestCorr(const arma::vec &colj, const arma::mat &U, const int coln,
                              const unsigned int k, double qCorr, double precScale) {
  // For a given column colj this computes the k highest absolute
  // correlations(obtained by robCorr) between colj and the columns
  // of the matrix U, and also reports which k columns are selected.
  // Assumes that coln gives the number of the column
  // Assumes colj is a vector with same length as the
  //  the columns of U, and that all normalizations have
  // already happened.
  // Assumes k <= d - 1 which is ensured before calling this function.
  
  kbestcorr result = { arma::zeros<arma::uvec>(k),arma::zeros<arma::vec>(k) };
  arma::vec allCorrs(U.n_cols, arma::fill::zeros);
  for (unsigned int i = 0; i < U.n_cols; i++) {
    allCorrs(i) = corrGKWLS(U.col(i), qCorr, colj, precScale);
  }
  arma::uvec selected = arma::sort_index(arma::abs(allCorrs), "descend");
  selected = selected(arma::find(selected != coln)); // removes correlation of j with j
  // selected has length d - 1
  selected = selected.head(k);  // always works since k <= d - 1
  result.corrs = allCorrs(selected);
  result.selected = selected;
  return result;
}

double DDC::deShrink(const arma::vec &colj, const arma::mat &Z,
                     const int coln, double qRegr, double precScale) {
  // Deshrinks the column colj by multiplying it by the robSlope
  // of column j of the matrix Z on colj.
  // Assumes that the first entry of colj is the number of the column.
  // Assumes the remainder of colj is a vector with same length as the
  // columns of Z, and that both columns were already centered.
  
  arma::vec zj = Z.col(coln);   // column with response(from Z)
  double a = slopeMedWLS(colj, zj, qRegr, precScale);
  // arma::vec deshrunk = a * colj;
  return(a);
}


// Code for Fastrobcor

void DDC::get_NN_2Set(double *data, double *query, int *D, int *ND, int *NQ,
                      int *K, double *EPS, int *SEARCHTYPE, int *USEBDTREE,
                      double *SQRAD, int *nn_index, double *distances) {
  const int d = *D;		// Number of Dimensions for points
  const int nd = *ND;		// Number of Data points
  const int nq= *NQ;		// Number of Query points
  const int k = *K;		// Maximum number of Nearest Neighbours
  
  const int searchtype = *SEARCHTYPE;
  const bool usebdtree = *USEBDTREE?true:false;
  
  const double error_bound = *EPS;	// enough said!
  const double sqRad = *SQRAD;		// Squared Radius for rad search
  
  ANNkd_tree	*the_tree;	// Search structure
  
  ANNpointArray data_pts 	= annAllocPts(nd,d);		// Allocate data points
  ANNidxArray nn_idx 		= new ANNidx[k];		// Allocate near neigh indices
  ANNdistArray dists 		= new ANNdist[k];		// Allocate near neighbor dists
  
  int *d_ptr = new int[d];
  int ptr = 0;
  
  // set up column offsets for query point matrix (to convert Row/Col major)
  for(int i = 0; i < d; i++)
  {
    d_ptr[i] = i*nd;
  }
  
  for(int i = 0; i < nd; i++) // now construct the points
  {
    for(int j = 0; j < d; j++)
    {
      data_pts[i][j]=data[ d_ptr[j]++ ];
    }
  }
  
  if(usebdtree){
    the_tree = new ANNbd_tree(	// Build search structure
      data_pts,			// The data points
      nd,					// Number of data points
      d);					// Dimension of space				
  } else {
    the_tree = new ANNkd_tree( data_pts, nd, d);
  }
  
  // set up offsets for query point matrix (to convert Row / Col major)
  for(int i = 0; i < d; i++)
  {
    d_ptr[i] = i*nq;
  }
  
  ANNpoint pq = annAllocPt(d);
  for(int i = 0; i < nq; i++)	// Run all query points against tree
  {
    // read coords of current query point
    for(int j = 0; j < d; j++)
    {
      pq[j]=query[ d_ptr[j]++ ];
    }
    
    switch(searchtype){
    case 1:
      the_tree->annkSearch(	// search
          pq,	// query point
          k,		// number of near neighbors
          nn_idx,		// nearest neighbors (returned)
          dists,		// distance (returned)
          error_bound);	// error bound			
      break;
      
    case 2:  // Priority search
      the_tree->annkPriSearch(pq, k, nn_idx, dists, error_bound);
      break;
      
    case 3: // Fixed radius search 
      the_tree->annkFRSearch(	pq,	sqRad, k, nn_idx, dists,error_bound);			
      break;
    }		
    
    for (int j = 0; j < k; j++)
    {
      distances[ptr] = ANN_ROOT(dists[j]);	// unsquare distance
      nn_index[ptr++]  = nn_idx[j];	// put indices in returned array
      //return C/C++ indices (not converted to R-indices with +1)
    }
  }
  // Do a little bit of memory management......
  annDeallocPt(pq);
  annDeallocPts(data_pts);
  delete [] nn_idx;
  delete [] dists;
  delete [] d_ptr;
  delete the_tree;
}

arma::vec DDC::transClassic(arma::vec y, const double precScale) {
  
  arma::uvec finiteinds = arma::find_finite(y);
  double yloc = mean(y(finiteinds));
  y = y - yloc; // still has NAs, will be zeroed later
  double rmsv = sqrt(mean(y(finiteinds) % y(finiteinds)));
  // rmsv
  y = y / rmsv;
  return(y);
}




DDC::fastRobCorout DDC::FastRobCorActual(const arma::mat &X, const arma::vec &locX,
                                         const arma::vec &scaleX,
                                         unsigned int k, int qdim,
                                         const unsigned int nCorr, int absCorr,
                                         int transFun,  double precScale,
                                         int treetype, int searchtype,  double radius,
                                         double eps, int includeSelf) {
  // Finds k(approximate) nearest neighbors of each column of X
  // in the sense of(absolute value of) robust correlation.
  // Assumes the data is such that checkDataSet(X) equals X.
  // X is the input data, and must be a matrix or a data frame.
  //   It must always be provided.We denote n = nrow(X), d = ncol(X).
  // k is the number of nearest neighbors sought, excluding the
  //   variable itself, so k is between 1 and d - 1.
  // absCorr: F means look for the k highest correlations,
  //          T means look for the k highest absolute correlations.
  // qdim 0 means do not apply dimension reduction.Otherwise qdim
  //      is the number of dimensions to which the transformed data is
  //      reduced before applying nn2.q lies between 1 and n, but
  //      should not be much higher than 20 (? ) for nn2 to work well
  //      when an apprimate algorithm is used(i.e.eps = 0).
  // robLoc   Function for the robust location of each column.
  // robScale Function for the robust scale of each column.
  // transFun Function tansforming each standardized column, which
  //          determines the type of robust correlation used, e.g.
  //          transRank yields Spearman rank correlation.
  //  treetype    standard 'kd' tree, or a 'bd' (box - decomposition, AMNSW98)
  //             tree which may perform better for larger point sets.
  // searchtype 'priority' visits cells in increasing order of distance
  //             from the query point, and hence, should converge more
  //             rapidly on the true nearest neighbor, but 'standard' is
  //             usually faster for exact searches.
  //       'radius' only searches for neighbors within a specified
  //             radius of the point, so will obtain AT MOST k neighbors.
  //             If there are no neighbors then nn.idx will contain 0
  //             and nn.dists will contain 1.340781e+154 for that point.
  // radius	    Radius of search for searchtype = 'radius'
  // eps	        Error bound : default of 0.0 implies exact nearest
  //             neighbour search.
  
  // First, we include all the functions that we need here and will
  // not use by themselves.
  
  
  // Turn data into a matrix if possible
  
  fastRobCorout result;
  arma::imat ngbrs ;
  arma::mat robcors;
  arma::mat Ufull = X;
  arma::mat U;
  
  
  ///////////////////////////////////////////////////////
  //    STEP 1: TRANSFORM COLUMNS TO REMOVE OUTLIERS   //
  ///////////////////////////////////////////////////////
  
  
  switch(transFun) {
  case 1: //transHuber
    for (unsigned int i = 0; i < Ufull.n_cols; i++) {
      arma::vec y = Ufull.col(i);
      arma::uvec finiteinds = arma::find_finite(y);
      arma::vec u = y(finiteinds) - locX(i);
      double ycap = 1.5 * scaleX(i);
      u.transform( [ycap](double val) { 
        if(std::abs(val) > ycap) {
          return(ycap * ((val > 0) - (val < 0)));
        }else{
          return(val);
        }
      });
      y(finiteinds) = u;
      Ufull.col(i) = y;
    }
    break;
  case 2: //transWrap
    for (unsigned int i = 0; i < Ufull.n_cols; i++) {
      arma::uvec finiteinds = arma::find_finite(Ufull.col(i));
      arma::vec u = Ufull.col(i) - locX(i);
      u = u / scaleX(i);
      arma::vec ufin = u(finiteinds);
      LocScaleEstimators::psiTanh(ufin);
      u(finiteinds) = ufin * scaleX(i);
      Ufull.col(i) = u;
    }
    break;
  case 3: // transRank
    for (unsigned int i = 0; i < Ufull.n_cols; i++) {
      arma::uvec finiteinds = arma::find_finite(Ufull.col(i));
      arma::vec u = Ufull.col(i);
      arma::vec ufin = u(finiteinds);
      u(finiteinds) = LocScaleEstimators::rank(ufin);
      u(finiteinds) = u(finiteinds) - arma::mean(u(finiteinds));
      Ufull.col(i) = u;
    }
    break;
  }
  
  /////////////////////////////////////////////////
  //    STEP 0: SELECT rows TO CALCULATE NGBRS   //
  /////////////////////////////////////////////////
  
  // 
  
  arma::uvec rowNAs(Ufull.n_rows, arma::fill::zeros);
  arma::uvec selectedInds;
  
  for (unsigned int i = 0; i < Ufull.n_rows; i++) {
    arma::uvec nonFiniteElem = arma::find_nonfinite(Ufull.row(i).t());
    rowNAs(i) = nonFiniteElem.size();
  }
 
  if ((unsigned int) std::count(rowNAs.begin(), rowNAs.end(), 0) >= nCorr) {
    // in this case we have enough NA-free rows
    // select nCorr rows randomly 
    arma::uvec indices = arma::find(rowNAs == 0);
    
    if (indices.size() == nCorr) {
      U = Ufull;
    } else {
      selectedInds = LocScaleEstimators::sample(indices, nCorr, false);
      U = Ufull.rows(selectedInds);
    }
  } else { // select nCorr rows with lowest number of NAs
    
    // Determine max number of missings in the nCorr rows with lowest nb of NAs
    selectedInds = arma::sort_index(rowNAs);
    unsigned int nbMisCut = rowNAs(selectedInds(nCorr - 1)); 
    
    // now decide how many rows have exactly nbMiscut missings
    // those are the ones we still have to sample to get nCorr rows in total
    
    arma::uvec sortNAs = rowNAs(selectedInds); // sorted nb of NAs, increasing
    auto it_lo = std::lower_bound(sortNAs.begin(),
                                  sortNAs.end(), nbMisCut); //iterator pointing to first element in range greater or equal to nbMisCut
    unsigned int nFixed = std::distance(sortNAs.begin(), it_lo); // index of that first element == number of 'fixed' rows, i.e. rows with #NAs < nbMisCut
    auto it_hi = std::upper_bound(sortNAs.begin(),
                                  sortNAs.end(), nbMisCut);//iterator pointing to first element in range GREATER than nbMisCut
    unsigned int nVariable = std::distance(sortNAs.begin(), it_hi);// first element with > nbMisCut NAs, which is equal to nb of rows with #NAs <= nbMiscut
    nVariable = nVariable - nFixed;

    // Now take all the rows with < nCorr missings, and sample the remaining
    // from the rows with precisely nCorr missings
    arma::uvec fixedSelection = selectedInds.head(nFixed);
    arma::uvec sampleFrom = selectedInds.subvec(nFixed, nFixed + nVariable - 1);
    arma::uvec variableSelection =  LocScaleEstimators::sample(sampleFrom, nCorr - nFixed, false);
    selectedInds = arma::join_cols(fixedSelection, variableSelection);
    
    U = Ufull.rows(selectedInds); // now we have selected nCorr rows
    
    // generate extra NAs so that each VARIABLE has same number of NAs!
    // this is needed to make the ANN procedure work after replacing
    // the NAs by zeroes, assuming the NAs are MCAR
    
    arma::uvec colFins(U.n_cols, arma::fill::zeros); // # finite elements in columns of U
    for (unsigned int i = 0; i < U.n_cols; i++) {
      arma::uvec finiteElem = arma::find_finite(U.col(i));
      colFins(i) = finiteElem.size();
    }

    
    unsigned int nbFinCut = arma::min(colFins); // min number of finite elements in column of our reduced dataset U
    // here we need a safety check: if nbFincut < 0.5 * nCorr, we consider the columns as standAlone and recalculate nbFincut
    if (nbFinCut < 0.5 * nCorr) {
      arma::uvec standAlonecols = arma::find(colFins < 0.5 * nCorr);
   
      if (standAlonecols.size() < colFins.size()) {
        nbFinCut = arma::min(colFins(DDC::vdiff(arma::regspace<arma::uvec>(0, colFins.size() - 1),
                                                standAlonecols)));
      }
    }
    
    for (unsigned int i = 0; i < U.n_cols; i++) {
      arma::vec colTemp = U.col(i);
      arma::uvec finiteElem = arma::find_finite(colTemp);
      if (finiteElem.size() > nbFinCut) { // generate extra NAs, not that it doesn't do so for standAlone cols
        arma::uvec makeNA = LocScaleEstimators::sample(finiteElem, finiteElem.size() - nbFinCut, false);
        colTemp(makeNA).fill(arma::datum::nan);
        U.col(i) = colTemp;
      }
    }

  }
  
  
  /////////////////////////////////
  //    STEP 2: STANDARDIZE U    //
  /////////////////////////////////
  
  // Set NA's to zero in u. This is needed for the relation between
  // correlation and euclidean distance.As a side effect, variables
  // with many NA's are less likely to be selected as neighbors.
  // it makes sense since the transformations above leave U with centered
  // variables
  
  U(arma::find_nonfinite(U)).zeros();
  // round(U, 2)
  
  // standardize the columns of U.This must be done the classical way.
  for (unsigned int i = 0; i < U.n_cols; i++) {
    U.col(i) = transClassic(U.col(i), precScale);
  }
  
  arma::mat tU = U.t(); // has d rows and n columns
  
  
  
  ///////////////////////////////////////
  //    STEP 3: DIMENSION REDUCTION    //
  ///////////////////////////////////////   
  
  // tU has d rows and n columns
  
  int lowdim = U.n_cols < U.n_rows ? U.n_cols : U.n_rows;
  
  if (qdim != 0) {
    lowdim = lowdim < qdim ? lowdim : qdim; 
  }
  
  if ( (unsigned int) lowdim < U.n_cols) {
    // reducing the dimension from n to lowdim :
    // pca through svd: can't use the built-in pca,
    // because we don't need centering
    arma::vec s; //holding singular values
    arma::mat Usvd, L, scores;
    
    arma::svd_econ(Usvd, s, L, tU, "both","dc");
    scores = tU * L;
    
    tU = scores.head_cols(lowdim); // now tU is d*lowdim
  }
  ///////////////////////////////////////////////////////////
  //    STEP 4: FIND k NEAREST NEIGHBORS OF EACH COLUMN    //
  ///////////////////////////////////////////////////////////  
  
  // Find approximate k nearest neighbors of each variable
  // we look for 2*k +1 nearest neighbors, and then select the k + 1 closest ones
  // 
  
  int ktemp = (int) ((2 * k + 1) < X.n_cols ? 2 * k + 1 : X.n_cols); // max number of nearest neighbours
  
  
  
  if (absCorr == 1) {
    arma::mat data = arma::join_cols(tU, -tU); // to deal with negative correlations
    int dtu = tU.n_cols;
    int ntu = tU.n_rows;
    int n2tu = data.n_rows;
    
    arma::ivec nn_index(X.n_cols * ktemp);
    arma::vec distances(X.n_cols * ktemp);
    get_NN_2Set(data.memptr(), tU.memptr(), &dtu, &n2tu,
                &ntu, &ktemp, &eps, &searchtype, &treetype, &radius,
                &nn_index[0], &distances[0]);
                
                ngbrs = arma::imat(&nn_index[0], ktemp, X.n_cols,true,false).t();
                ngbrs = ngbrs.tail_cols(ktemp - 1); // exclude self
                
                //correct row
                ngbrs.elem(arma::find(ngbrs > (X.n_cols - 1))) -= X.n_cols;
  } else {
    int dtu = tU.n_cols;
    int ntu = tU.n_rows;
    arma::ivec nn_index(X.n_cols * ktemp);
    arma::vec distances(X.n_cols * ktemp);
    
    get_NN_2Set(tU.memptr(), tU.memptr(), &dtu, &ntu,
                &ntu, &ktemp, &eps, &searchtype, &treetype, &radius,
                &nn_index[0], &distances[0]);
                
                ngbrs = arma::imat(&nn_index[0],ktemp,X.n_cols,true,false).t();
                ngbrs = ngbrs.tail_cols(ktemp - 1); //exclude self
  }
  
  // now select the k nearest neighbors among the 2k calculated ones
  // after recalculation of the 2k true correlations on FULL dataset
  
  
  
  robcors = arma::zeros<arma::mat>(ngbrs.n_rows, ngbrs.n_cols);
  
  for (unsigned int i = 0; i < X.n_cols; i++) {
    arma::ivec tempngbrs = ngbrs.row(i).t(); // ngbrs of variable colNumber
    arma::vec temprobcors(tempngbrs.size(), arma::fill::zeros);
    arma::vec coli = Ufull.col(i);
    arma::uvec finiteIndsi = arma::find_finite(coli);
    
    for (unsigned int j = 0; j < tempngbrs.size(); j++) {
      arma::vec colj = Ufull.col(tempngbrs(j));
      arma::uvec finiteIndsj = arma::find_finite(colj);
      arma::uvec finiteIndsij = DDC::vinter(finiteIndsi, finiteIndsj);
      temprobcors(j) = arma::as_scalar(arma::cor(coli(finiteIndsij),
                                       colj(finiteIndsij)));
    }
    arma::uvec highestCors = arma::sort_index(arma::abs(temprobcors), "descend");
  
    ngbrs.row(i) = tempngbrs(highestCors).t();
    robcors.row(i) = temprobcors(highestCors).t();
  }
  ngbrs = ngbrs.head_cols(k);
  robcors = robcors.head_cols(k);
  
  if (includeSelf == 1) {
    // If you really want to include column j among its own neighbors :
    ngbrs = join_rows(arma::regspace<arma::ivec>(0, (X.n_cols - 1), X.n_cols), ngbrs); 
    robcors = join_rows(arma::ones<arma::vec>(X.n_cols), robcors);
  }
  
  result.ngbrs = arma::conv_to<arma::umat>::from(ngbrs);
  result.robcorrs = robcors;
  
  
  return(result);
}
