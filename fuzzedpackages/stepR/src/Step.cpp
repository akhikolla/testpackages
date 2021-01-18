#include "Step.h"
#include <cstdio>

/***************
* class Step
* virtual class allowing to fit step functions to data
* Thomas Hotz, 2007-2011
***************/

/*************
* constructor for n data points
****************/
Step::Step(unsigned int n) : N(n) {}
Step::Step(unsigned int n, double* xlb, double* xub) : N(n), lb(xlb), ub(xub) {}

/*************
* cost
* calculate cost of a block without considering boundary conditions, simply by setting the boundaries to -Inf and Inf, resp.
****************/
// double Step::cost(unsigned int startIndex, unsigned int endIndex) const {
//   return costBound(startIndex, endIndex, R_NegInf, R_PosInf);
// }

/*************
* costBound
* calculate cost of a block with boundary conditions, needs to be implemented for each derived class
****************/
double Step::costBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const {
  error("Step::costBound has to be overwritten!");
  return R_NaReal;
}

/*************
* estBound
* calculate eatimate of a block with boundary conditions, needs to be implemented for each derived class
****************/
double Step::estBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const {
  error("Step::estBound has to be overwritten!");
  return R_NaReal;
}

/*************
* findCandidate
* finds a candidate jump between two existing jumps
*
* in:
* prev : the previous jump
* next : the next jump
* 
* out:
* the candidate jump
****************/
Jump Step::findCandidate(const Jump& prev, const Jump& next) const {
  double pim; // the potential improvement
  int cind = -1; // the index of the currently best improvement
  double cim = -1; // the currently best improvement
  double fullcost = cost(prev.rightIndex + 1, next.rightIndex); // the cost of the entire block
  
  if(next.rightIndex - prev.rightIndex < 2) {
    error("No room left for candidate!");
  } else {
    for(int i = prev.rightIndex + 1; i < next.rightIndex; i++) {
      pim = fullcost - ( cost(prev.rightIndex + 1, i) + cost(i + 1, next.rightIndex) );
      if(pim > cim) {
        // improvement
        cim = pim;
        cind = i;
      }
    }
/*    if(cim <= -1) {
//       Rprintf(" fullcost = %f, cost(%d:%d) = %f, cost(%d:%d) = %f, pim = %f\n", fullcost, prev.rightIndex + 2, next.rightIndex, cost(prev.rightIndex + 1, next.rightIndex - 1), next.rightIndex + 1, next.rightIndex + 1, cost(next.rightIndex, next.rightIndex), pim);
      error("No potential candidate found between %d and %d!", prev.rightIndex + 2, next.rightIndex + 1); // R style indices
    }*/
  }
/*  Rprintf("   Found a candidate\n");
  Rprintf("   rightIndex = %d, improve = %f\n", cind, cim);*/
  return Jump(NA_INTEGER, cind, cim);
}

/*************
* forward
* compute forward selection for block signals
*
* in:
* maxBlocks : an integer giving the maximal number of blocks
*
* out:
* SEXP : data.frame comprising rightIndex, number, depth and improve of the candidates
****************/
SEXP Step::forward(unsigned int maxBlocks) const {
  // check maxNum
  if(maxBlocks < 1) error("there must be at least one block allowed");
  if(maxBlocks > N) error("there may not be more than N blocks");
  
  // initialize computation
  Jump before = Jump(); // "before" the observed data
  Jump prev = before; // the previous jump, initially "before" the observed data
  Jump next = Jump(0, N - 1, 0); // the following jump, initially the "end" of the observed data
  BinTree<Jump> bt(next); // start with the "end", i.e. with the block comprising all data
  unsigned int num = 1; // number of jumps selected
  
  // allocate storage, freed by R
  int* number = (int*) R_alloc(maxBlocks, sizeof(int));
  int* depth = (int*) R_alloc(maxBlocks, sizeof(int));
  int* rightIndex = (int*) R_alloc(maxBlocks, sizeof(int));
  double* improve = (double*) R_alloc(maxBlocks, sizeof(double));
  
  if(N > 1) {
    bt.addLeft(findCandidate(prev, next)); // find potential candidate
    
    // compute tree
    Jump cand; // the currently best candidate
    
    while(num < maxBlocks) {
      // find next jump among potential candidates
      bt.first();
      cand = before; // "before" has negative improvement
      do {
        if(bt.getValue().number == NA_INTEGER && bt.getValue().improve >= cand.improve) {
          // jump has not yet been selected but has larger improvement: new candidate
          cand = bt.getValue();
        }
      } while(bt.next());
      
      if(cand.rightIndex < 0) break; //error("No candidate found!");
      
      // find best candidate again
      bt.first();
      while(bt.getValue().rightIndex != cand.rightIndex) {
        if(!bt.next()) error("Could not find candidate %d again!", cand.rightIndex);
      }
      
      // select best candidate
      cand.number = num;
      bt.setValue(cand);
//       Rprintf(" num = %d, rightIndex = %d\n", num, cand.rightIndex);
      num++;
      
      // add two more potential candidates (note that there will never be two potential candidates next to each other)
      if(bt.previous()) {
        prev = bt.getValue();
        bt.next(); // move back to selected
      } else {
        prev = before; // we were at the beginning
      }
      if(cand.rightIndex - prev.rightIndex > 1) {
        bt.addLeft(findCandidate(prev, cand)); // add potential candidate to the left
        bt.right(); // move back to selected
      }
      bt.next();
      next = bt.getValue(); // add other candidate left of next jump to balance tree
      if(next.rightIndex - cand.rightIndex > 1) {
        bt.addLeft(findCandidate(cand, next)); // add potential candidate to the right
      }
    }
  }
  
  // extract results
  double totalcost;
  flattenTree(&bt, number, depth, rightIndex, improve, totalcost);
  
  // return data.frame
  SEXP ret = allocVector(VECSXP, 4);
  PROTECT(ret);
  
  SEXP names = allocVector(STRSXP, 4);
  PROTECT(names);
  SET_STRING_ELT(names, 0, mkChar("rightIndex"));
  SET_STRING_ELT(names, 1, mkChar("number"));
  SET_STRING_ELT(names, 2, mkChar("depth"));
  SET_STRING_ELT(names, 3, mkChar("improve"));
  ret = namesgets(ret, names);
  
  SEXP rownames = allocVector(STRSXP, num);
  PROTECT(rownames);
  char buffer [8];
  for(unsigned int i = 0; i < num; i++) {
    std::sprintf(buffer, "%d", i + 1);
    SET_STRING_ELT(rownames, i, mkChar(buffer));
  }
  setAttrib(ret, R_RowNamesSymbol, rownames);
    
  SEXP sclass = allocVector(STRSXP, 1);
  PROTECT(sclass);
  SET_STRING_ELT(sclass, 0, mkChar("data.frame"));
  classgets(ret, sclass);
  
  SEXP ri = allocVector(INTSXP, num);
  SET_VECTOR_ELT(ret, 0, ri);
  int *xri = INTEGER(ri);
  
  SEXP n = allocVector(INTSXP, num);
  SET_VECTOR_ELT(ret, 1, n);
  int *xn = INTEGER(n);
  
  SEXP d = allocVector(INTSXP, num);
  SET_VECTOR_ELT(ret, 2, d);
  int *xd = INTEGER(d);
  
  SEXP im = allocVector(REALSXP, num);
  SET_VECTOR_ELT(ret, 3, im);
  double *xim = REAL(im);
  
  for(unsigned int i = 0; i < num; i++) {
    xn[i] = number[i];
    xd[i] = depth[i];
    xri[i] = rightIndex[i] + 1;  // turn from C-style index into R-style index
    xim[i] = improve[i];
//     Rprintf(" number = %d, depth = %d, rightIndex = %d\n", number[i], depth[i], rightIndex[i]);
  }
  
  SEXP scost = allocVector(REALSXP, 1);
  PROTECT(scost);
  REAL(scost)[0] = totalcost;
  setAttrib(ret, install("cost"), scost);
  
  UNPROTECT(5);
  
  return ret;

}

/*************
* flattenTree
* flattens SBinTree and stores values in vectors
*
* in:
* bt : the result tree
* number, depth, rightIndex, improve : arrays of length maxNum
*
* out:
* number : order in which jumps have been selected
* depth : depth of the node
* rightIndex : rightIndex of the block
* improve : the improvement this jump brought about
* totalcost : the total cost if all candidates were used
****************/
void Step::flattenTree(BinTree<Jump>* bt, int* number, int* depth, int* rightIndex, double* improve, double &totalcost) const {
//   Rprintf("Flattening Tree\n");
  // go to first leaf
  bt->first();
  Jump value = bt->getValue();
//   Rprintf(" depth = %d, rightIndex = %d\n", bt->depth(), value.rightIndex);
  int i = 0;
  int leftIndex = 0;
  totalcost = 0;
  if(value.number != NA_INTEGER) { // if indeed selected
    number[i] = value.number;
    depth[i] = bt->depth();
    rightIndex[i] = value.rightIndex;
    improve[i] = value.improve;
    totalcost += cost(leftIndex, rightIndex[i]);
    leftIndex = rightIndex[i] + 1;
    i++;
  }
  
  while(bt->next()) {
    value = bt->getValue();
//     Rprintf(" number = %d, depth = %d, rightIndex = %d\n", value.number, bt->depth(), value.rightIndex);
    if(value.number != NA_INTEGER) { // if indeed selected
      number[i] = value.number;
      depth[i] = bt->depth();
      rightIndex[i] = value.rightIndex;
      improve[i] = value.improve;
      totalcost += cost(leftIndex, rightIndex[i]);
      leftIndex = rightIndex[i] + 1;
      i++;
    }
  }
}

/*************
* path
* compute solution path of Potts minimisers
*
* in:
* maxBlocks : an integer giving the maximal number of blocks
*
* out:
* SEXP : a list conprimising the right indices of the solution path, and the associated costs
****************/
SEXP Step::path(unsigned int maxBlocks) const {
  // allocate storage
  TriArrayFF<double> D(N); // the cost(m,n) for m <= s <= n (has first index running fastest, see below)
  TriArray<double> B(N - 1); // the minimal values: B[k, n] is the cost of the optimal fit up to time n+1 using k+1 jumps
  TriArray<int> r(N - 1); // the index of the optimal (last) jump corresponding to B[k, n], i.e.
                           // B[k, n] = B[k - 1, r[k, n] - 1] + D[r[k, n] + 1, n + 1]
  TriArrayFF<int> p(N - 1); // the solution path, i.e. p[i, k] is the (i+1)th jump in the solution having k+1 jumps
  
  // check maxBlocks
  if(maxBlocks < 1) error("there must be at least one block allowed");
  if(maxBlocks > N) error("there may not be more than N blocks");
  
  // compute D
  for(unsigned int n = 0; n < N; n++) {
    D(0, n) = cost(0, n); // no jump before 0
    for(unsigned int m = 1; m <= n; m++) {
      D(m, n) = cost(m, n);
    }
  }
  
  // compute B, r
  unsigned int candR; // candidate r
  double candB; // candidate B
  for(unsigned int n = 0; n < N - 1; n++) {
    // initialise if one jump up to n (k = 0)
    r(0, n) = 0;
    B(0, n) = D(0, 0) + D(1, n + 1);
    for(candR = 1; candR <= n; candR++) {
      candB = D(0, candR) + D(candR + 1, n + 1); // here we want first index of D running fastest!
      if(candB < B(0, n)) {
        r(0, n) = candR;
        B(0, n) = candB;
      }
    }
  }
  for(unsigned int k = 1; k < maxBlocks - 1; k++) {
    // two and more jumps
    for(unsigned int n = k; n < N - 1; n++) {
      // search for candidate r
      r(k, n) = k; // initialise
      B(k, n) = B(k - 1, k - 1) + D(k + 1, n + 1); // candR is a candidate for the kth jump and therefore >= k
      for(candR = k + 1; candR <= n; candR++) {
        candB = B(k - 1, candR - 1) + D(candR + 1, n + 1); // here we want first index of D running fastest!
        if(candB < B(k, n)) {
          r(k, n) = candR;
          B(k, n) = candB;
        }
      }
    }
  }
  
  // compute solution
  unsigned int num = maxBlocks; // the number of blocks that have been found
  for(unsigned int k = 0; k < maxBlocks - 1; k++){
    if(B(k, N - 2) == R_PosInf) { // k jumps are impossible
      num = k + 1;
      break;
    }
    p(k, k) = r(k, N - 2); // the last of the k jumps (up to N - 1) is given by r[k, N - 2]
    for(unsigned int i = 1; i <= k; i++) { // we update the (k - i + 1)th jump out of k + 1 (going backwards)
      p(k - i, k) = r(k - i, p(k - i + 1, k) - 1); // get the optimal (last) jump up to the latest jump found
    }
  }
  
  // return result
  SEXP ret = allocVector(VECSXP, 2); // return list
  PROTECT(ret);
  
  SEXP names = allocVector(STRSXP, 2);
  PROTECT(names);
  SET_STRING_ELT(names, 0, mkChar("path")); // contains list of vectors comprising right Indices
  SET_STRING_ELT(names, 1, mkChar("cost"));
  ret = namesgets(ret, names);
  
  SEXP path = allocVector(VECSXP, num);
  PROTECT(path);
  SET_VECTOR_ELT(ret, 0, path);
  
  SEXP retCost = allocVector(REALSXP, num);
  PROTECT(retCost);
  SET_VECTOR_ELT(ret, 1, retCost);
  double *xcost = REAL(retCost);
  
  SEXP solutionk = allocVector(INTSXP, 1); // vector of this solution's right indices 
  PROTECT(solutionk);
  SET_VECTOR_ELT(path, 0, solutionk);
  int *xsolutionk = INTEGER(solutionk);
  
  xcost[0] = D(0, N - 1); // cost without a jump, i.e. one block
  xsolutionk[0] = N; // only jump in this case is at the end
  
  for(unsigned int k = 1; k < num; k++){
    solutionk = allocVector(INTSXP, k + 1); // vector of this solution's right indices 
    SET_VECTOR_ELT(path, k, solutionk);
    xsolutionk = INTEGER(solutionk);
    
    xcost[k] = B(k - 1, N - 2);
    
    for(unsigned int i = 0; i < k; i++) {
      xsolutionk[i] = p(i, k - 1) + 1;  // turn from C-style index into R-style index
    }
    xsolutionk[k] = N; // add right end of last block, turn from C-style index into R-style index
  }
  
  UNPROTECT(5);
  
  return(ret);
}

/*************
* bounded
* compute optimal solution with minimal jumps fulfilling bounds
*
* out:
* SEXP : a list giving the solution
****************/
SEXP Step::bounded(Bounds& B) const {
  // allocate storage
  unsigned int s = 0; // number of jumps
  double* const J = (double*) R_alloc(N, sizeof(double)); // cost for optimal solution with s jumps over [0, ..., k], indexed by k
  int* const L = (int*) R_alloc(N, sizeof(int)); // last jump for optimal solution with s jumps over [0, ..., k], indexed by k
  double* const V = (double*) R_alloc(N, sizeof(double)); // estimate on last constant interval [l+1, ..., k], indexed by k
  double curJ = R_PosInf; // the cost for current solution with last jump at l
  double curD = R_PosInf; // cost for constant estimate on interval [l+1, ..., k]
  unsigned int* const K = (unsigned int*) R_alloc(N + 2, sizeof(unsigned int)); // maximal k s.t. s jumps are sufficient for feasible solution over [0, ... ,k], index s runs from -2 to n-1, requires offset of 2:
  unsigned int const Koffset = 2;
  int* const KL = (int*) R_alloc(N - 1, sizeof(int)); // minimal k s.t. s jumps are unnecessary for feasible solution over [0, ... ,k-1], index s runs from 1 to n-1, requires offset of -1:
  int const KLoffset = -1;
  unsigned int k; // index for interval [0, ..., k]
  int l; // index of last jump, i.e. right index of second but last block
  
  // inititalize
  K[-2 + Koffset] = 0; K[-1 + Koffset] = 0; // for "negative number of jumps"
  for(k = 0; k < N; k++) { // find constant solution over [0, ..., k]
    J[k] = R_PosInf;
    #ifdef DEBUGbounded
    Rprintf("s = %d, k = %d\n", 0, k);
    #endif
    for(l = k - 1; l > 0; l--) { // precompute bounds on [l, k] for l > 0
      B.current(l, k);
    }
    curD = costBound(0, k, B.current(0, k));
    #ifdef DEBUGbounded
    Rprintf("  l = %d, curD = %4.3e, J[l] = %4.3e\n", -1, curD, 0);
    #endif
    if(curD == R_PosInf) break; // no constant solution on [0, ..., k] possible
    curJ = curD;
    #ifdef DEBUGbounded
    Rprintf("  curJ = %4.3e, curJ < J[k] = %d\n", curJ, curJ < J[k]);
    #endif
    if(curJ < J[k]) { // improvement
      J[k] = curJ;
      L[k] = -1;
      V[k] = estBound(0, k, B.current(0, k));
    }
    K[0 + Koffset] = k;
    #ifdef DEBUGbounded
    Rprintf("s = %d, k = %d feasible\n", 0, k);
    #endif
  }
  if(K[s + Koffset] != N - 1) { // found no feasible solution on [0, ..., N-1]
    // calculations with at least one jump
    for(s = 1; s < N; s++) { // try to find solution with s jumps
      for(k = K[s - 1 + Koffset] + 1; k < N; k++) { // find solution over [0, ..., k]
        J[k] = R_PosInf;
        #ifdef DEBUGbounded
        Rprintf("s = %d, k = %d\n", s, k);
        #endif
        for(l = k - 1; l >= (int) K[s - 1 + Koffset] + 1; l--) { // precompute bounds on [l, k] for l > K[s - 1 + Koffset]
          B.current(l, k);
        }
        #ifdef DEBUGbounded
        Rprintf("K[s - 1 + Koffset] = %d, K[s - 2 + Koffset] = %d\n", K[s - 1 + Koffset], K[s - 2 + Koffset]);
        #endif
        for(l = (int) K[s - 1 + Koffset]; l >= (int) K[s - 2 + Koffset]; l--) { // try for last jump at l, i.e. right index of second but last block
          curD = costBound(l + 1, k, B.current(l + 1, k));
          #ifdef DEBUGbounded
          Rprintf("  l = %d, curD = %4.3e, J[l] = %4.3e\n", l, curD, J[l]);
          #endif
          if(curD == R_PosInf) break; // no constant solution on [l+1, ..., k] possible
          curJ = J[l] + curD;
          if(curJ < J[k]) { // improvement
            J[k] = curJ;
            L[k] = l;
            V[k] = estBound(l + 1, k, B.current(l + 1, k));
          }
        }
        if(J[k] == R_PosInf) break; // no feasible solution with s jumps on [0, ..., k]
        if(k == K[s - 1 + Koffset] + 1) KL[s + KLoffset] = l + 1;
        K[s + Koffset] = k;
        #ifdef DEBUGbounded
        Rprintf("s = %d, k = %d feasible\n", s, k);
        #endif
      }
      if(K[s + Koffset] == N - 1) break; // found a feasible solution on [0, ..., N-1]
    }
  }
  
  // return result
  SEXP ret = allocVector(VECSXP, 4); // return list
  PROTECT(ret);
  
  // create list which can be turned into a data.frame
  SEXP names = allocVector(STRSXP, 4);
  PROTECT(names);
  SET_STRING_ELT(names, 0, mkChar("rightEnd")); // contains list of vectors comprising right Indices
  SET_STRING_ELT(names, 1, mkChar("value"));
  SET_STRING_ELT(names, 2, mkChar("endLeftBound"));
  SET_STRING_ELT(names, 3, mkChar("endRightBound"));
  namesgets(ret, names);
/*  SEXP clas = allocVector(STRSXP, 1);
  SET_STRING_ELT(clas, 0, mkChar("data.frame"));
  classgets(ret, clas);*/
  
  SEXP rightEnd = allocVector(INTSXP, s + 1);
  PROTECT(rightEnd);
  SET_VECTOR_ELT(ret, 0, rightEnd);
  int *xrightEnd = INTEGER(rightEnd);
  
  // return cost as attribute
  SEXP xJ = ScalarReal(J[N-1]);
  PROTECT(xJ);
  SEXP costString = install("cost");
  PROTECT(costString);
  setAttrib(ret, costString, xJ);
  
  SEXP value = allocVector(REALSXP, s + 1);
  PROTECT(value);
  SET_VECTOR_ELT(ret, 1, value);
  double *xvalue = REAL(value);
  
  SEXP jumpLeftBound = allocVector(INTSXP, s + 1);
  PROTECT(jumpLeftBound);
  SET_VECTOR_ELT(ret, 2, jumpLeftBound);
  int *xjumpLeftBound = INTEGER(jumpLeftBound);
  
  SEXP jumpRightBound = allocVector(INTSXP, s + 1);
  PROTECT(jumpRightBound);
  SET_VECTOR_ELT(ret, 3, jumpRightBound);
  int *xjumpRightBound = INTEGER(jumpRightBound);
  
  k = N - 1;
  for(int i = s; i >= 0; i--) {
    #ifdef DEBUGbounded
    Rprintf("i = %d, k = %d\n", i, k);
    #endif
    xrightEnd[i] = k + 1;  // turn from C-style index into R-style index
    xvalue[i] = V[k];
    if(i == (int) s) {
      xjumpLeftBound[i] = N;
      xjumpRightBound[i] = N;
    } else {
      xjumpLeftBound[i] = KL[i] + 1;  // turn from C-style index into R-style index
      xjumpRightBound[i] = K[i + Koffset] + 1;  // turn from C-style index into R-style index
    }
    k = L[k];
  }
  
  UNPROTECT(8);
  
  return(ret);
}


// C wrapper
extern "C" {

/*************
* confBand
* function to be called from R
* computes pointwise confidence band
*
* in:
* confLeft : an integer vector, the (R-style) index of the left end of each jump's confidence interval, ending with n, i.e. one more than number of jumps
* confRight : an integer vector, the (R-style) index of the left end of each jump's confidence interval, starting with 0, i.e. one more than number of jumps
* start : for every possible left index where intervals with this left index start in the list of intervals (increasing in left indices), NA if none
* rightIndex : right indices of the intervals in the list, increasing for each left index
* lower : the lower bounds for the estimator at the respective interval
* upper : the upper bounds for the estimator at the respective interval
*
* note:
* confidence intervals for jumps may not overlap
*
* out:
* a list conprimising the lower and upper confidence bands
****************/
SEXP confBand(SEXP confLeft, SEXP confRight, SEXP start, SEXP rightIndex, SEXP lower, SEXP upper) {
  int* cl = INTEGER(confLeft); // R-style indices!
  int* cr = INTEGER(confRight); // R-style indices!
  LUBound cb; // current bound
  int k; // right index for interval
  int l; // left index for interval
  int b; // index for block
  int minl; // minimal left index for this block
  
  // check lengths
  if(length(confLeft) < 1) error("there must be at least one block");
  if(length(confLeft) != length(confRight)) error("confLeft must have same length as confRight (number of blocks)");
  if(cl[length(confLeft) - 1] != length(start)) error("confLeft must end with n, i.e. length of start");
  if(cr[0] != 0) error("confRight must start with 0");
  if(length(lower) != length(upper)) error("lower must have same length as upper");
  if(length(upper) != length(rightIndex)) error("upper must have same length as rightIndex");
  
  Bounds B = Bounds(length(start), INTEGER(start), length(lower), INTEGER(rightIndex), REAL(lower), REAL(upper));

  // allocate result
  SEXP ret = allocVector(VECSXP, 2); // return list
  PROTECT(ret);
  
  // create list which can be turned into a data.frame
  SEXP names = allocVector(STRSXP, 2);
  PROTECT(names);
  SET_STRING_ELT(names, 0, mkChar("lower")); // contains list of vectors comprising right Indices
  SET_STRING_ELT(names, 1, mkChar("upper"));
  namesgets(ret, names);

  SEXP low = allocVector(REALSXP, length(start));
  PROTECT(low);
  SET_VECTOR_ELT(ret, 0, low);
  double *xlow = REAL(low);
  
  SEXP up = allocVector(REALSXP, length(start));
  PROTECT(up);
  SET_VECTOR_ELT(ret, 1, up);
  double *xup = REAL(up);
  
  // loop over blocks, note that cl, cr are R-style indices!
  for(b = 0; b < length(confLeft); b++) {
    
    // precompute bounds on "interior" of block, i.e. between right end of left jump's confidence interval and left end of right jumps's confidence interval
    // as well as on left extension if any
    for(k = cr[b]; k < cl[b]; k++) {
      if(b > 0) {
	minl = cl[b-1];
      } else {
	minl = cr[b];
      }
      for(l = k - 1; l >= minl; l--) { // precompute bounds on [l, k] for l >= cr[b], note that [k,k] had been computed at inititalization
        #ifdef DEBUGbounded
        Rprintf("interior: b = %d, l = %d, k = %d\n", b, l, k);
        #endif
	B.current(l, k);
      }
    }
    // save result on interior
    cb = B.current(cr[b], cl[b] - 1);
    #ifdef DEBUGbounded
    Rprintf("saving: b = %d, l = %d, k = %d\n", b, cr[b], cl[b] - 1);
    #endif
    for(k = cr[b]; k < cl[b]; k++) {
      xlow[k] = cb.lower;
      xup[k] = cb.upper;
    }
    // extend boundaries from this block into left confidence interval (if any)
    if(b > 0) {
      for(l = cr[b] - 1; l >= cl[b-1]; l--) {
	cb = B.current(l, cl[b] - 1);
        #ifdef DEBUGbounded
	Rprintf("saving: b = %d, l = %d, k = %d\n", b, l, cl[b] - 1);
        #endif
	xlow[l] = fmin2(xlow[l], cb.lower); // jump either befor or after
	xup[l] = fmax2(xup[l], cb.upper); // jump either befor or after
      }
    }
    
    // extend boundaries from this block into right confidence interval (if any)
    if(b + 1 < length(confLeft)) {
      for(k = cl[b]; k < cr[b+1]; k++) {
	for(l = k - 1; l >= cr[b]; l--) { // precompute bounds on [l, k] for l >= cr[b], note that [k,k] had been computed at inititalization
	  #ifdef DEBUGbounded
	  Rprintf("extend right: b = %d, l = %d, k = %d\n", b, l, k);
	  #endif
	  B.current(l, k);
	}
	cb = B.current(cr[b], k);
        #ifdef DEBUGbounded
	Rprintf("saving: b = %d, l = %d, k = %d\n", b, cr[b], k);
        #endif
	xlow[k] = cb.lower; // fmin2(xlow[k], cb.lower);
	xup[k] = cb.upper; // fmax2(xup[k], cb.upper);
      }
    }
  }
  
  // return result
  UNPROTECT(4);
  return(ret);
}

} // end C wrapper

