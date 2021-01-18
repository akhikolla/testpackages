#include "Forest.h"
#include "globals.h"

int Forest::trainRF(std::vector<std::shared_ptr<Tree> >& trees,
                    const arma::umat& mat1Z,
                    const arma::mat& mat1f,
                    const arma::field<arma::umat>& mat2Zf,
                    const arma::umat& range0,
                    const arma::umat& ids,
                    const arma::uvec& e) {
  int n = mat1Z.n_rows;
  for(size_t i = 0; i != NUM_TREE; i++) {
    arma::field<arma::umat>  mat2Zsub(K);
    for(size_t k = 0; k != K; k++) {
      int nr = mat2Zf(k).n_rows;
      mat2Zsub(k) = mat2Zf(k).rows( ids( ids.n_rows*i + find(ids.col(i) >= n-nr) ) + nr - n );
    }
    //arma::umat mat1Zsub = mat1Z.rows( ids.col(i) );
    //arma::mat mat1fsub = mat1f.cols( ids.col(i) );
    trees.push_back( train( mat1Z.rows( ids.col(i) ),
                            mat1f.cols( ids.col(i) ),
                            mat2Zsub, range0,
                            e( ids.col(i) )) );
  }
  return 1;
}

std::shared_ptr<Tree> Forest::train(const arma::umat& mat1Z,
                                    const arma::mat& mat1f,
                                    const arma::field<arma::umat>& mat2Zf,
                                    const arma::umat& range0,
                                    const arma::uvec& e) const
{
  int n = mat1Z.n_rows;
  int P = range0.n_cols;
  arma::ucube ranges = arma::zeros<arma::ucube>(MAX_NODE, P, 2);
  arma::uvec left_childs = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec right_childs = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec split_vars = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec split_values = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec isLeaf = arma::zeros<arma::uvec>(MAX_NODE);
  ranges.row(0) = range0.t();
  arma::field<arma::uvec> nodeSampleY(MAX_NODE);
  nodeSampleY(0) = arma::regspace<arma::uvec>(0, n-1);
  arma::field<arma::uvec> nodeSample(K,MAX_NODE);
  for(size_t k = 0; k < K; k++) {
    nodeSample(k,0) = arma::regspace<arma::uvec>(0, mat2Zf(k).n_rows-1);
  }
  if(spCriterion == 1) {
    size_t ndcount = 0;
    size_t countsp = 0;
    int end = 0;
    while(end == 0) {
      end = split_DICON(mat1Z, mat1f, mat2Zf,
                        left_childs, right_childs,
                        split_vars, split_values, isLeaf,
                        ranges, nodeSampleY, nodeSample,
                        countsp, ndcount, e);
      if(ndcount >= MAX_NODE - 2) {
        isLeaf( arma::find(left_childs == 0) ).ones();
        break;
      }
    }
    arma::uvec nonEmpty = arma::regspace<arma::uvec>(0, ndcount);
    std::shared_ptr<Tree> tr(new Tree(left_childs(nonEmpty), right_childs(nonEmpty),
                                      split_vars(nonEmpty), split_values(nonEmpty), isLeaf(nonEmpty)));
    return tr;
  } else {
    arma::mat fmat = arma::zeros<arma::mat>(K, MAX_NODE);
    arma::umat Smat = arma::zeros<arma::umat>(K, MAX_NODE);
    size_t ndcount = 0;
    size_t countsp = 0;
    int end = 0;
    while(end == 0) {
      end = split_ICON(mat1Z, mat1f, mat2Zf,
                       left_childs, right_childs,
                       split_vars, split_values, isLeaf,
                       fmat, Smat, 
                       ranges, nodeSampleY, nodeSample,
                       countsp, ndcount, e);
      if(ndcount >= MAX_NODE - 2) {
        isLeaf( arma::find(left_childs == 0) ).ones();
        break;
      }
    }
    arma::uvec nonEmpty = arma::regspace<arma::uvec>(0, ndcount);
    //Rcpp::Rcout << "2";
    std::shared_ptr<Tree> tr(new Tree(left_childs(nonEmpty),
                                      right_childs(nonEmpty),
                                      split_vars(nonEmpty),
                                      split_values(nonEmpty),
                                      isLeaf(nonEmpty)));
    return tr;
  }
}


int Forest::split_DICON(const arma::umat& mat1Z,
                        const arma::mat& mat1f,
                        const arma::field<arma::umat>& mat2Zf, // dat
                        arma::uvec& left_childs,
                        arma::uvec& right_childs,
                        arma::uvec& split_vars,
                        arma::uvec& split_values,
                        arma::uvec& isLeaf, // tree
                        arma::ucube& ranges,
                        arma::field<arma::uvec>& nodeSampleY,
                        arma::field<arma::uvec>& nodeSample,
                        size_t& countsp,
                        size_t& ndcount,
                        const arma::uvec& e) const {
  int end = 0;
  int varsp = -1;
  int cutsp = 0;
  size_t nd = countsp;
  while(varsp == -1 && countsp <= ndcount) {
    nd = countsp;
    arma::ivec bestSp = find_split_DICON(nd,
                                         mat1Z, mat1f, mat2Zf,
                                         ranges, nodeSampleY, nodeSample, e);
    varsp = bestSp(1);
    cutsp = bestSp(2);
    if(varsp == -1) {
      isLeaf(nd) = 1;
      while(countsp <= ndcount) {
        countsp++;
        if(isLeaf(countsp) == 0) break;
      }
    }
  }
  int n = mat1Z.n_rows;
  if(varsp != -1) {
    split_vars(nd) = varsp;
    split_values(nd) = cutsp;
    arma::uword ndc1 = ndcount + 1;
    arma::uword ndc2 = ndcount + 2;
    left_childs(nd) = ndc1;
    right_childs(nd) = ndc2;
    arma::uvec nodeSampleYnd = std::move(nodeSampleY(nd));
    arma::uvec zvarspsub = mat1Z( varsp*n + nodeSampleYnd );
    nodeSampleY(ndc1) = nodeSampleYnd( arma::find(zvarspsub <=cutsp) );
    nodeSampleY(ndc2) = nodeSampleYnd( arma::find(zvarspsub >cutsp) );
    for(size_t k = 0; k < K; k++) {
      arma::uvec nodeSampleknd = std::move(nodeSample(k,nd));
      arma::uvec zvarspsub = mat2Zf(k)( varsp*mat2Zf(k).n_rows + nodeSampleknd );
      nodeSample(k, ndc1) = nodeSampleknd( arma::find(zvarspsub<=cutsp) );
      nodeSample(k, ndc2) = nodeSampleknd( arma::find(zvarspsub>cutsp) );
      // if(nodeSample(k, ndc1).size() < MIN_SPLIT1)
      // {isLeaf(ndc1) = 1;}
      // if(nodeSample(k, ndc2).size() < MIN_SPLIT1)
      // {isLeaf(ndc2) = 1;}
    }
    //if( arma::sum( e(nodeSampleY(ndc1))) < MIN_SPLIT2 )
    //if( nodeSampleY(ndc1).size() < MIN_SPLIT2 && isLeaf(ndc1) == 1)
    //if( arma::sum( e(nodeSampleY(ndc1))) < MIN_SPLIT2 &&  isLeaf(ndc1) )
    if(nodeSample(0, ndc1).size() < MIN_SPLIT1) {
      isLeaf(ndc1) = 1;
    } else {
      isLeaf(ndc1) = 0;
    }
    if(nodeSample(0, ndc2).size() < MIN_SPLIT1) {
      isLeaf(ndc2) = 1;
    } else {
      isLeaf(ndc2) = 0;
    }
    // if(nodeSample(0, ndc1).size() < MIN_SPLIT1)
    // {isLeaf(ndc1) = 1;}
    // if(nodeSample(0, ndc2).size() < MIN_SPLIT1)
    // {isLeaf(ndc2) = 1;}
    // if( isL1 == 1 && isLeaf(ndc1) == 1)
    // {
    //   isLeaf(ndc1) = 1;
    // }else{
    //   isLeaf(ndc1) = 0;
    // }
    //
    // if( isL2 == 1 && isLeaf(ndc2) == 1)
    // {
    //   isLeaf(ndc2) = 1;
    // }else{
    //   isLeaf(ndc2) = 0;
    // }
    //arma::sum( e(nodeSampleY(ndc1))) < MIN_SPLIT2 &&
    // if( nodeSample(0, ndc1).size() < MIN_SPLIT1 )
    // {
    //   isLeaf(ndc1) = 1;
    // }
    // // arma::sum( e(nodeSampleY(ndc2))) < MIN_SPLIT2 &&
    // if( nodeSample(0, ndc2).size() < MIN_SPLIT1)
    // {
    //   isLeaf(ndc2) = 1;
    // }
    ranges.row(ndc1) = ranges.row(nd);
    ranges.row(ndc2) = ranges.row(nd);
    ranges(ndc2,varsp,0) = cutsp+1;
    ranges(ndc1,varsp,1) = cutsp;
    ndcount += 2;
    while(countsp <= ndcount) {
      countsp++;
      if(isLeaf(countsp) == 0) break;
    }
  } else {
    end = 1;
  }
  return end;
}

arma::ivec Forest::find_split_DICON(size_t nd,
                                    const arma::umat& mat1Z,
                                    const arma::mat& mat1f,
                                    const arma::field<arma::umat>& mat2Zf, // dat
                                    const arma::ucube& ranges,
                                    const arma::field<arma::uvec>& nodeSampleY,
                                    const arma::field<arma::uvec>& nodeSample,
                                    const arma::uvec& e) const {
  int P = mat1Z.n_cols;
  int n = mat1Z.n_rows;
  int varsp = -1;
  int cutsp = 0;
  double dICONmax = 0;
  double dICONTemp = 0;
  arma::uvec spSet = arma::shuffle( arma::regspace<arma::uvec>(0,P-1) );
  for(auto p : spSet.head(mtry)) {
    arma::uvec indY = nodeSampleY(nd)( sort_index( mat1Z(p*n + nodeSampleY(nd)) ));
    arma::field<arma::uvec> indp(K);
    arma::uvec SRSum = arma::zeros<arma::uvec>(K);
    for(size_t k = 0; k < K; k++) {
      arma::uvec zpsub = mat2Zf(k)( p*mat2Zf(k).n_rows + nodeSample(k,nd) );
      indp(k) = nodeSample(k,nd)(sort_index(zpsub));
      SRSum(k) = zpsub.size();
    }
    arma::vec fLSum = arma::zeros<arma::vec>(K);
    arma::uvec SLSum = arma::zeros<arma::uvec>(K);
    arma::vec fRSum = sum(mat1f.cols(indY),1);
    int j = 0;
    arma::uvec jv = arma::zeros<arma::uvec>(K);
    int nj = indY.size();
    int nel = 0;
    // int nelr = arma::sum(e( indY ) );
    arma::uvec rangeCut = arma::regspace<arma::uvec>(ranges(nd, p, 0),
                                                     ranges(nd, p, 1));
     //arma::vec den = (fLSum + fRSum)%(SLSum + SRSum);
     //den( arma::find(den == 0) ).ones();
    for(auto cu : rangeCut) {
      while(j < nj) {
        int indYj = indY(j);
        if( mat1Z(indYj, p) == cu) {
          arma::vec df = mat1f.col( indYj );
          fLSum = fLSum + df;
          fRSum = fRSum - df;
          nel += e( indYj );
          j++;
        } else {
          break;
        }
      }
      for(size_t k = 0; k < K; k++) {
        arma::uvec indpk = indp(k);
        while(jv(k) < indpk.size()) {
          if(mat2Zf(k)( indpk( jv(k) ) , p) == cu) {
            SLSum(k)++;
            SRSum(k)--;
            jv(k)++;
          } else {
            break;
          }
        }
      }
      // if(  (SLSum(0) < MIN_NODE1 || SRSum(0) < MIN_NODE1) )
      //if( (SLSum(0) < MIN_NODE1 || SRSum(0) < MIN_NODE1) && ( nel < MIN_NODE2 || nelr-nel < MIN_NODE2) ) //|| (SLSum(0) > 9*SRSum(0) || SRSum(0) > 9*SLSum(0))
      //if( (( nel < MIN_NODE2 || nelr-nel < MIN_NODE2) && (SLSum.min() < MIN_NODE1 || SRSum.min() < MIN_NODE1) )  )
      if( (SLSum(0) < MIN_NODE1 || SRSum(0) < MIN_NODE1) ) {
        dICONTemp = 0;
      } else {
        //arma::vec hL = fLSum / SLSum;
        //arma::vec hR = fRSum / SRSum;
        //hL(arma::find(SLSum == 0)).zeros();
        //hR(arma::find(SRSum == 0)).zeros();
        //dICONTemp = abs(arma::sum( (hL - hR) ));
        dICONTemp = arma::sum( abs(fLSum%SRSum - fRSum%SLSum) );
      }
      if(dICONTemp>dICONmax) {
        dICONmax = dICONTemp;
        varsp = p;
        cutsp = cu;
      }
    }
  }
  arma::ivec vecsp(3);
  if(varsp == -1) {
    vecsp(0) = 0;
    vecsp(1) = -1;
    vecsp(2) = 0;
  } else {
    vecsp(0) = 1;
    vecsp(1) = varsp;
    vecsp(2) = cutsp;
  }
  return vecsp;
}


int Forest::split_ICON(const arma::umat& mat1Z,
                       const arma::mat& mat1f,
                       const arma::field<arma::umat>& mat2Zf, // dat
                       arma::uvec& left_childs,
                       arma::uvec& right_childs,
                       arma::uvec& split_vars,
                       arma::uvec& split_values,
                       arma::uvec& isLeaf,
                       arma::mat& fmat,
                       arma::umat& Smat,// tree
                       arma::ucube& ranges,
                       arma::field<arma::uvec>& nodeSampleY,
                       arma::field<arma::uvec>& nodeSample,
                       size_t& countsp,
                       size_t& ndcount,
                       const arma::uvec& e) const {
  int end = 0;
  int varsp = -1;
  int cutsp = 0;
  size_t nd = countsp;
  while(varsp == -1 && countsp <= ndcount) {
    nd = countsp;
    arma::ivec bestSp(3);
    bestSp = find_split_ICON(nd,
                             mat1Z, mat1f, mat2Zf, isLeaf,
                             ranges, nodeSampleY, nodeSample,
                             fmat, Smat, ndcount, e);
    varsp = bestSp(1);
    cutsp = bestSp(2);
    if(varsp == -1) {
      isLeaf(nd) = 1;
      while(countsp <= ndcount) {
        countsp++;
        if(isLeaf(countsp) == 0) break;
      }
    }
  }
  int n = mat1Z.n_rows;
  if(varsp != -1) {
    split_vars(nd) = varsp;
    split_values(nd) = cutsp;
    arma::uword ndc1 = ndcount + 1;
    arma::uword ndc2 = ndcount + 2;
    left_childs(nd) = ndc1;
    right_childs(nd) = ndc2;
    arma::uvec nodeSampleYnd = std::move(nodeSampleY(nd));
    arma::uvec zvarspsub = mat1Z( varsp*n + nodeSampleYnd );
    nodeSampleY(ndc1) = nodeSampleYnd( arma::find(zvarspsub <=cutsp) );
    nodeSampleY(ndc2) = nodeSampleYnd( arma::find(zvarspsub >cutsp) );
    for(size_t k = 0; k < K; k++) {
      arma::uvec nodeSampleknd = std::move(nodeSample(k,nd));
      arma::uvec zvarspsub = mat2Zf(k)( varsp*mat2Zf(k).n_rows + nodeSampleknd );
      nodeSample(k, ndc1) = nodeSampleknd( arma::find(zvarspsub<=cutsp) );
      nodeSample(k, ndc2) = nodeSampleknd( arma::find(zvarspsub>cutsp) );
    }
    if(nodeSample(0, ndc1).size() < MIN_SPLIT1) isLeaf(ndc1) = 1;
    if(nodeSample(0, ndc2).size() < MIN_SPLIT1) isLeaf(ndc2) = 1;
    ranges.row(ndc1) = ranges.row(nd);
    ranges.row(ndc2) = ranges.row(nd);
    ranges(ndc2,varsp,0) = cutsp+1;
    ranges(ndc1,varsp,1) = cutsp;
    ndcount = ndcount+2;
    while(countsp <= ndcount) {
      countsp++;
      if(isLeaf(countsp) == 0) break;
    }
  } else {
    end = 1;
  }
  return end;
}

arma::ivec Forest::find_split_ICON(size_t nd,
                                   const arma::umat& mat1Z,
                                   const arma::mat& mat1f,
                                   const arma::field<arma::umat>& mat2Zf, // dat
                                   const arma::uvec& isLeaf,
                                   const  arma::ucube& ranges,
                                   const arma::field<arma::uvec>& nodeSampleY,
                                   const arma::field<arma::uvec>& nodeSample,
                                   arma::mat& fmat,
                                   arma::umat& Smat,
                                   int ndcount,
                                   const arma::uvec& e) const {
  int P = mat1Z.n_cols;
  int n = mat1Z.n_rows;
  int varsp = -1;
  int cutsp = 0;
  double dICONmax = 0;
  double dICONTemp = 0;
  arma::mat fmatTerm = fmat.cols(arma::find(isLeaf == 1));
  arma::umat SmatTerm = Smat.cols(arma::find(isLeaf == 1));
  arma::uvec spSet = arma::shuffle( arma::regspace<arma::uvec>(0,P-1) );
  for(auto p : spSet.head(mtry)) {
    arma::uvec indY = nodeSampleY(nd)( sort_index( mat1Z(p*n + nodeSampleY(nd)) ));
    arma::field<arma::uvec> indp(K);
    arma::uvec SRSum = arma::zeros<arma::uvec>(K);
    for(size_t k = 0; k < K; k++) {
      arma::uvec zpsub = mat2Zf(k)( p*mat2Zf(k).n_rows + nodeSample(k,nd) );
      indp(k) = nodeSample(k,nd)(sort_index(zpsub));
      SRSum(k) = zpsub.size();
    }
    arma::vec fLSum = arma::zeros<arma::vec>(K);
    arma::uvec SLSum = arma::zeros<arma::uvec>(K);
    arma::vec fRSum = sum(mat1f.cols(indY),1);
    int j = 0;
    arma::uvec jv = arma::zeros<arma::uvec>(K);
    int nj = indY.size();
    arma::uvec rangeCut = arma::regspace<arma::uvec>(ranges(nd, p, 0),
                                                     ranges(nd, p, 1));
    int nel = 0;
    // int nelr = arma::sum(e( indY ) );
    for(auto cu : rangeCut) {
      while(j < nj) {
        int indYj = indY(j);
        size_t z = mat1Z(indYj, p);
        if(z == cu) {
          arma::vec df = mat1f.col( indYj );
          fLSum = fLSum + df;
          fRSum = fRSum - df;
          nel += e( indYj );
          j++;
        } else {
          break;
        }
      }
      for(size_t k = 0; k < K; k++) {
        arma::uvec indpk = indp(k);
        while(jv(k) < indpk.size()) {
          if(mat2Zf(k)( indpk( jv(k) ) , p) == cu) {
            SLSum(k)++;
            SRSum(k)--;
            jv(k)++;
          } else {
            break;
          }
        }
      }
      if( (SLSum(0) < MIN_NODE1 || SRSum(0) < MIN_NODE1) ) {
        dICONTemp = 0;
      } else {
        dICONTemp = arma::sum(abs(fLSum%SRSum - fRSum%SLSum))/2;
        for(size_t i = 0; i < fmatTerm.n_cols; i++) {
          if(i != nd ) {
            for(size_t k = 0; k < K; k++) {
              if( (fLSum(k)+fRSum(k))*SmatTerm(k,i) <=  fmatTerm(k,i)*(SLSum(k)+SRSum(k)) ) {
                if( fmatTerm(k,i)*SLSum(k) <  fLSum(k)*SmatTerm(k,i) ) 
                {// lambdaR < lambdand <= lambdai < lambdaL
                  dICONTemp = dICONTemp +  fLSum(k)*SmatTerm(k,i) - fmatTerm(k,i)*SLSum(k) ;
                }
                if( fmatTerm(k,i)*SRSum(k) <  fRSum(k)*SmatTerm(k,i))
                {// lambdaL < lambdand <= lambdai < lambdaR
                  dICONTemp = dICONTemp + fRSum(k)*SmatTerm(k,i) - fmatTerm(k,i)*SRSum(k) ;
                }
              }else{
                if( fmatTerm(k,i)*SLSum(k) >  fLSum(k)*SmatTerm(k,i) )
                {// lambdaL < lambdai < lambdand < lambdaR
                  dICONTemp = dICONTemp + fmatTerm(k,i)*SLSum(k) - fLSum(k)*SmatTerm(k,i) ;
                }
                if( fmatTerm(k,i)*SRSum(k) >  fRSum(k)*SmatTerm(k,i) )
                {// lambdaR < lambdai < lambdand < lambdaL
                  dICONTemp = dICONTemp + fmatTerm(k,i)*SRSum(k) - fRSum(k)*SmatTerm(k,i) ;
                }
              }
            }
          }
        }
      }
      if(dICONTemp>dICONmax) {
        dICONmax = dICONTemp;
        varsp = p;
        cutsp = cu;
        fmat.col(ndcount + 1) = fLSum;
        fmat.col(ndcount + 2) = fRSum;
        Smat.col(ndcount + 1) = SLSum;
        Smat.col(ndcount + 2) = SRSum;
      }
    }
  }
  arma::ivec vecsp(3);
  if(varsp == -1) {
    vecsp(0) = 0;
    vecsp(1) = -1;
    vecsp(2) = 0;
  } else {
    vecsp(0) = 1;
    vecsp(1) = varsp;
    vecsp(2) = cutsp;
  }
  return vecsp;
}
