#include <memory>
#include <RcppArmadillo.h>

#include "TreeGrow.h"
#include "globals.h"


double TreeGrow::get_ICONTest(const arma::uvec& isLeafTemp,
                              const arma::mat& fmat, const arma::umat& Smat,
                              const arma::mat& fmat2, const arma::umat& Smat2) const
{
  arma::uvec leafTemp = arma::find(isLeafTemp==1);
  int numLeafTemp = leafTemp.n_elem;
  arma::vec icon = arma::zeros<arma::vec>(K);
  arma::vec fSum = arma::zeros<arma::vec>(K);
  arma::vec SSum = arma::zeros<arma::vec>(K);
  for(size_t k = 0; k != K; k++) {
    for(int i = 0; i < numLeafTemp; i++) {
      int li = leafTemp(i);
      double fi = fmat(k,li);
      double Si = Smat(k,li);
      double fi2 = fmat2(k,li);
      double Si2 = Smat2(k,li);
      for(int j = 0; j < i; j++) {
        int lj = leafTemp(j);
        if(fi * Smat(k,lj) > fmat(k,lj) * Si) {
          icon(k) += fi2 * Smat2(k,lj);
        } else {
          icon(k) += fmat2(k,lj) * Si2;
        }
      }
      fSum(k) += fi2;
      SSum(k) += Si2;
      icon(k) += fi2 * Si2 / 2;
    }
  }
  icon = icon / (fSum % SSum);
  arma::uvec zeroDen = arma::find(fSum == 0 || SSum == 0);
  icon ( zeroDen ).ones();
  icon ( zeroDen ) = icon ( zeroDen )/2; // set icon(t_k) to 0.5 if denominator is 0
  return arma::sum(icon) / K;
}

std::shared_ptr<Tree> TreeGrow::trainCV(const arma::umat& mat1Z,
                                        const arma::mat& mat1f,
                                        const arma::field<arma::umat>& mat2Zf,
                                        const arma::umat& range0,
                                        const arma::uvec& e) const
{
  arma::mat fmat = arma::zeros<arma::mat>(K, MAX_NODE);
  arma::umat Smat = arma::zeros<arma::umat>(K, MAX_NODE);
  std::shared_ptr<Tree> tr = grow(mat1Z, mat1f, mat2Zf, range0, fmat, Smat, e);
  if (NUM_FOLD > 1) {
    const arma::uvec& isl = tr->get_isLeaf();
    uint numLeaf = arma::sum(isl);
    arma::field<arma::uvec> nodeSetList(numLeaf);
    arma::vec iconAll(numLeaf);
    tr->findOptimalSizekSubtree(fmat, Smat, iconAll, nodeSetList, numLeaf);
    // Initialization
    arma::uvec sizeTree = arma::regspace<arma::uvec>(1,iconAll.n_elem);
    arma::vec beta(iconAll.n_elem);
    Tree::findBeta(iconAll, beta, sizeTree); 
    arma::vec iconBeta = prune(beta, mat1Z, mat1f, mat2Zf, range0, e);
    uint qo = iconBeta.index_max();
    // Rcpp::Rcout << "iconBeta" << iconBeta << std::endl;
    // Rcpp::Rcout << "qo" << qo << std::endl;
    arma::uvec nodeSetFinal = nodeSetList(sizeTree(qo)-1);
    tr->cut(nodeSetFinal);
  }
  return tr;
}

std::shared_ptr<Tree> TreeGrow::grow(const arma::umat& mat1Z,
                                     const arma::mat& mat1f,
                                     const arma::field<arma::umat>& mat2Zf,
                                     const arma::umat& range0,
                                     arma::mat& fmat,
                                     arma::umat& Smat,
                                     const arma::uvec& e) const
{
  int n = mat1Z.n_rows;
  int P = mat1Z.n_cols;
  arma::ucube ranges = arma::zeros<arma::ucube>(MAX_NODE, P, 2);
  arma::uvec left_childs = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec right_childs = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec split_vars = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec split_values = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec isLeaf = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec parents = arma::zeros<arma::uvec>(MAX_NODE);
  ranges.row(0) = range0.t();
  arma::field<arma::uvec> nodeSampleY(MAX_NODE);
  nodeSampleY(0) = arma::regspace<arma::uvec>(0, n-1);
  arma::field<arma::uvec> nodeSample(K, MAX_NODE);
  for(size_t k = 0; k < K; k++) {
    nodeSample(k,0) = arma::regspace<arma::uvec>(0, mat2Zf(k).n_rows-1);
  }
  size_t ndcount = 0;
  size_t countsp = 0;
  int end = 0;
  while(end == 0) {
    end = split(mat1Z, mat1f, mat2Zf,
                left_childs, right_childs,
                split_vars, split_values, isLeaf,
                parents, fmat, Smat, 
                ranges, nodeSampleY, nodeSample,
                countsp, ndcount, e);
    if(ndcount >= MAX_NODE - 2) {
      isLeaf( arma::find(left_childs == 0) ).ones();
      break;
    }
  }
  arma::uvec nonEmpty = arma::regspace<arma::uvec>(0, ndcount);
  std::shared_ptr<Tree> tr(new Tree(left_childs(nonEmpty),
                                    right_childs(nonEmpty),
                                    split_vars(nonEmpty),
                                    split_values(nonEmpty),
                                    isLeaf(nonEmpty),
                                    parents(nonEmpty)  ));
  if( ndcount + 1 <=  MAX_NODE-1 ) {
    arma::uvec Empty = arma::regspace<arma::uvec>(ndcount+1, MAX_NODE-1);
    fmat.shed_cols(Empty);
    Smat.shed_cols(Empty);
  }
  if(ndcount > 1) {
    fmat.col(0) = fmat.col(1) + fmat.col(2);
    Smat.col(0) = Smat.col(1) + Smat.col(2);
  } else {
    // no need to prune if there is only one node
  }
  return tr;
}



std::shared_ptr<Tree> TreeGrow::grow(const arma::umat& mat1Z,
                                     const arma::mat& mat1f,
                                     const arma::field<arma::umat>& mat2Zf,
                                     const arma::umat& mat1ZVal,
                                     const arma::mat& mat1fVal,
                                     const arma::field<arma::umat>& mat2ZfVal,
                                     const arma::umat& range0,
                                     arma::mat& fmat,
                                     arma::umat& Smat,
                                     arma::mat& fmat2,
                                     arma::umat& Smat2,
                                     const arma::uvec& e) const
{
  int n = mat1Z.n_rows;
  int P = mat1Z.n_cols;
  arma::ucube ranges = arma::zeros<arma::ucube>(MAX_NODE, P, 2);
  arma::uvec left_childs = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec right_childs = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec split_vars = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec split_values = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec isLeaf = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec parents = arma::zeros<arma::uvec>(MAX_NODE);
  ranges.row(0) = range0.t();
  arma::field<arma::uvec> nodeSampleY(MAX_NODE);
  nodeSampleY(0) = arma::regspace<arma::uvec>(0, n-1);
  arma::field<arma::uvec> nodeSample(K,MAX_NODE);
  for(size_t k = 0; k < K; k++) {
    nodeSample(k,0) = arma::regspace<arma::uvec>(0, mat2Zf(k).n_rows-1);
  }
  int nVal = mat1ZVal.n_rows;
  arma::field<arma::uvec> nodeSampleYVal(MAX_NODE);
  nodeSampleYVal(0) = arma::regspace<arma::uvec>(0, nVal-1);
  arma::field<arma::uvec> nodeSampleVal(K,MAX_NODE);
  for(size_t k = 0; k < K; k++) {
    nodeSampleVal(k,0) = arma::regspace<arma::uvec>(0, mat2ZfVal(k).n_rows-1);
  }
  size_t ndcount = 0;
  size_t countsp = 0;
  int end = 0;
  while(end == 0) {
    end = split(mat1Z, mat1f, mat2Zf,mat1ZVal, mat1fVal, mat2ZfVal,
                left_childs, right_childs,
                split_vars, split_values, isLeaf,
                parents, fmat, Smat, fmat2, Smat2, 
                ranges, nodeSampleY, nodeSample, nodeSampleYVal, nodeSampleVal,
                countsp, ndcount, e);
    if(ndcount >= MAX_NODE - 2) {
      isLeaf( arma::find(left_childs == 0) ).ones();
      break;
    }
  }
  arma::uvec nonEmpty = arma::regspace<arma::uvec>(0, ndcount);
  std::shared_ptr<Tree> tr(new Tree(left_childs(nonEmpty),
                                    right_childs(nonEmpty),
                                    split_vars(nonEmpty),
                                    split_values(nonEmpty),
                                    isLeaf(nonEmpty),
                                    parents(nonEmpty)  ));
  if( ndcount + 1 <=  MAX_NODE-1 ) {
    arma::uvec Empty = arma::regspace<arma::uvec>(ndcount+1, MAX_NODE-1);
    fmat.shed_cols(Empty);
    Smat.shed_cols(Empty);
    fmat2.shed_cols(Empty);
    Smat2.shed_cols(Empty);
  }


  if(ndcount > 1)
  {
    fmat.col(0) = fmat.col(1) + fmat.col(2);
    Smat.col(0) = Smat.col(1) + Smat.col(2);
    fmat2.col(0) = fmat2.col(1) + fmat2.col(2);
    Smat2.col(0) = Smat2.col(1) + Smat2.col(2);
  }else{
    // no need to prune if there is only one node
  }
  return tr;
}

arma::vec TreeGrow::prune(arma::vec& beta,
                          const arma::umat& mat1Z,
                          const arma::mat& mat1f,
                          const arma::field<arma::umat>& mat2Zf,
                          const arma::umat& range0,
                          const arma::uvec& e) const
{

  // create folds
  int n = mat1Z.n_rows;
  arma::uvec s = arma::shuffle(arma::regspace<arma::uvec>(0,n-1));
  int quotient = n / NUM_FOLD;
  int remainder = n % NUM_FOLD;
  arma::uvec foldstart = arma::regspace<arma::uvec>(0, quotient, n-1);
  arma::uvec foldend = arma::regspace<arma::uvec>(quotient-1, quotient, n-1);
  arma::uvec::iterator it = foldstart.end();
  arma::uvec::iterator it2 = foldend.end();
  it--;
  it2--;
  while(remainder > 0 ) {
    (*it) = *it + remainder;
    (*it2) = *it2 + remainder;
    it--;
    it2--;
    remainder--;
  }
  // create folds end
  arma::mat iconVal(beta.n_elem, NUM_FOLD);
  for(size_t l = 0; l != NUM_FOLD; l++) {
    // get training/validation set
    arma::uvec valid = s(arma::regspace<arma::uvec>(foldstart(l), foldend(l)));
    arma::uvec trainid;
    if(l == 0) {
      trainid = s(arma::regspace<arma::uvec>(foldend(l) + 1, n - 1));
    }
    else if(l == NUM_FOLD-1) {
      trainid = s(arma::regspace<arma::uvec>(0, foldstart(l) - 1));
    } else {
      arma::uvec trainid1 = arma::regspace<arma::uvec>(0, foldstart(l) - 1);
      arma::uvec trainid2 = arma::regspace<arma::uvec>(foldend(l) + 1, n - 1);
      trainid = s(join_cols(trainid1, trainid2));
    }
    arma::field<arma::umat>  mat2Ztrain(K);
    arma::field<arma::umat>  mat2Zval(K);
    for(size_t k = 0; k != K; k++) {
      arma::umat mat2Zfk = mat2Zf(k);
      arma::ivec idk = arma::conv_to<arma::ivec>::from(trainid + mat2Zfk.n_rows - n);
      arma::ivec idkp = idk( find(idk>=0) );
      arma::uvec idkp2 = arma::conv_to<arma::uvec>::from(idkp);
      mat2Ztrain(k) = mat2Zfk.rows( idkp2  );
      arma::ivec idkv = arma::conv_to<arma::ivec>::from(valid + mat2Zfk.n_rows - n);
      arma::ivec idkpv = idkv( find(idkv>=0) );
      arma::uvec idkpv2 = arma::conv_to<arma::uvec>::from(idkpv);
      mat2Zval(k) = mat2Zfk.rows( idkpv2  );
    }
    // arma::umat mat1Ztrain = mat1Z.rows( trainid );
    // arma::mat mat1ftrain = mat1f.cols( trainid );
    // arma::umat mat1Zval = mat1Z.rows( valid );
    // arma::mat mat1fval = mat1f.cols( valid );
    // arma::uvec etrain = e(trainid);
    // get training/validation set end
    arma::mat fmat = arma::zeros<arma::mat>(K, MAX_NODE);
    arma::umat Smat = arma::zeros<arma::umat>(K, MAX_NODE);
    arma::mat fmat2 = arma::zeros<arma::mat>(K, MAX_NODE);
    arma::umat Smat2 = arma::zeros<arma::umat>(K, MAX_NODE);
    TreeGrow tg(K, spCriterion, MAX_NODE, MIN_NODE1, MIN_SPLIT1);
    std::shared_ptr<Tree> trl = tg.grow(mat1Z.rows( trainid ),
                                        mat1f.cols( trainid ),
                                        mat2Ztrain,
                                        mat1Z.rows( valid ),
                                        mat1f.cols( valid ),
                                        mat2Zval,
                                        range0,
                                        fmat, Smat, fmat2, Smat2, e(trainid));
    const arma::uvec& il = trl->get_isLeaf();
    uint numLeaf = arma::sum(il);
    arma::field<arma::uvec> nodeSetList(numLeaf);
    arma::vec iconAll(numLeaf);
    if (numLeaf > 1) 
      trl->findOptimalSizekSubtree(fmat, Smat, iconAll, nodeSetList, numLeaf);
    arma::vec sizeTree = arma::regspace<arma::vec>(1,iconAll.n_elem);
    arma::uvec isLeafTemp(il.n_elem);
    for(size_t j = 0; j < beta.n_elem; j++) {
      arma::vec iconbetajAll = iconAll -  beta(j) * arma::pow(sizeTree,1);
      uint opt = iconbetajAll.index_max();
      isLeafTemp.zeros();
      isLeafTemp( nodeSetList(opt) ).ones();
      iconVal(j, l) = TreeGrow::get_ICONTest(isLeafTemp, fmat, Smat, fmat2, Smat2);
    }
  }
  return arma::sum(iconVal, 1)/NUM_FOLD;
}



int TreeGrow::split(const arma::umat& mat1Z,
                    const arma::mat& mat1f,
                    const arma::field<arma::umat>& mat2Zf, // dat
                    arma::uvec& left_childs,
                    arma::uvec& right_childs,
                    arma::uvec& split_vars,
                    arma::uvec& split_values,
                    arma::uvec& isLeaf,
                    arma::uvec& parents,
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
    if( spCriterion == 1 ) {
      bestSp = find_split_DICON(nd,
                                mat1Z, mat1f, mat2Zf,
                                ranges, nodeSampleY, nodeSample,
                                fmat, Smat, ndcount, e);
    } else {
      bestSp = find_split_ICON(nd,
                               mat1Z, mat1f, mat2Zf, isLeaf,
                               ranges, nodeSampleY, nodeSample,
                               fmat, Smat, ndcount, e);
    }
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
  if(varsp != -1) {
    split_vars(nd) = varsp;
    split_values(nd) = cutsp;
    arma::uword ndc1 = ndcount + 1;
    arma::uword ndc2 = ndcount + 2;
    left_childs(nd) = ndc1;
    right_childs(nd) = ndc2;
    parents(ndc1) = nd;
    parents(ndc2) = nd;
    arma::uvec nodeSampleYnd = std::move(nodeSampleY(nd));
    arma::uvec zvarspsub = mat1Z( varsp*mat1Z.n_rows + nodeSampleYnd );
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

int TreeGrow::split(const arma::umat& mat1Z,
                    const arma::mat& mat1f,
                    const arma::field<arma::umat>& mat2Zf,
                    const arma::umat& mat1ZVal,
                    const arma::mat& mat1fVal,
                    const arma::field<arma::umat>& mat2ZfVal,// dat
                    arma::uvec& left_childs,
                    arma::uvec& right_childs,
                    arma::uvec& split_vars,
                    arma::uvec& split_values,
                    arma::uvec& isLeaf,
                    arma::uvec& parents,
                    arma::mat& fmat,
                    arma::umat& Smat,
                    arma::mat& fmat2,
                    arma::umat& Smat2,// tree
                    arma::ucube& ranges,
                    arma::field<arma::uvec>& nodeSampleY,
                    arma::field<arma::uvec>& nodeSample,
                    arma::field<arma::uvec>& nodeSampleYVal,
                    arma::field<arma::uvec>& nodeSampleVal,
                    size_t& countsp,
                    size_t& ndcount,
                    const arma::uvec& e) const {
  int end = 0;
  int varsp = -1;
  int cutsp = 0;
  size_t nd = countsp;
  if( spCriterion == 1) {
    while(varsp == -1 && countsp <= ndcount) {
      nd = countsp;
      arma::ivec bestSp(3);
      bestSp = find_split_DICON(nd,
				mat1Z, mat1f, mat2Zf,
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
  }
  if( spCriterion != 1) {
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
  }         
  int n = mat1Z.n_rows;
  if(varsp != -1) {
    split_vars(nd) = varsp;
    split_values(nd) = cutsp;
    arma::uword ndc1 = ndcount + 1;
    arma::uword ndc2 = ndcount + 2;
    left_childs(nd) = ndc1;
    right_childs(nd) = ndc2;
    parents(ndc1) = nd;
    parents(ndc2) = nd;
    arma::uvec nodeSampleYnd = std::move(nodeSampleY(nd));
    arma::uvec zvarspsub = mat1Z( varsp*n + nodeSampleYnd );
    nodeSampleY(ndc1) = nodeSampleYnd( arma::find(zvarspsub <= cutsp) );
    nodeSampleY(ndc2) = nodeSampleYnd( arma::find(zvarspsub > cutsp) );
    arma::uvec nodeSampleYndVal = std::move(nodeSampleYVal(nd));
    arma::uvec zvarspsubVal = mat1ZVal( varsp*mat1ZVal.n_rows + nodeSampleYndVal );
    nodeSampleYVal(ndc1) = nodeSampleYndVal( arma::find(zvarspsubVal <=cutsp) );
    nodeSampleYVal(ndc2) = nodeSampleYndVal( arma::find(zvarspsubVal >cutsp) );
    fmat2.col(ndc1) = arma::sum( mat1fVal.cols( nodeSampleYVal(ndc1) ), 1);
    fmat2.col(ndc2) = arma::sum( mat1fVal.cols( nodeSampleYVal(ndc2) ), 1);
    for(size_t k = 0; k < K; k++) {
      arma::uvec nodeSampleknd = std::move(nodeSample(k,nd));
      arma::uvec zvarspsub = mat2Zf(k)( varsp*mat2Zf(k).n_rows + nodeSampleknd );
      nodeSample(k, ndc1) = nodeSampleknd( arma::find(zvarspsub <= cutsp) );
      nodeSample(k, ndc2) = nodeSampleknd( arma::find(zvarspsub > cutsp) );
    }
    if(nodeSample(0, ndc1).size() < MIN_SPLIT1) isLeaf(ndc1) = 1;
    if(nodeSample(0, ndc2).size() < MIN_SPLIT1) isLeaf(ndc2) = 1;
    for(size_t k = 0; k < K; k++) {
      arma::uvec nodeSamplekndVal = std::move(nodeSampleVal(k,nd));
      arma::uvec zvarspsubVal = mat2ZfVal(k)( varsp*mat2ZfVal(k).n_rows + nodeSamplekndVal );
      nodeSampleVal(k, ndc1) = nodeSamplekndVal( arma::find(zvarspsubVal<=cutsp) );
      nodeSampleVal(k, ndc2) = nodeSamplekndVal( arma::find(zvarspsubVal>cutsp) );
      Smat2(k,ndc1) = nodeSampleVal(k, ndc1).n_elem;
      Smat2(k,ndc2) = nodeSampleVal(k, ndc2).n_elem;
    }
    ranges.row(ndc1) = ranges.row(nd);
    ranges.row(ndc2) = ranges.row(nd);
    ranges(ndc2,varsp,0) = cutsp+1;
    ranges(ndc1,varsp,1) = cutsp;
    ndcount += 2;

    while(countsp <= ndcount) {
      countsp++;
      if(isLeaf(countsp) == 0)
      {break;}
    }
  } else {
    end = 1;
  }
  return end;
}



arma::ivec TreeGrow::find_split_DICON(size_t nd,
                                      const arma::umat& mat1Z,
                                      const arma::mat& mat1f,
                                      const arma::field<arma::umat>& mat2Zf, // dat
                                      const arma::ucube& ranges,
                                      const arma::field<arma::uvec>& nodeSampleY,
                                      const arma::field<arma::uvec>& nodeSample,
                                      arma::mat& fmat,
                                      arma::umat& Smat,
                                      size_t ndcount,
                                      const arma::uvec& e) const {
  int P = mat1Z.n_cols;
  int n = mat1Z.n_rows;
  int varsp = -1;
  int cutsp = 0;
  double dICONmax = 0;
  double dICONTemp = 0;
  for(int p = 0; p < P; p++) {
    // arma::uvec zp = mat1Z.col(p);
    // arma::uvec zpsub = zp(nodeSampleY(nd));
    // arma::uvec ind = arma::regspace<arma::uvec>(0,n); // should be n-1
    // arma::uvec indsub = ind(nodeSampleY(nd));
    // arma::uvec indY = indsub(sort_index(zpsub));
    arma::uvec indY = nodeSampleY(nd)( sort_index( mat1Z(p*n + nodeSampleY(nd)) ));
    //indY.print();
    arma::field<arma::uvec> indp(K);
    arma::uvec SRSum = arma::zeros<arma::uvec>(K);
    for(size_t k = 0; k < K; k++) {
      // arma::umat mat2zp = mat2Zf(k);
      // arma::uvec zp = mat2zp.col(p);
      // arma::uvec zpsub = zp(nodeSample(k,nd));
      // arma::uvec ind = arma::regspace<arma::uvec>(0,zp.size() );
      // arma::uvec indsub = ind(nodeSample(k,nd));
      // indp(k) = indsub(sort_index(zpsub));
      arma::uvec zpsub = mat2Zf(k)( p*mat2Zf(k).n_rows + nodeSample(k,nd) );
      indp(k) = nodeSample(k,nd)(sort_index(zpsub));
      SRSum(k) = zpsub.size();
    }
    arma::vec fLSum = arma::zeros<arma::vec>(K);
    arma::uvec SLSum = arma::zeros<arma::uvec>(K);
    arma::vec fRSum = arma::sum(mat1f.cols(indY),1);
    int j = 0;
    arma::uvec jv = arma::zeros<arma::uvec>(K);
    int nj = indY.size();
    arma::uvec rangeCut = arma::regspace<arma::uvec>(ranges(nd, p, 0), ranges(nd, p, 1));
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
        dICONTemp = arma::sum(abs(fLSum%SRSum - fRSum%SLSum) );
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

arma::ivec TreeGrow::find_split_ICON(size_t nd,
                                     const arma::umat& mat1Z,
                                     const arma::mat& mat1f,
                                     const arma::field<arma::umat>& mat2Zf, // dat
                                     const arma::uvec& isLeaf,
                                     arma::ucube& ranges,
                                     arma::field<arma::uvec>& nodeSampleY,
                                     arma::field<arma::uvec>& nodeSample,
                                     arma::mat& fmat,
                                     arma::umat& Smat,
                                     size_t ndcount,
                                     const arma::uvec& e) const {
  int P = mat1Z.n_cols;
  int n = mat1Z.n_rows;
  int varsp = -1;
  int cutsp = 0;
  double dICONmax = 0;
  double dICONTemp = 0;
  arma::mat fmatTerm = fmat.cols(arma::find(isLeaf == 1));
  arma::umat SmatTerm = Smat.cols(arma::find(isLeaf == 1));     
  for(int p = 0; p < P; p++) {
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
    arma::vec fRSum = (sum(mat1f.cols(indY),1));
    int j = 0;
    arma::uvec jv = arma::zeros<arma::uvec>(K);
    int nj = indY.size();
    int nel = 0;
    // int nelr = arma::sum( e(indY) );
    arma::uvec rangeCut = arma::regspace<arma::uvec>(ranges(nd, p, 0), ranges(nd, p, 1));
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
      //if( (SLSum.min() < MIN_NODE1 || SRSum.min() < MIN_NODE1) && ( nel < MIN_NODE2 || nelr-nel < MIN_NODE2) )
      if(  (SLSum(0) < MIN_NODE1 || SRSum(0) < MIN_NODE1) ) {
        dICONTemp = 0;
      } else {
        dICONTemp = arma::sum(abs(fLSum%SRSum - fRSum%SLSum))/2;
        for(size_t i = 0; i < fmatTerm.n_cols; i++) {
          if( i != nd ) {
            for(size_t k = 0; k < K; k++) {
              //  DOUBLE CHECK !!!
              if( (fLSum(k) + fRSum(k)) * SmatTerm(k,i) <= fmatTerm(k,i) * (SLSum(k) + SRSum(k)) ) {
                if( fmatTerm(k,i)*SLSum(k) <  fLSum(k)*SmatTerm(k,i) )
                {// lambdaR < lambdand <= lambdai < lambdaL
                  dICONTemp = dICONTemp +  fLSum(k)*SmatTerm(k,i) - fmatTerm(k,i)*SLSum(k) ;
                }
                if( fmatTerm(k,i)*SRSum(k) <  fRSum(k)*SmatTerm(k,i) )
                {// lambdaL < lambdand <= lambdai < lambdaR
                  dICONTemp = dICONTemp + fRSum(k)*SmatTerm(k,i) - fmatTerm(k,i)*SRSum(k) ;
                }
              } else {
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
