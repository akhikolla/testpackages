#include "TreePrediction.h"
#include "globals.h"

void TreePrediction::transformZ(const arma::mat& z,
                                arma::umat& z2,
                                const arma::mat& dat,
                                const arma::uvec& e,
                                const arma::vec& breaks,
                                const arma::uvec& disc)
{
  int P = z.n_rows;
  int n = e.n_elem;
  arma::vec::const_iterator bb = breaks.begin();
  arma::vec::const_iterator be = breaks.end();
  for(int p = 0; p < P; p++) {
    if(disc(p) == 0) {
      arma::uvec ind = arma::cumsum(arma::regspace<arma::uvec>(0,n-1));
      arma::vec zp = dat.col(p);
      // BE CAREFUL
      int j = 0;
      double z2ecdf;
      for(int i = 0; i < n; i++) {
        if(e(i) == 1) {
          arma::vec zpref = arma::sort(zp.elem( ind ));
          z2ecdf = (std::lower_bound(zpref.begin(), zpref.end(), z(p,i) )- zpref.begin())/(n-i+0.0)  ;
          z2(p,j) = std::distance(bb, std::upper_bound(bb, be, z2ecdf ) ) + 1;
          j++;
        }
        ind = ind + 1;
        if(ind.n_elem > 0){ ind.shed_row(0);}
      }
    } else {
      z2.row(p) = arma::conv_to<arma::urowvec>::from(z.row(p));
    }
  }
}

void TreePrediction::transformZH(const arma::mat& z,  // z on tg
                                 const arma::vec& tg, //grid of time point
                                 arma::umat& z2,
                                 const arma::mat& dat, // dat is predictor matrix
                                 const arma::vec& y,
                                 const arma::uvec& e,
                                 const arma::vec& breaks,
                                 const arma::uvec& disc)
{
  int P = z.n_rows;
  int n = e.n_elem;
  int NG = tg.n_elem;
  arma::vec::const_iterator bb = breaks.begin();
  arma::vec::const_iterator be = breaks.end();
  for(int p = 0; p < P; p++) {
    if(disc(p) == 0) {
      arma::uvec ind = arma::cumsum(arma::regspace<arma::uvec>(0,n-1));
      arma::vec zp = dat.col(p);
      // BE CAREFUL
      int j = 0;
      double z2ecdf;
      for(int i = 0; i < n-1; i++) {
        if( y(i) <= tg(j) && tg(j) < y(i+1) ) {
          arma::vec zpref = arma::sort(zp.elem( ind ));
          z2ecdf = (std::lower_bound(zpref.begin(), zpref.end(), z(p,j) )- zpref.begin())/(n-i+0.0)  ;
          z2(p,j) = std::distance(bb, std::upper_bound(bb, be, z2ecdf ) ) + 1;
          j++;
          if(j == NG) break;
        }
        ind = ind + 1;
        if(ind.n_elem > 0){ ind.shed_row(0);}
      }
    } else {
      z2.row(p) = arma::conv_to<arma::urowvec>::from(z.row(p));
    }
  }
}


TreePrediction::TreePrediction(const arma::umat& zy,
                               const arma::field<arma::umat>& zt,
                               const arma::uvec& vars,
			       const arma::uvec& values,
			       const arma::uvec& lcs,
			       const arma::uvec& rcs,
			       const arma::uvec& il)
{
  arma::uvec ndy(zy.n_cols);
  for(size_t it = 0; it != zy.n_cols; it++) {
    arma::uvec zyi = (zy.col(it));
    int isl = 0;
    int varsp = 0;
    size_t cutsp = 0;
    int k = 0;
    while(isl == 0) {
      varsp = vars(k);
      cutsp = values(k);
      if(zyi(varsp) > cutsp) {
        k = rcs(k);
      } else {
        k = lcs(k);
      }
      isl = il(k);
    }
    ndy(it) = k;
  }
  size_t nT = zt.size();
  //arma::field<arma::uvec> ndst(nT);
  int nNd = arma::accu(il);
  //Rcpp::Rcout << arma::find(il == 1);
  arma::uvec tnd = arma::regspace<arma::uvec>(0, il.n_elem-1);
  //Rcpp::Rcout << MAX_NODE;
  arma::uvec tnd2 = tnd(arma::find(il == 1));
  arma::umat ndsz = arma::zeros<arma::umat>(nNd, nT);
  arma::uvec tnd3 = arma::zeros<arma::uvec>(il.n_elem);
  tnd3( tnd2 ) = arma::regspace<arma::uvec>(0, nNd-1);
  for(arma::uword c=0; c < nT; c++) {
    arma::umat m = zt(c);
    int it_end = m.n_cols;
    arma::uvec ndt(it_end);
    for(int it = 0; it != it_end; it++) {
      arma::uvec zti = (m.col(it));
      int isl = 0;
      int varsp = 0;
      size_t cutsp = 0;
      int k = 0;
      while(isl == 0) {
        varsp = vars(k);
        cutsp = values(k);
        if(zti(varsp) > cutsp) {
          k = rcs(k);
        } else {
          k = lcs(k);
        }
        isl = il(k);
      }
      ndt(it) = k;
      //Rcpp::Rcout << k;
      ndsz(tnd3(k),c)++;
    }
    //ndst(c) = ndt;
  }
  this->nodeLabel = ndy;
  this->nodeSize = ndsz;
  this->tnd3 = tnd3;
}


arma::vec TreePrediction::getSurvival(const arma::umat& zt2,
                                      const arma::vec& y,
                                      const arma::uvec& e,
                                      const arma::uvec& vars,
				      const arma::uvec& values,
				      const arma::uvec& lcs,
				      const arma::uvec& rcs,
				      const arma::uvec& il)
{
  arma::vec y2 = y( arma::find(e == 1));
  int NE = y2.n_elem;
  arma::vec w2 = arma::zeros<arma::vec>(NE);
  for(int i = 0; i < NE; i++) {
    arma::uvec zY = zt2.col(i);
    int isl = 0;
    int varsp = 0;
    size_t cutsp = 0;
    int k = 0;
    while(isl == 0) {
      varsp = vars(k);
      cutsp = values(k);
      if(zY(varsp) > cutsp) {
        k = rcs(k);
      } else {
        k = lcs(k);
      }
      isl = il(k);
    }
    size_t nd = k;
    if(nodeLabel(i) == nd) {
      w2(i) = w2(i) + 1.0/nodeSize(tnd3(nd),i);
    }
  }
  arma::vec w = arma::zeros<arma::vec>(y.n_elem);
  w( arma::find(e == 1)) = w2;
  return exp(-cumsum(w));
}

arma::vec TreePrediction::getSurvival(const arma::umat& zt2,
                                      const arma::vec& y,
                                      const arma::uvec& e,
                                      const arma::umat& nodeSize0,
                                      const arma::uvec& nodeLabel0,
                                      const arma::uvec& tnd30,
                                      const arma::umat& treeMat)
{
  arma::vec y2 = y( arma::find(e == 1));
  arma::uvec vars = treeMat.col(0);
  arma::uvec values = treeMat.col(1);
  arma::uvec lcs = treeMat.col(2);
  arma::uvec rcs = treeMat.col(3);
  arma::uvec il = treeMat.col(4);
  int NE = y2.n_elem;
  arma::vec w2 = arma::zeros<arma::vec>(NE);
  for(int i = 0; i < NE; i++) {
    arma::uvec zY = zt2.col(i);
    int isl = 0;
    int varsp = 0;
    size_t cutsp = 0;
    size_t k = 0;
    while(isl == 0) {
      varsp = vars(k);
      cutsp = values(k);
      if(zY(varsp) > cutsp) {
        k = rcs(k);
      } else {
        k = lcs(k);
      }
      isl = il(k);
    }
    // size_t nd = k;
    if(nodeLabel0(i) == k) {
      w2(i) = w2(i) + 1.0/nodeSize0(tnd30(k),i);
    }
  }
  arma::vec w = arma::zeros<arma::vec>(y.n_elem);
  w( arma::find(e == 1)) = w2;
  return exp(-cumsum(w));
}


arma::vec TreePrediction::getHazard(const arma::umat& ztvec,
                                    const arma::vec& tg,
                                    const arma::vec& y,
                                    const arma::uvec& e,
                                    const arma::mat& fy2,
                                    const double h,
                                    const arma::umat& nodeSize0,
                                    const arma::uvec& nodeLabel0,
                                    const arma::uvec& tnd30,
                                    const arma::umat& treeMat)
{
  arma::vec y2 = y( arma::find(e == 1));
  arma::uvec vars = treeMat.col(0);
  arma::uvec values = treeMat.col(1);
  arma::uvec lcs = treeMat.col(2);
  arma::uvec rcs = treeMat.col(3);
  arma::uvec il = treeMat.col(4);
  int NG = tg.n_elem;
  int NE = y2.n_elem;
  arma::vec hz = arma::zeros<arma::vec>(NG);
  for(int j = 0; j < NG; j++) {
    arma::uvec zj = ztvec.col(j);
    int isl = 0;
    int varsp = 0;
    size_t cutsp = 0;
    size_t k = 0;
    while(isl == 0) {
      varsp = vars(k);
      cutsp = values(k);
      if(zj(varsp) > cutsp) {
        k = rcs(k);
      } else {
        k = lcs(k);
      }
      isl = il(k);
    }
    arma::vec v0 = arma::zeros<arma::vec>(NE);
    arma::vec v1 = arma::zeros<arma::vec>(NE);
    for(int i = 0; i < NE; i++) {
      if( y2(i) >= tg(j) - h && y2(i) <= tg(j) + h && nodeLabel0(i) == k) {
        //hz(j) = hz(j) + fy2(i,j)/nodeSize0(tnd30(k),i);
        v0(i)++;
      }
      v1(i) += nodeSize0(tnd30(k),i);
    }
    v1(arma::find(v1 == 0)).ones();
    hz(j) = arma::sum( fy2.col(j) %  v0 / v1 );
  }
  return hz;
}

TreePrediction::TreePrediction(const Data2* dat2,
                               const arma::uvec& vars,
			       const arma::uvec& values,
			       const arma::uvec& lcs,
			       const arma::uvec& rcs,
			       const arma::uvec& il)
{
  arma::umat zy = dat2->get_zy();
  arma::field<arma::umat> zt = dat2->get_zt();
  arma::uvec ndy(zy.n_cols);
  for(size_t it = 0; it != zy.n_cols; it++) {
    arma::uvec zyi = zy.col(it);
    int isl = 0;
    int varsp = 0;
    size_t cutsp = 0;
    int k = 0;
    while(isl == 0) {
      varsp = vars(k);
      cutsp = values(k);
      if(zyi(varsp) > cutsp) {
        k = rcs(k);
      } else {
        k = lcs(k);
      }
      isl = il(k);
    }
    ndy(it) = k;
  }
  size_t nT = zt.size();
  //arma::field<arma::uvec> ndst(nT);
  int nNd = arma::accu(il);
  arma::uvec tnd = arma::regspace<arma::uvec>(0, il.n_elem-1);
  arma::uvec tnd2 = tnd(arma::find(il == 1));
  arma::umat ndsz = arma::zeros<arma::umat>(nNd, nT);
  arma::uvec tnd3 = arma::zeros<arma::uvec>(il.n_elem);
  tnd3.elem( tnd2 ) = arma::regspace<arma::uvec>(0, nNd-1);
  for(arma::uword c=0; c < nT; c++) {
    arma::umat m = zt(c);
    int it_end = m.n_cols;
    arma::uvec ndt(it_end);
    for(int it = 0; it != it_end; it++) {
      arma::uvec zti = (m.col(it));
      int isl = 0;
      int varsp = 0;
      size_t cutsp = 0;
      int k = 0;
      while(isl == 0) {
        varsp = vars(k);
        cutsp = values(k);
        if(zti(varsp) > cutsp) {
          k = rcs(k);
        } else {
          k = lcs(k);
        }
        isl = il(k);
      }
      ndt(it) = k;
      ndsz(tnd3(k),c)++;
    }
  }
  this->nodeLabel = ndy;
  this->nodeSize = ndsz;
  this->tnd3 = tnd3;
}
