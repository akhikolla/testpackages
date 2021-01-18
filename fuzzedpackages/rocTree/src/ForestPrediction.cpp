#include "ForestPrediction.h"
#include "globals.h"
#include <algorithm>


void ForestPrediction::transformZ(const arma::mat& z,
                                  arma::umat& z2,
                                  const arma::mat& matX,
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
      arma::vec zp = matX.col(p);
      // BE CAREFUL
      int j = 0;
      for(int i = 0; i < n; i++) {
        if(e(i) == 1) {
          arma::vec zpref = arma::sort(zp.elem( ind ));
          double z2ecdf = (std::lower_bound(zpref.begin(), zpref.end(), z(p,i) ) - zpref.begin())/(n-i+0.0);
          z2(p,j) = std::distance(bb, std::upper_bound(bb, be, z2ecdf ) ) + 1 ;
          j++;
        }
        ind = ind + 1;
        if(ind.n_elem > 0) ind.shed_row(0);
      }
    } else {
      z2.row(p) = arma::conv_to<arma::urowvec>::from(z.row(p));
    }
  }
}

void ForestPrediction::transformZH(const arma::mat& z,  // z on tg
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
        if(ind.n_elem > 0) ind.shed_row(0);
      }
    } else {
      z2.row(p) = arma::conv_to<arma::urowvec>::from(z.row(p));
    }
  }
}

// The new version
arma::vec ForestPrediction::getSurvival(const arma::umat& zt2,
                                        const arma::vec& y,
                                        const arma::uvec& e,
                                        const arma::field<arma::uvec>&& nodeSizeB0,
                                        const arma::umat&& nodeLabelB0,
                                        const arma::field<arma::uvec>&& tnd3B0,
                                        const arma::umat&& ids,
                                        const arma::field<arma::umat>&& trees)
{
  arma::uvec e1 = arma::find(e == 1);
  arma::vec y2 = y( e1 );
  arma::uvec idAll = arma::regspace<arma::uvec>(0, y.n_elem-1);
  arma::uvec idY2 = idAll( e1 );
  uint NE = y2.n_elem;
  uint NUM_TREE = trees.size();
  arma::vec w0 = arma::zeros<arma::vec>(NE);
  arma::vec w1 = arma::zeros<arma::vec>(NE);
  for(size_t b = 0; b != NUM_TREE; b++) {
    // take the b th tree
    arma::umat tt = (trees(b));
    arma::uvec vars =  tt.col(0);
    arma::uvec values = tt.col(1);
    arma::uvec lcs = tt.col(2);
    arma::uvec rcs = tt.col(3);
    arma::uvec il = tt.col(4);
    arma::uvec tnd3 = (tnd3B0(b));
    arma::uvec zY;
    // take the ids with uncensored survival times
    // IDs have been sorted when sampling is done; need to sort if not sorted before
    arma::uvec idb = ids.col(b);
    //arma::uvec idb = arma::sort(ids.col(b));
    arma::uvec idbe = idb( arma::find( e(idb) == 1 )  );
    int j = 0;
    int nidbe = idbe.n_elem;
    for(size_t i = 0; i < NE; i++) {
      int count = 0;
      while(j < nidbe) {
        if( idbe(j) == idY2(i) ) {
          count++;
          j++;
        } else {
          break;
        }
      }
      zY = zt2.col(i);
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
      arma::uvec nbi = nodeSizeB0(b + NUM_TREE * i);
      if(nodeLabelB0(b,i) == k) {
	w0(i) += count ;
      }
      w1(i) += nbi(tnd3(k));
    }   
  }
  arma::vec w = arma::zeros<arma::vec>(y.n_elem);
  w(e1) = w0/w1;
  w.replace(arma::datum::nan, 0);
  return exp(-cumsum(w));
}



arma::vec ForestPrediction::getHazard(const arma::umat& ztvec,
                                      const arma::vec& tg,
                                      const arma::vec& y,
                                      const arma::uvec& e,
                                      const arma::mat& Kmat,
                                      const double h,
                                      const arma::field<arma::uvec>&& nodeSizeB0,
                                      const arma::umat&& nodeLabelB0,
                                      const arma::field<arma::uvec>&& tnd3B0,
                                      const arma::umat&& ids,
                                      const arma::field<arma::umat>&& trees)
{
  arma::uvec e1 = arma::find(e == 1);
  arma::vec y2 = y( e1 );
  arma::uvec idAll = arma::regspace<arma::uvec>(0, y.n_elem-1);
  arma::uvec idY2 = idAll( e1 );
  uint NE = y2.n_elem;
  uint NUM_TREE = trees.size();
  uint NG = tg.n_elem;
  arma::mat v0 = arma::zeros<arma::mat>(NE, NG);
  arma::mat v1 = arma::zeros<arma::mat>(NE, NG);
  for(size_t b = 0; b != NUM_TREE; b++) {
    // take the b th tree
    arma::umat tt = (trees(b));
    arma::uvec vars =  tt.col(0);
    arma::uvec values = tt.col(1);
    arma::uvec lcs = tt.col(2);
    arma::uvec rcs = tt.col(3);
    arma::uvec il = tt.col(4);
    arma::uvec tnd3 = tnd3B0(b);
    // take the ids with uncensored survival times
    arma::uvec idb = ids.col(b); // IDs have been sorted when sampling is done; need to sort if not sorted before
    //arma::uvec idb = arma::sort(ids.col(b));
    arma::uvec idbe = idb( arma::find( e(idb) == 1 )  );
    int nidbe = idbe.n_elem;
    for(size_t l = 0; l < NG; l++) {
      arma::uvec zl = ztvec.col(l);
      int isl = 0;
      int varsp = 0;
      size_t cutsp = 0;
      size_t k = 0;
      while(isl == 0) {
        varsp = vars(k);
        cutsp = values(k);
        if(zl(varsp) > cutsp) {
          k = rcs(k);
        } else {
          k = lcs(k);
        }
        isl = il(k);
      }
      int j = 0;
      for(size_t i = 0; i < NE; i++) {
        int count = 0;
        while(j < nidbe) {
          if( idbe(j) == idY2(i) ) {
            count++;
            j++;
          } else {
            break;
          }
        }
	if( y2(i) >= tg(l) - h && y2(i) <= tg(l) + h && nodeLabelB0(b,i) == k) {
          v0(i,l) += count;
        }
        arma::uvec nbi = nodeSizeB0(b+NUM_TREE*i);
        v1(i,l) += nbi(tnd3(k));
      }
    }
  }
  arma::vec hz = arma::zeros<arma::vec>(NG);
  v1(arma::find(v1 == 0)).ones();
  for(size_t l = 0; l < NG; l++) {
    hz(l) = arma::sum(Kmat.col(l) % v0.col(l) / v1.col(l)) ;
  }
  return hz;
}


// The old version
arma::vec ForestPrediction::getSurvival2(const arma::umat& zt2,
                                        const arma::vec& y,
                                        const arma::uvec& e,
                                        const arma::field<arma::uvec>&& nodeSizeB0,
                                        const arma::umat&& nodeLabelB0,
                                        const arma::field<arma::uvec>&& tnd3B0,
                                        const arma::umat&& ids,
                                        const arma::field<arma::umat>&& trees)
{
  arma::uvec e1 = arma::find(e == 1);
  arma::vec y2 = y( e1 );
  arma::uvec idAll = arma::regspace<arma::uvec>(0, y.n_elem-1);
  arma::uvec idY2 = idAll( e1 );
  uint NE = y2.n_elem;
  uint NUM_TREE = trees.size();
  arma::vec w2 = arma::zeros<arma::vec>(NE);
  for(size_t b = 0; b != NUM_TREE; b++) {
    // take the b th tree
    arma::umat tt = (trees(b));
    arma::uvec vars =  tt.col(0);
    arma::uvec values = tt.col(1);
    arma::uvec lcs = tt.col(2);
    arma::uvec rcs = tt.col(3);
    arma::uvec il = tt.col(4);
    arma::uvec tnd3 = (tnd3B0(b));
    arma::uvec zY;
    for(size_t i = 0; i < NE; i++) {
      // need to add this back if using sample-spliting
      //if(std::find(ids.begin_col(b), ids.end_col(b), idY2(i) ) != ids.end_col(b))
      zY = zt2.col(i);
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
      arma::uvec nbi = nodeSizeB0(b+NUM_TREE*i);
      double den = nbi(tnd3(k));
      if(nodeLabelB0(b,i) == k) {
	w2(i) += (1.0/NUM_TREE)/den;
      }
    }
  }
  arma::vec w = arma::zeros<arma::vec>(y.n_elem);
  w( e1 ) = w2;
  //return w;
  return exp(-cumsum(w));
}

ForestPrediction::ForestPrediction(const arma::umat& zy,
                                   const arma::field<arma::umat>& zt,
                                   const arma::umat& ids,
                                   const std::vector<std::shared_ptr<Tree> >& trees,
                                   arma::uword n)

{
  size_t nT = zt.size();
  uint NUM_TREE = trees.size();
  arma::umat ndy(NUM_TREE, nT);
  arma::field<arma::uvec> ndsz(NUM_TREE, nT);
  arma::field<arma::uvec> tnd3B(NUM_TREE);
  std::vector<std::shared_ptr<Tree> >::const_iterator it;
  int i = 0;
  for(it = trees.begin(); it != trees.end(); it++, i++) {
    // take the i th tree
    std::shared_ptr<Tree> tt = *it;
    arma::uvec vars =  tt->get_split_vars();
    arma::uvec lcs = tt->get_left_childs();
    arma::uvec rcs = tt->get_right_childs();
    arma::uvec values = tt->get_split_values();
    arma::uvec il = tt->get_isLeaf();
    int nNd = arma::accu(il);
    //arma::uvec tnd = arma::regspace<arma::uvec>(0, il.n_elem-1);
    arma::uvec tnd2 = (arma::find(il == 1));
    arma::uvec tnd3 = arma::zeros<arma::uvec>(il.n_elem);
    tnd3.elem( tnd2 ) = arma::regspace<arma::uvec>(0, nNd-1);
    for(size_t j = 0; j != nT; j++) {
      arma::uvec zyj = zy.col( j );
      int isl = 0;
      int varsp = 0;
      size_t cutsp = 0;
      size_t k = 0;
      while(isl == 0) {
        varsp = vars(k);
        cutsp = values(k);
        if(zyj(varsp) > cutsp) {
          k = rcs(k);
        } else {
          k = lcs(k);
        }
        isl = il(k);
      }
      ndy(i,j) = k;
    }
    arma::uvec idi = ids.col(i);
    for(arma::uword c=0; c < nT; c++) {
      arma::uvec ndszic = arma::zeros<arma::uvec>(nNd);
      arma::umat m = zt(c);
      //arma::umat m(zt(c).memptr(), zt(c).n_rows, zt(c).n_cols,false);
      //arma::ivec idc = arma::conv_to<arma::ivec>::from(idi + m.n_cols - n);
      //arma::ivec idcp = idc( find(idc>=0) );
      arma::uvec idcp = idi( arma::find(idi >= n-m.n_cols) ) + m.n_cols - n;
      int j_end = idcp.n_elem;
      if(j_end > 0) {
        for(int j = 0; j != j_end; j++) {
          arma::uvec zti = m.col( idcp(j) );
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
          //Rcpp::Rcout << k;
          ndszic(tnd3(k))++;
        }
      }
      //Rcpp::Rcout << ndszic;
      ndsz(i,c) = ndszic;
    }
    tnd3B(i) = tnd3;
  }
  this->nodeLabelB = ndy;
  this->nodeSizeB = ndsz;
  this->tnd3B = tnd3B;
}

// arma::vec ForestPrediction::getSurvival(const arma::umat& zt2,
//                                         const Data2* dat2,
//                                         const arma::umat& ids,
//                                         const std::vector<std::shared_ptr<Tree>>& trees)
// {
//   arma::vec y = dat2->get_Y();
//   arma::uvec e = dat2->get_e();
//   arma::vec y2 = y( arma::find(e == 1));
//   arma::uvec idAll = arma::regspace<arma::uvec>(0, y.n_elem-1);
//   arma::uvec idY2 = idAll( arma::find(e == 1));
//   uint NE = y2.n_elem;
//   uint NUM_TREE = trees.size();
//   arma::vec w2 = arma::zeros<arma::vec>(NE);
//
//
//   std::vector<std::shared_ptr<Tree>>::const_iterator it;
//   int b = 0;
//   for(it = trees.begin(); it != trees.end(); it++)
//   {
//     // take the i th tree
//     std::shared_ptr<Tree> tt = *it;
//     const arma::uvec& vars =  tt->get_split_vars();
//     const arma::uvec& lcs = tt->get_left_childs();
//     const arma::uvec& rcs = tt->get_right_childs();
//     const arma::uvec& values = tt->get_split_values();
//     const arma::uvec& il = tt->get_isLeaf();
//
//     arma::uvec tnd3 = tnd3B(b);
//
//     for(int i = 0; i < NE; i++)
//     {
//       arma::uvec zY = (zt2.col(i));
//
//       int isl = 0;
//       int varsp = 0;
//       int cutsp = 0;
//       int k = 0;
//
//       while(isl == 0)
//       {
//         varsp = vars(k);
//         cutsp = values(k);
//         if(zY(varsp) > cutsp)
//         {
//           k = rcs(k);
//         }else{
//           k = lcs(k);
//         }
//         isl = il(k);
//       }
//
//       int nd = k;
//       arma::uvec nbi = nodeSizeB(b,i);
//       arma::uword den = nbi(tnd3(nd));
//       if(nodeLabelB(b,i) == nd && den > 0 && arma::sum(idY2(i) == ids.col(b))>0 )//&&
//       {
//         w2(i) = w2(i) + (1.0/NUM_TREE)/den;
//       }
//     }
//     b++;
//   }
//
//
//   arma::vec w = arma::zeros<arma::vec>(y.n_elem);
//   w( arma::find(e == 1)) = w2;
//
//   //return w;
//   return exp(-cumsum(w));
// }
//
//
//
// ForestPrediction::ForestPrediction(const Data2* dat2,
//                                    const arma::umat& ids,
//                                    const std::vector<std::shared_ptr<Tree>>& trees,
//                                    arma::uword n)
//
// {
//   arma::umat zy = dat2->get_zy();
//   arma::field<arma::umat> zt = dat2->get_zt();
//
//   int nT = zt.size();
//   uint NUM_TREE = trees.size();
//   arma::umat ndy(NUM_TREE, nT);
//   arma::field<arma::uvec> ndsz(NUM_TREE, nT);
//   arma::field<arma::uvec> tnd3B(NUM_TREE);
//
//   std::vector<std::shared_ptr<Tree>>::const_iterator it;
//   int i = 0;
//   for(it = trees.begin(); it != trees.end(); it++)
//   {
//     // take the i th tree
//     std::shared_ptr<Tree> tt = *it;
//     arma::uvec vars =  tt->get_split_vars();
//     arma::uvec lcs = tt->get_left_childs();
//     arma::uvec rcs = tt->get_right_childs();
//     arma::uvec values = tt->get_split_values();
//     arma::uvec il = tt->get_isLeaf();
//     int nNd = arma::accu(il);
//     arma::uvec tnd = arma::regspace<arma::uvec>(0, il.n_elem-1);
//     arma::uvec tnd2 = tnd(arma::find(il == 1));
//     arma::uvec tnd3 = arma::zeros<arma::uvec>(il.n_elem);
//     tnd3.elem( tnd2 ) = arma::regspace<arma::uvec>(0, nNd-1);
//
//     for(int j = 0; j != nT; j++)
//     {
//       arma::uvec zyj = (zy.col( j  ));
//
//       int isl = 0;
//       int varsp = 0;
//       int cutsp = 0;
//       int k = 0;
//       while(isl == 0)
//       {
//         varsp = vars(k);
//         cutsp = values(k);
//         if(zyj(varsp) > cutsp)
//         {
//           k = rcs(k);
//         }else{
//           k = lcs(k);
//         }
//         isl = il(k);
//       }
//       ndy(i,j) = k;
//     }
//
//     arma::uvec idi = ids.col(i);
//
//     for(arma::uword c=0; c < nT; c++)
//     {
//       arma::uvec ndszic = arma::zeros<arma::uvec>(nNd);
//       arma::umat m = zt(c);
//       arma::ivec idc = arma::conv_to<arma::ivec>::from(idi + m.n_cols - n);
//       arma::ivec idcp = idc( find(idc>=0) );
//       int j_end = idcp.n_elem;
//
//
//       if(j_end > 0)
//       {
//         for(int j = 0; j != j_end; j++)
//         {
//           arma::uvec zti = (m.col( idcp(j) ));
//
//           int isl = 0;
//           int varsp = 0;
//           int cutsp = 0;
//           int k = 0;
//           while(isl == 0)
//           {
//             varsp = vars(k);
//             cutsp = values(k);
//             if(zti(varsp) > cutsp)
//             {
//               k = rcs(k);
//             }else{
//               k = lcs(k);
//             }
//             isl = il(k);
//           }
//
//           //Rcpp::Rcout << k;
//           ndszic(tnd3(k))++;
//         }
//       }
//       //Rcpp::Rcout << ndszic;
//       ndsz(i,c) = ndszic;
//     }
//
//
//     tnd3B(i) = tnd3;
//
//     i++;
//
//   }
//
//
//
//   this->nodeLabelB = ndy;
//   this->nodeSizeB = ndsz;
//   this->tnd3B = tnd3B;
//
// }

// Time-independent way
// void ForestPrediction::transformZ0(const arma::mat& z,
//                                    arma::umat& z2,
//                                    const arma::mat& matX,
//                                    const arma::uvec& e,
//                                    const arma::vec& breaks,
//                                    const arma::uvec& disc)
// {
//   int P = z.n_rows;
//   int n = e.n_elem;
//   arma::vec::const_iterator bb = breaks.begin();
//   arma::vec::const_iterator be = breaks.end();
//
//   for(int p = 0; p < P; p++)
//   {
//     if(disc(p) == 0)
//     {
//       arma::uvec ind = arma::cumsum(arma::regspace<arma::uvec>(0,n-1));
//       arma::vec zp = matX.col( p );
//       arma::vec zpref = arma::sort(zp( ind ));
//       // BE CAREFUL
//       int j = 0;
//       for(int i = 0; i < n; i++)
//       {
//         if(e(i) == 1)
//         {
//           double z2ecdf = (std::lower_bound(zpref.begin(), zpref.end(), z(p,i) )- zpref.begin())/(n+0.0)  ;
//           z2(p,j) = std::distance(bb, std::upper_bound(bb, be, z2ecdf ) ) + 1 ;
//           j++;
//         }
//       }
//     }
//     else{
//       z2.row(p) = arma::conv_to<arma::urowvec>::from(z.row(p));
//     }
//   }
// }

// Averaging CHF
// arma::vec ForestPrediction::getSurvival(const arma::umat& zt2,
//                                          const arma::vec& y,
//                                          const arma::uvec& e,
//                                          const arma::field<arma::uvec>&& nodeSizeB0,
//                                          const arma::umat&& nodeLabelB0,
//                                          const arma::field<arma::uvec>&& tnd3B0,
//                                          const arma::umat&& ids,
//                                          const arma::field<arma::umat>&& trees)
// {
//   arma::uvec e1 = arma::find(e == 1);
//   arma::vec y2 = y( e1 );
//   arma::uvec idAll = arma::regspace<arma::uvec>(0, y.n_elem-1);
//   arma::uvec idY2 = idAll( e1 );
//   uint NE = y2.n_elem;
//   uint NUM_TREE = trees.size();
//
//   arma::vec Hes = arma::zeros<arma::vec>(y.n_elem);
//
//   for(int b = 0; b != NUM_TREE; b++)
//   {
//     // take the b th tree
//     arma::umat tt = (trees(b));
//     arma::uvec vars =  tt.col(0);
//     arma::uvec values = tt.col(1);
//     arma::uvec lcs = tt.col(2);
//     arma::uvec rcs = tt.col(3);
//     arma::uvec il = tt.col(4);
//     arma::uvec tnd3 = (tnd3B0(b));
//
//     arma::uvec zY;
//
//     // take the ids with uncensored survival times
//     arma::uvec idb = ids.col(b); // IDs have been sorted when sampling is done; need to sort if not sorted before
//     //arma::uvec idb = arma::sort(ids.col(b));
//     arma::uvec idbe = idb( arma::find( e(idb) == 1 )  );
//
//     int j = 0;
//     int nidbe = idbe.n_elem;
//
//     arma::vec wb0 = arma::zeros<arma::vec>(NE);
//     arma::vec wb1 = arma::zeros<arma::vec>(NE);
//
//     for(int i = 0; i < NE; i++)
//     {
//       int count = 0;
//       while(j < nidbe)
//       {
//         if( idbe(j) == idY2(i) )
//         {
//           count++;
//           j++;
//         }else{
//           break;
//         }
//       }
//
//
//       {
//         zY = zt2.col(i);
//
//         int isl = 0;
//         int varsp = 0;
//         int cutsp = 0;
//         int k = 0;
//
//         while(isl == 0)
//         {
//           varsp = vars(k);
//           cutsp = values(k);
//           if(zY(varsp) > cutsp)
//           {
//             k = rcs(k);
//           }else{
//             k = lcs(k);
//           }
//           isl = il(k);
//         }
//
//         arma::uvec nbi = nodeSizeB0(b+NUM_TREE*i);
//
//         if(nodeLabelB0(b,i) == k)
//         {
//           wb0(i) += count ;
//         }
//         wb1(i) += nbi(tnd3(k));
//       }
//
//     }
//     arma::vec wb = arma::zeros<arma::vec>(y.n_elem);
//     wb(e1) = wb0/wb1;
//     wb.replace(arma::datum::nan, 0);
//     Hes += arma::cumsum(wb)/NUM_TREE;
//   }
//
//   return exp(-Hes);
// }
