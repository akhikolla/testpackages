#include "Tree.h"
#include "globals.h"

const arma::uvec& Tree::get_split_vars() const
{
  return split_vars;
}

const arma::uvec& Tree::get_split_values() const
{
  return split_values;
}

const arma::uvec& Tree::get_left_childs() const
{
  return left_childs;
}

const arma::uvec& Tree::get_right_childs() const
{
  return right_childs;
}

const arma::uvec& Tree::get_isLeaf() const
{
  return isLeaf;
}

const arma::uvec& Tree::get_parents() const
{
  return parents;
}

void Tree::setzero(uint i, uint ndcount) {
  uint lid = left_childs(i);
  uint rid = right_childs(i);
  if(lid <= ndcount && isLeaf(lid) == 0) {
    setzero(lid,ndcount);
  }
  if(rid <= ndcount && isLeaf(rid) == 0) {
    setzero(rid,ndcount);
  }
  left_childs(i) = 0;
  right_childs(i) = 0;
  split_values(i) = 0;
  split_vars(i) = 0;
}

void Tree::cut(arma::uvec& nodeTerminal)
{
  arma::uvec isLeaf2 = isLeaf;
  isLeaf2.zeros();
  isLeaf2( nodeTerminal ).ones();
  int ndcount = isLeaf.n_elem-1;
  int ndcount2 = nodeTerminal(nodeTerminal.n_elem-1);
  for(int i = 0; i <= ndcount2; i++) {
    if(isLeaf2(i) == 1 && isLeaf(i) == 0) {
      setzero(i, ndcount);
    }
  }
  arma::uvec nonEmpty = arma::regspace<arma::uvec>(0, ndcount2);
  right_childs = right_childs.elem(nonEmpty);
  left_childs = left_childs.elem(nonEmpty);
  split_values = split_values.elem(nonEmpty);
  split_vars = split_vars(nonEmpty);
  isLeaf = isLeaf2(nonEmpty);
}

double Tree::get_ICONTrain(const arma::uvec& isLeafTemp,
                           const arma::mat& fmat,
                           const arma::umat& Smat)
{
  arma::mat fmatTemp = fmat.cols( arma::find(isLeafTemp == 1) );
  arma::umat SmatTemp = Smat.cols( arma::find(isLeafTemp == 1) );
  int numLeafTemp = arma::sum(isLeafTemp == 1);
  uint K = fmat.n_rows;
  //SmatTemp.print("a");
  arma::vec icon = arma::zeros<arma::vec>(K);
  for(size_t k = 0; k != K; k++) {
    for(int i = 0; i < numLeafTemp; i++) {
      for(int j = 0; j < i; j++) {
        double fi = fmatTemp(k,i);
        double fj = fmatTemp(k,j);
        int Sj = SmatTemp(k,j);
        int Si = SmatTemp(k,i);
        if(fi*Sj > fj*Si) {
          icon(k) = icon(k) + fi*Sj;
        } else {
          icon(k) = icon(k) + fj*Si;
        }
      }
      icon(k) = icon(k) + fmatTemp(k,i)*SmatTemp(k,i)*0.5;
    }
  }
  icon = icon / (arma::sum(fmatTemp, 1) % arma::sum(SmatTemp, 1));
  return arma::sum(icon)/K;
}

//
// void Tree::findOptimalSizekSubtree(arma::mat& fmat, arma::umat& Smat,
//                                    arma::vec& iconAll, arma::field<arma::uvec>& nodeSetList, uint numLeaf)
// {
//   arma::uvec nodeID = arma::regspace<arma::uvec>(0, isLeaf.n_elem-1);
//   arma::uvec isLeafTemp = isLeaf;
//   arma::uvec isLeafTemp2 = isLeafTemp;
//   arma::uvec isLeafTemp3= isLeafTemp;
//   iconAll(numLeaf-1) = get_ICONTrain(isLeaf, fmat, Smat);
//   nodeSetList(numLeaf-1) = nodeID(arma::find(isLeafTemp == 1) );
//   while(numLeaf >= 2)
//   {
//     int cutNd;
//     double iconMax = 0;
//     double iconl = 0;
//     arma::uvec nodeTermTemp = nodeID(arma::find(isLeafTemp == 1) );
//
//     isLeafTemp2 = isLeafTemp;
//     isLeafTemp3 = isLeafTemp;
//
//     if(nodeTermTemp.n_elem > 1)
//     {
//       for(int l = 0; l < nodeTermTemp.n_elem; l++)
//       {
//         isLeafTemp2 = isLeafTemp3;
//         cutNd = nodeTermTemp(l);
//
//         uint pid = parents(cutNd);
//         if( arma::sum( nodeTermTemp == right_childs(pid) )>0 && arma::sum( nodeTermTemp == left_childs(pid))>0 )
//         {
//           isLeafTemp2( pid ) = 1;
//           isLeafTemp2( right_childs(pid) ) = 0;
//           isLeafTemp2( left_childs(pid) ) = 0;
//
//           iconl = get_ICONTrain(isLeafTemp2, fmat, Smat);
//           if(iconl >= iconMax) // = is used for the case with only one terminal nodes
//           {
//             iconMax = iconl;
//             nodeSetList(numLeaf-2) = nodeID(arma::find(isLeafTemp2 == 1) );
//             isLeafTemp = isLeafTemp2;
//           }
//         }
//
//       }
//     }else{
//     }
//
//     iconAll(numLeaf-2) = iconMax;
//     numLeaf--;
//   }
//
// }



void Tree::findOptimalSizekSubtree(arma::mat& fmat, arma::umat& Smat,
                                   arma::vec& iconAll, arma::field<arma::uvec>& nodeSetList,
				   uint numLeaf)
{
  arma::uvec nodeID = arma::regspace<arma::uvec>(0, isLeaf.n_elem-1);
  arma::uvec isLeafTemp = arma::zeros<arma::uvec>(isLeaf.n_elem);
  isLeafTemp(0) = 1;
  iconAll(0) = get_ICONTrain(isLeafTemp, fmat, Smat);
  nodeSetList(0) = nodeID(arma::find(isLeafTemp == 1) );
  arma::uvec isLeafTemp2 = isLeafTemp;
  arma::uvec isLeafTemp3 = isLeafTemp;
  size_t i = 1;
  while(i < numLeaf - 1) {
    double iconMax = 0;
    double iconl = 0;
    isLeafTemp2 = isLeafTemp;
    isLeafTemp3 = isLeafTemp;
    arma::uvec nodeTermTemp = nodeID(arma::find(isLeafTemp == 1) );
    for(size_t l = 0; l < nodeTermTemp.n_elem; l++ ) {
      int spNd = nodeTermTemp(l);
      if(isLeaf(spNd) == 0) {
        isLeafTemp2 = isLeafTemp3;
        int lid = left_childs(spNd);
        int rid = right_childs(spNd);
        isLeafTemp2(spNd) = 0;
        isLeafTemp2(lid) = 1;
        isLeafTemp2(rid) = 1;
        iconl = get_ICONTrain(isLeafTemp2, fmat, Smat);
        if(iconl > iconMax) {
          iconMax = iconl;
          nodeSetList(i) = nodeID(arma::find(isLeafTemp2 == 1) );
          isLeafTemp = isLeafTemp2;
        }
      }
    }
    iconAll(i) = iconMax;
    i++;
  }
  iconAll(i) = get_ICONTrain(isLeaf, fmat, Smat);
  nodeSetList(i) = nodeID(arma::find(isLeaf == 1) );
}

void Tree::findBeta(arma::vec& iconAll, arma::vec& beta, arma::uvec& sizeTree)
{
  arma::vec alpha(iconAll.n_elem);
  int L = iconAll.n_elem;
  size_t q = 1;
  alpha(0) = 0;
  sizeTree(0) = L;
  while( L > 1 ) {
    arma::vec iconSmallerTree = iconAll.head( L-1 );
    //arma::vec alphaTT = (iconAll(L-1) - iconSmallerTree)/(  arma::regspace<arma::uvec>(L-1,-1,1)  ); // L - 1 to L - (L-1)
    // arma::vec alphaTT = (iconAll(L-1) - iconSmallerTree)/(pow(L,1) - arma::pow(arma::regspace<arma::vec>(1,L-1),1)  ); // L - 1 to L - (L-1)
    arma::vec LL = arma::regspace(L - 1, 1);
    arma::vec alphaTT = (iconAll(L-1) - iconSmallerTree) / LL;
    // Rcpp::Rcout << "alphaTT" << std::endl;
    // Rcpp::Rcout << alphaTT << std::endl; 
    alpha(q) = alphaTT.min();
    sizeTree(q) = alphaTT.index_min() + 1;
    L = sizeTree(q);
    q++;
  }
  if(q < iconAll.n_elem) {
    sizeTree.shed_rows(q, iconAll.n_elem-1);
    alpha.shed_rows(q, iconAll.n_elem-1);
    beta.shed_rows(q, iconAll.n_elem-1);
  }
  for(size_t i = 0; i < alpha.n_elem; i++) {
    if(i < alpha.n_elem - 1) {
      beta(i) = sqrt( alpha(i)*alpha(i+1)  );
    } else {
      beta(i) = alpha(i);
    }
  }
}
