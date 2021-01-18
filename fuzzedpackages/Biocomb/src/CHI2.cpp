#include <Rcpp.h>
#include <queue>

using namespace Rcpp;
using namespace std;


IntegerVector top_i_pq(NumericVector v, int n)
{
  typedef pair<double, int> Elt;
  priority_queue< Elt, vector<Elt>, greater<Elt> > pq;
  vector<int> result;

  for (int i = 0; i != v.size(); ++i) {
    if (((int)pq.size()) < n)
      pq.push(Elt(v[i], i));
    else {
      Elt elt = Elt(v[i], i);
      if (pq.top() < elt) {
        pq.pop();
        pq.push(elt);
      }
    }
  }

  result.reserve(pq.size());
  while (!pq.empty()) {
    result.push_back(pq.top().second + 1);
    pq.pop();
  }

  return wrap(result);
}

NumericVector in_range(NumericVector x, double low, double high) {
  return x[(x >= low) & (x < high)];
}


//algorithm code---------------
// [[Rcpp::export]]
List fun1_chi(DataFrame data,CharacterVector classI)
{
  NumericVector ind;
  NumericVector vrem;
  NumericVector var1;
  CharacterVector class1;
  CharacterVector cl;
  int index;
  int indexS;
  int nFeat;
  int nCase;

  NumericVector index2;
  NumericVector index22;

  nFeat=data.length();
  nCase=data.nrows();


  List int_list(nFeat);
  List alist;
  List ttt;

  index2=seq(0,(nCase-1));
  alist=List::create();
  for (int i=0; i<nFeat; i++)
  {
    alist.erase(alist.begin(),alist.end());
    vrem=data(i);
    ind=top_i_pq(vrem, nCase);
    for(int jj=0;jj<nCase;jj++)
    {
    ind[jj]=ind[jj]-1;
    }

    vrem=vrem[ind];
    cl=classI[ind];

    indexS=0;
    index=0;

    while (index<nCase) {
      if(vrem[index]==vrem[indexS])
      {
        index=index+1;
      }
      else
      {
        index22=in_range(index2,indexS,index);

        var1=ind[index22];
        class1=cl[index22];
        ttt =	List::create(_["var1"]=var1, _["class1"]=class1);

        alist.push_back(ttt);
        indexS=index;
      }
    }
    index22=in_range(index2,indexS,index);
    var1=ind[index22];
    class1=cl[index22];
    ttt =	List::create(_["var1"]=var1, _["class1"]=class1);
    alist.push_back(ttt);

    int_list[i]=alist;
  }

  return int_list;
}

double check_stat(NumericMatrix mat_int)
{
  int dm_col=mat_int.ncol();
  int dm_row=mat_int.nrow();
  double result;

  NumericVector row_sum;
  NumericVector col_sum;
  NumericMatrix res(dm_row,dm_col);

  for(int i=0; i<dm_col; i++)
  {
    col_sum.push_back(sum(mat_int(_,i)));
  }
  for(int i=0; i<dm_row; i++)
  {
    row_sum.push_back(sum(mat_int(i,_)));
  }
  for(int i=0;i<dm_row;i++)
  {
    for(int j=0;j<dm_col;j++)
    {
      res(i,j)=row_sum[i]*col_sum[j];
      res(i,j)=res(i,j)/sum(mat_int);
      if(res(i,j)==0)
      {
        res(i,j)=0.1;
      }
      res(i,j)=pow(mat_int(i,j)-res(i,j),2)/res(i,j);

    }
  }

  result=sum(res);
  return(result);
}

// [[Rcpp::export]]
List fun2_chi(List int_l,NumericMatrix mat_int)
{
  int dm;
  int till;
  NumericVector chi_vrem;
  CharacterVector outChar;
  CharacterVector vv1;
  CharacterVector names;

  List int_list=clone(int_l);
  List vrem;
  List vv;
  int ik;
  int il;

  int indexS;

  dm=mat_int.ncol();
  names=colnames(mat_int);

  till=int_list.size();

  List chi_stat(till);

  for (int i=0; i<till; i++)
  {
    vrem=int_list[i];
    chi_vrem.erase(chi_vrem.begin(),chi_vrem.end());

    for (int j=0; j<(vrem.size()-1); j++)
    {
      for(int ij=0; ij<2; ij++)
      {
        if((j>0)&&(ij==0))
        {
          mat_int(ij,_)=mat_int(ij+1,_);
        }
        else
        {
          vv=vrem[j+ij];
          vv1=vv[1];

          outChar=vv1.sort();
          indexS=0;
          mat_int(ij,_)=rep(0,dm);
          for(ik=0; ik<outChar.size();ik++)
          {
            if(outChar[ik]!=outChar[indexS])
            {
              for(il=0; il<dm;il++)
              {
                if(outChar[indexS]==names[il])
                {
                  mat_int(ij,il)=ik-indexS;
                }
              }
              indexS=ik;
            }
          }
          for(il=0; il<dm;il++)
          {
            if(outChar[indexS]==names[il])
            {
              mat_int(ij,il)=ik-indexS;//ik-indexS+1
            }
          }
        }
      }


      chi_vrem.push_back(check_stat(mat_int));

    }
    chi_stat[i]=chi_vrem;
  }
  return chi_stat;
}

//
// NumericVector vecunion(NumericVector v1,NumericVector v2) {
//   // Rcpp supports STL-style iterators
//
//   NumericVector dest1;
//
//   std::set_union(v1.begin(), v1.end(),
//                  v2.begin(), v2.end(),
//                  std::back_inserter(dest1));
//
//
//   // we want the value so dereference
//   return dest1;
// }

Rcpp::List removeElementList(Rcpp::List x, int j)
{
  IntegerVector idx = seq_len(x.length());
  return(x[idx != j]);
}

int vecminInd(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference
  return it - x.begin();
}

// [[Rcpp::export]]
List fun3_chi(List chi_s,List int_l,DataFrame datain, double chi_value, NumericMatrix mat_int)
{
  int d1;
  int dm;
  int ik;
  int il;

  int irow;
  int icol;
  int check;
  int calc;

  List vv1;
  List vv2;
  List vv;
  NumericVector a1;
  NumericVector a2;
  CharacterVector b1;
  CharacterVector b2;

  NumericVector vrem_chi;
  List vrem;
  CharacterVector outChar;
  CharacterVector names;
  CharacterVector vv_char;
  CharacterVector vv_ch;

  NumericVector vec_index;
  NumericVector vdata;

   List chi_stat = clone(chi_s);
   List int_list=clone(int_l);
  DataFrame data=clone(datain);

  d1=data.size();
  dm=mat_int.ncol();
  names=colnames(mat_int);

  int flag=0;

  for(irow=0; irow<d1; irow++)
  {
    calc=0;
    while(true)
    {
      //select the min element of chi.stat

      vrem_chi=chi_stat[irow];
      check=vrem_chi.size();

      if(check==0)
      {
        break;
      }
      //icol=which_min(vrem_chi); //which interval
      double minv=1000;
      int index3=-1;
      for(int id=0;id<vrem_chi.size();id++)
      {
        if(vrem_chi[id]<(minv-0.0001))
        {

          minv=vrem_chi[id];
          index3=id;
        }
      }
      icol=index3;

        if(vrem_chi[icol]>chi_value)
        {
          break;
        }


        vrem=int_list[irow];

        vv1=vrem[icol];
        vv2=vrem[icol+1];
        a1=vv1[0];
        a2=vv2[0];
        b1=vv1[1];
        b2=vv2[1];


        for(int id=0; id<a2.size();id++)
        {
          a1.push_back(a2[id]);
        }

        for(int id=0; id<b2.size();id++)
        {
          b1.push_back(b2[id]);
        }

        vv1[0]=a1;
        vv1[1]=b1;
        vrem[icol]=vv1; //vrem.insert(icol,vv1);


          if(icol!=(vrem.size()-2))
          {
            for( int i=(icol+1); i<(vrem.size()-1); i++)
            {
              vrem[i]=vrem[i+1];
              if(i==(vrem.size()-2)) break;
                vrem_chi[i]=vrem_chi[i+1];
            }
          }


          vrem.erase(vrem.size()-1);
          vrem_chi.erase(vrem.size()-1);

        //new chi values

          for(int j=icol-1; j<=icol; j++)
          {
            if((j>-1)&&(j<(vrem.size()-1)))
            {
              for(int ij=0; ij<2; ij++)
              {
                vv=vrem[j+ij];
                vv_ch=vv[1];
                //without sort
                mat_int(ij,_)=rep(0,dm);
                for(ik=0;ik<vv_ch.size();ik++)
                {
                  for(il=0;il<names.size();il++)
                  {
                    if(vv_ch[ik]==names[il])
                    {
                      mat_int(ij,il)=mat_int(ij,il)+1;
                      break;
                    }
                  }
                }
              }
              vrem_chi[j]=check_stat(mat_int);
            }
          }

          vv1=vrem[icol];
          vec_index=vv1[0];
          vdata=data(irow);

          vdata[vec_index]=vdata[vec_index[0]];
          data(irow)=vdata;

          int_list[irow]=vrem;
          chi_stat[irow]=vrem_chi;

          calc++;
    }
    if(flag==1)
      break;
  }
  List ttt=  List::create(_["data"]=data, _["chi_stat"]=chi_stat, _["int_list"]=int_list);
  return ttt;
}

// NumericVector myRow(DataFrame& x,int nrow) {
//
//
// int nCols=x.size();
//   NumericVector y(nCols);
//   for (int j=0; j<nCols;j++) {
//     NumericVector column = x[j] ;
//     y[j] = column[nrow] ;
//   }
//   return y;
// }

double pdistC(NumericVector x, NumericVector y) {

  int n = y.size();
  double out=0;

  for(int i = 0; i < n; i++) {
    out =out + pow((y[i] - x[i]),2.0);
  }
  return out;
}

NumericVector lapply1(NumericMatrix input, int index) {

  int n = input.nrow();
  NumericVector out(n);

  for(int i = 0; i < n; i++) {
    out[i] = pdistC(input(i,_),input(index,_));
  }

  return out;
}

NumericMatrix DFtoNM(DataFrame x,int nRows) {
  //int nRows=x.nrows();
  NumericMatrix y(nRows,x.size());
  for (int i=0; i<x.size();i++) {
    y(_,i)=NumericVector(x(i));
  }
  return y;
}

CharacterVector new_class(CharacterVector classI,NumericVector out)
{
  int i;
  out = out.sort();

  for(i=(out.length()-1);i>=0;i--)
  {
    classI.erase(out[i]);
  }

  return classI;
}

NumericMatrix row_erase (NumericMatrix& x, NumericVector& rowID) {

  rowID = rowID.sort();

  NumericMatrix x2(Dimension(x.nrow()- rowID.size(), x.ncol()));
  int iter = 0;
  int del = 1; // to count deleted elements
  int index=0;

  for (int i = 0; i < x.nrow(); i++) {

    if (i != rowID[del - 1]) {
      x2.row(iter) = x.row(i);
      iter++;
    } else {
      del++;
    }

    if((del-1)==(rowID.size()))
    {
      index=i+1;
      break;
    }
  }

  if(index!=x.nrow())
  {
  for(int i=index;i<x.nrow();i++)
  {
    x2.row(iter+(i-index)) = x.row(i);
  }
  }
  return x2;
}

// [[Rcpp::export]]
double check_incons(DataFrame data,DataFrame vrem_nom,CharacterVector cl)
{
  int i;
  double incons;
  int ik;
  int indexS;

  NumericVector vrem;
  DataFrame data1=clone(data);
  CharacterVector classI=clone(cl);
  DataFrame vrem_nominal=clone(vrem_nom);

  CharacterVector classVr;
  CharacterVector outChar;
  NumericVector counts;
  NumericMatrix data_check;


  if(vrem_nominal.size()>0)
  {
    for(i=0; i<vrem_nominal.size();i++)
    {
      data1.push_back(vrem_nominal(i));
    }
  }

  int num_row=data.nrows();

  data_check=DFtoNM(data1,num_row);

  int input;
  input=num_row;

  NumericVector out;
  incons=0;
  //int calc=0;
  while(input>1)
  {
    vrem=lapply1(data_check, 0);
    out.erase(out.begin(),out.end());
    for(i=0; i<vrem.length(); i++)
    {
      if(vrem[i]==0)
      {
        out.push_back(i);
      }
    }

    int maxc=0;
    int vc;
    if(out.length()>1)
    {
      classVr=classI[out];

      outChar=classVr.sort();
      indexS=0;
      for(ik=0; ik<outChar.size();ik++)
      {
        if(outChar[ik]!=outChar[indexS])
        {
          vc=ik-indexS;
          if(vc>maxc) maxc=vc;
          indexS=ik;
        }
      }

      vc=ik-indexS;
      if(vc>maxc) maxc=vc;

      incons=incons+(out.length()-maxc);
    }

    input=input-out.length();

    if(input<2) break;
    classI=new_class(classI,out);
    data_check=row_erase(data_check,out);
  }

  incons=incons/num_row;

  return incons;
}

//[[Rcpp::export]]
DataFrame fun4_chi(List chi_s,List int_l,DataFrame datain,DataFrame vrem_nominal,NumericVector chi_attr, NumericVector sig_attr,CharacterVector cl,NumericMatrix mat_int, double threshold,int df,double step,int delta, int shag)
{
  int irow;
  int check;
  int icol;
  int ik;
  int il;
  int dm;
  double incons;

  List vv1;
  List vv2;
  List vv;
  NumericVector a1;
  NumericVector a2;
  CharacterVector b1;
  CharacterVector b2;

  NumericVector vec_index;
  NumericVector vdata;
  NumericVector vrem_chi;
  NumericVector data1;
  NumericVector data11;

  List vrem;
  List int_list1;
  List int_list11;
  NumericVector chi_stat1;
  NumericVector chi_stat11;

  CharacterVector outChar;
  CharacterVector names;
  CharacterVector vv_char;
  CharacterVector vv_ch;

  //int index3;

  DataFrame data=clone(datain);
  List int_list=clone(int_l);
  List chi_stat=clone(chi_s);

  int d1=data.size();
  dm=mat_int.ncol();
  names=colnames(mat_int);

  NumericVector flag(d1,1.0);
  //further
  int calc=0;
  while(true)
  {
    if(max(flag)==0) break;
    for(irow=0; irow<d1; irow++)
    {
      if(flag[irow]==1.0)
      {
        vrem=int_list[irow];
        vrem_chi=chi_stat[irow];
        //for the recovery
        data11=data(irow);
        data1=clone(data11);
        int_list11=int_list[irow];
        int_list1=clone(int_list11);
        chi_stat11=chi_stat[irow];
        chi_stat1=clone(chi_stat11);

        while(true)
        {
          //select the min element of chi.stat
          check=vrem_chi.size();
          if(check==0)
          {
            flag[irow]=0;
            break;
          }
            //icol=which_min(vrem_chi); //which interval
            //icol=vecminInd(vrem_chi);
            double minv=1000;
          int index3=-1;
          for(int id=0;id<vrem_chi.size();id++)
          {
            if(vrem_chi[id]<(minv-0.0001))
            {
              minv=vrem_chi[id];
              index3=id;
            }
          }
          icol=index3;
            if(vrem_chi[icol]>chi_attr[irow])
            {
              break;
            }
            vrem=int_list[irow];

            vv1=vrem[icol];
            vv2=vrem[icol+1];
            a1=vv1[0];
            a2=vv2[0];
            b1=vv1[1];
            b2=vv2[1];


            for(int id=0; id<a2.size();id++)
            {
              a1.push_back(a2[id]);
            }
            for(int id=0; id<b2.size();id++)
            {
              b1.push_back(b2[id]);
            }

            vv1[0]=a1;
            vv1[1]=b1;
            vrem[icol]=vv1; //vrem.insert(icol,vv1);


              if(icol!=(vrem.size()-2))
              {
                for( int i=(icol+1); i<(vrem.size()-1); i++)
                {
                  vrem[i]=vrem[i+1];
                  if(i==(vrem.size()-2)) break;
                  vrem_chi[i]=vrem_chi[i+1];
                }
              }


              vrem.erase((vrem.size()-1));
              vrem_chi.erase((vrem.size()-1));

              for(int j=icol-1; j<=icol; j++)
              {
                if((j>-1)&&(j<(vrem.size()-1)))
                {
                  for(int ij=0; ij<2; ij++)
                  {
                    vv=vrem[j+ij];
                    vv_ch=vv[1];

                    //without sort
                    mat_int(ij,_)=rep(0,dm);
                    for(ik=0;ik<vv_ch.size();ik++)
                    {
                      for(il=0;il<names.size();il++)
                      {
                        if(vv_ch[ik]==names[il])
                        {
                          mat_int(ij,il)=mat_int(ij,il)+1;
                          break;
                        }
                      }
                    }
                  }
                  vrem_chi[j]=check_stat(mat_int);
                }
              }

              vv1=vrem[icol];
              vec_index=vv1[0];
              vdata=data(irow);
              vdata[vec_index]=vdata[vec_index[0]];
              data(irow)=vdata;

              int_list[irow]=vrem;
              chi_stat[irow]=vrem_chi;

                calc++;
        }
        //check.incons
        incons=check_incons(data,vrem_nominal,cl);

          if(incons<=threshold)
          {
            sig_attr[irow]=sig_attr[irow]-step;


            chi_attr[irow]=R::qchisq(1-sig_attr[irow], df,0,1);
          }
          else
          {
            int_list[irow]=int_list1;
            chi_stat[irow]=chi_stat1;
            data(irow)=data1;
            flag[irow]=0;
          }
      }
    }
    if(shag==delta)
    {
      step=step*0.1;
      delta=delta+9;
    }
    shag=shag+1;
}
  return data;
}

