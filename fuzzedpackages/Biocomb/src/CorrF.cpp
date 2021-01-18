#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

double calc_Gain(NumericVector m1,NumericVector m2,bool symm)
{
  int dd;
  int ik;
  int il;
  int im;
  int len_v1;
  int len_v2;
  double entropyF1=0.0;
  double entropyF2=0.0;
  NumericVector values1;
  NumericVector values2;

  double vrem;

  dd=m1.length();
  //table
  values1=unique(m1);
  len_v1=values1.size();
  NumericVector fq1(len_v1,0.0);
  for(ik=0;ik<dd;ik++)
  {
    for(il=0;il<len_v1;il++)
    {
      if(m1[ik]==values1[il])
      {
        fq1[il]=fq1[il]+1;
        break;
      }
    }
  }


  for(ik=0; ik<len_v1;ik++)
  {
    fq1[ik]=fq1[ik]/dd;
    if(fq1[ik]!=0)
    {
      entropyF1=entropyF1-fq1[ik]*log(fq1[ik]);
    }

  }

  //F2
  dd=m2.length();
  //table
  values2=unique(m2);
  len_v2=values2.size();
  NumericVector fq2(len_v2,0.0);
  for(ik=0;ik<dd;ik++)
  {
    for(il=0;il<len_v2;il++)
    {
      if(m2[ik]==values2[il])
      {
        fq2[il]=fq2[il]+1;
        break;
      }
    }
  }

  for(ik=0; ik<len_v2;ik++)
  {
    fq2[ik]=fq2[ik]/dd;
    if(fq2[ik]!=0)
    {
      entropyF2=entropyF2-fq2[ik]*log(fq2[ik]);
    }

  }

  double entropyF12=0.0;
  NumericMatrix fq(len_v1,len_v2);
  double fq0;

  for(ik=0;ik<dd;ik++)
  {
    for(il=0;il<len_v1;il++)
    {
      if(m1[ik]==values1[il])
      {
        for(im=0;im<len_v2;im++)
        {
          if(m2[ik]==values2[im])
          {
            fq(il,im)=fq(il,im)+1;
            break;
          }
        }
        if(im<len_v2) break;
      }
    }
  }


  for(ik=0;ik<len_v2;ik++)
  {
    vrem=0;
    for(il=0; il<len_v1;il++)
    {
      fq0=fq(il,ik)/sum(fq(_,ik));

      if(fq0!=0)
      {
        vrem=vrem-fq0*log(fq0);
      }
    }
    entropyF12=entropyF12+vrem*fq2[ik];
  }

  double entropy;

  entropy=entropyF1-entropyF12;
    if(symm)
    {
      if((entropyF1+entropyF2)==0)
      {
        entropy=0;
      }
      else
      {
        entropy=2*entropy/(entropyF1+entropyF2);
      }
    }
     return(entropy);
}

NumericVector CalcFeature(DataFrame m3, NumericVector subset, int i)
{
  int ik;
  double out;
  int len;
  len=subset.length();
  NumericVector vrem;
  for(ik=i;ik<len;ik++)
  {
    out=calc_Gain(m3(subset[ik]),m3(subset[i]),true);
    vrem.push_back(out);
  }
   return vrem;
}

double cfs_forward(NumericVector subset, DataFrame m3)
{
 NumericVector vrem;
 double sum;
   sum=0.0;
 double out;
   out=0.0;
 double CFS;
 int i;
 int j;
 int len;
 len=subset.length();
 int ncl;
 ncl=m3.size();
 NumericVector cl;
 cl=m3(ncl-1);

  for(i=0; i<len; i++)
  {
    sum=sum+calc_Gain(m3(subset[i]),cl,true);
    if(len==1)
    {
      out=0.0;
    }
    else
    {
      vrem=CalcFeature(m3, subset, i);

      for(j=1;j<vrem.length();j++)
      {
        out=out+vrem[j];

      }
    }
  }
  out=2*out; //perhaps without len
    CFS=sum/sqrt(out+len); //according to the formula

    return(CFS);
}

// [[Rcpp::export]]
NumericVector forward_path(NumericVector features, DataFrame m3)
{

  NumericVector index;

  double out;
  int atr;
  atr=features.length();
  index=seq(0,(atr-1));

  double maxv;
  int icol=0;
  double maxv_old=-1;

  NumericVector subset(atr,-1.0);
  int pos=0;

  while(atr>0)
  {
    maxv=-1;
    NumericVector tt(subset.begin(),subset.begin()+pos+1);

  for(int i=0;i<index.length();i++)
  {
    tt[pos]=index[i];

    out=cfs_forward(tt, m3);
    if(out>maxv)
    {
      maxv=out;
      icol=i;
    }
  }
  if(maxv>maxv_old)
  {
    subset[pos]=index[icol];
  index.erase(icol);
  atr=index.length();
  pos++;
  maxv_old=maxv;
  }
  else
  {
    break;
  }
  }
  return subset;
}
