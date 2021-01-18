#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


//algorithm code---------------

List SearchF(int lvl, double dvalue, List s_data, int count, NumericVector seq, int countVrem)
{
  List res;
  List ttt;
  int num;
  
  
  
  NumericVector vrem=s_data[seq[lvl+1]-1];
  
  num=vrem.size();
   //cout << " num (lvl)= " << num<<","<<(lvl+1)<<"\n";
   //cout << " vrem = " << vrem[0]<<"\n";
   //cout << " dvalue = " << dvalue<<"\n";
    //cout << " Sp = " << Sp[lvl]<<"\n";
  for (int i=0; i<num; i++)
  {
    
    if (vrem[i] > dvalue)
    {
      if (lvl == (s_data.size() - 2))
      {
        count = count + 1;
        countVrem = countVrem + 1;
        //cout << " count = " << count<<"\n";
        
        
      }
      else
      {
        //S[lvl+1].push_back(i/num);
        res=SearchF(lvl + 1, vrem[i], s_data, count, seq, countVrem);
        count=res[0];
        countVrem=res[1];
      }
    }
  }
  
 
  ttt =	List::create(_["count"]=count, _["countVrem"]=countVrem);
  //cout << " countOut= " << count<<"\n";
  //cout << " countVrem= " << countVrem<<"\n";
  
  return ttt;
}

// [[Rcpp::export]]
List CalcGene(List s_data, NumericMatrix seqAll, int prodValue)
{
  int lvl;
  int count;
  int countVrem;
  int num;
  double HUM;
  int tlen;
  
  
  List out;
  List ttt;
 
  NumericVector vrem;
  tlen=seqAll.ncol();
  NumericVector seq(tlen);
  NumericVector seqMax(tlen);
  double HUMmax;
  HUMmax=-1;
  //cout << " Sp " << Sp[0]<<"\n";
  //cout << " seq " << seq[0]<<","<<seq[1]<<"\n";
  for (int j=0; j< seqAll.nrow();j++)
  {
  
    for(int k=0;k<seqAll.ncol();k++)
    {
      seq[k]=seqAll(j,k);
    }
    //cout << " seq " << seq[0]<<seq[1]<<"\n";
    //cout << "s_data.size()=" << s_data.size()<<"\n";
    
    
  lvl = 0;
  count = 0;
  countVrem=0;
  vrem=s_data[seq[lvl]-1];
  num=vrem.size();
  
  //cout << " num= " << num<<"\n";
   //cout << " j= " << j<<"\n";
  for (int i=0; i< num;i++)
  {
    countVrem=0;
    //cout << " startcount " << count<<"\n";
    out=SearchF(lvl, vrem[i], s_data ,count, seq, countVrem);
    
    count=out[0];
    countVrem=out[1];
    //cout << " endcount " << count<<"\n";
  }

  HUM = ((double)count) / prodValue;
  //cout << "HUM = " << HUM<<"\n";

  if(HUM>HUMmax)
  {
    HUMmax=HUM;
    for(int i=0;i<tlen;i++)
    {
    seqMax[i]=seq[i];
    }
  }
}
//cout << "HUMmax = " << HUMmax<<"\n";
  ttt=  List::create(_["HUM"]=HUMmax, _["seq"]=seqMax);
  return ttt;

}


// [[Rcpp::export]]
List CalcROC(List s_data, NumericVector seq, NumericVector thresholds)
{
  std::vector<double> Sp;
  std::vector<double> Sn;
  std::vector<double> S3;
  double count;
  List ttt;
  double accuracy;
  double optSn;
  double optSp;
  double optS3;
  double optThre;
  double optThre1;
  double optThre2;
  double maxValue;
  double vrem3;
  
  NumericVector vrem;
  
  if(s_data.size()==3)
  {
     maxValue=0;
    for (int i=0; i< thresholds.size();i++)
    {
        count=0;
        vrem=s_data[seq[0]-1];
        for(int j=0;j<vrem.size();j++)
        {
          if(vrem[j]<thresholds[i])
          {
            count++;
          }
          else
          break;
        }
        Sn.push_back(((double)(count))/vrem.size());
        for(int j=i;j< thresholds.size();j++)
        {
          if(j!=i)
          {
            Sn.push_back(Sn.back());
          }
          count=0;
          vrem=s_data[seq[1]-1];
          
          for(int k=0;k<vrem.size();k++)
          {
            if((vrem[k]<thresholds[j])&&(vrem[k]>=thresholds[i]))
            {
              count++;
            }
          }
         
          Sp.push_back(((double)(count))/vrem.size());
          count=0;
          vrem=s_data[seq[2]-1];
         
          for(int k=0;k<vrem.size();k++)
          {
            if(vrem[k]>=thresholds[j])
            {
              count++;
            }
          }
          
          S3.push_back(((double)(count))/vrem.size());
          vrem3=Sp.back()+Sn.back()+S3.back();
          //cout << " vrem3 = " << vrem3<<"\n";
          //cout << " Sn = " << vSn<<"\n";
          //cout << " Sp = " << vSp<<"\n";
          //cout << " S3 = " << vS3<<"\n";
          if(vrem3>maxValue)
          {
            maxValue=vrem3;
            optThre1=thresholds[i];
            optThre2=thresholds[j];
            optSn=Sn.back();
            optSp=Sp.back();
            optS3=S3.back();
          }
          if(j==i)
          {
            for(int k=0;k<i;k++)
            {
            Sn.push_back(Sn.back());
            Sp.push_back(Sp.back());
            S3.push_back(S3.back());
            }
          }
        }
        
        
    }
    accuracy=(optSn + optSp + optS3) / 3;
  }
  else
  {
    maxValue=0;
    for(int i=0;i< thresholds.size();i++)
    {
        
        count=0;
        vrem=s_data[seq[0]-1];
        for(int k=0;k<vrem.size();k++)
        {
            if(vrem[k]<thresholds[i])
            {
              count++;
            }
            else
            break;
        }
        Sn.push_back(((double)(count))/vrem.size());
        count=0;
        vrem=s_data[seq[1]-1];
        for(int k=0;k<vrem.size();k++)
          {
            if(vrem[k]>=thresholds[i])
            {
              count++;
            }
        }
        Sp.push_back(((double)(count))/vrem.size());
        if((Sp.back()+Sn.back())>maxValue)
        {
          maxValue=Sp.back()+Sn.back();
          optThre=thresholds[i];
          optSn=Sn.back();
          optSp=Sp.back();
        }
      }
      accuracy=(optSn + optSp) / 2;
  }
if(s_data.size()==2)
{
ttt=  List::create(_["Sp"]=Sp,_["Sn"]=Sn,_["accuracy"]=accuracy,_["optSp"]=optSp,_["optSn"]=optSn,_["optThre"]=optThre);
}
else
{
ttt=  List::create(_["Sp"]=Sp,_["Sn"]=Sn,_["S3"]=S3,_["accuracy"]=accuracy,_["optSp"]=optSp,_["optSn"]=optSn,_["optS3"]=optS3,_["optThre1"]=optThre1,_["optThre2"]=optThre2);
}

  return ttt; 
  
}
