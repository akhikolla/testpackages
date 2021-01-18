#include <algorithm>
#include <Rcpp.h>
#include <ctime>
using namespace Rcpp;

IntegerVector samplepoint(int size_point, int size_sample ) {
  IntegerVector indexsp(size_sample);
  IntegerVector indexpt(size_point);
  int i,k;
  
  if (size_point<=size_sample)
  {indexsp=seq_len(size_point)-1 ; }
  else
  {
    RNGScope scope;
    indexpt=seq_len(size_point)-1;
    for (i = 0; i < size_sample; i++)
    {
      k=(as<int>((runif(1,0,size_point-i))));
      indexsp(i)=indexpt(std::min(k,size_point-1-i));
      indexpt(k)=indexpt(size_point-1-i);
    }
  }
  return(indexsp);
}

List selection_pt(NumericVector pt_x, NumericVector pt_y, double coord_x,double coord_y,float bandwith) {
  
  int i,taille;
  double d2;
  double band2=pow(bandwith,2);
  std::vector<double> dist2;
  std::vector<double> pt_ex_x;
  std::vector<double> pt_ex_y;
  taille=0;
  
  
  for (i=0;i<pt_x.length();i++)
  {
    d2=pow(pt_x[i]-coord_x,2)+pow(pt_y[i]-coord_y,2);
    if (d2 <= band2)
    {
      dist2.push_back(d2);
      pt_ex_x.push_back(pt_x[i]);
      pt_ex_y.push_back(pt_y[i]);
      taille++;
    }
  }
  
  
  return List::create(
    _["dist2"] = dist2,
    _["pt_ex_x"] = pt_ex_x,
    _["pt_ex_y"] = pt_ex_y,
    _["taille"] = taille
  );
}



List sample_selection_pt(List selec ,double bandwith,int taille_sample) {
  
  int i,taille,taille_sp;
  double d2,pt_xi,pt_yi,band2;
  std::vector<double> dist2,dist2_sp;
  std::vector<double> pt_ex_x,pt_ex_sp_x;
  std::vector<double> pt_ex_y,pt_ex_sp_y;
  taille=0;
  band2=pow(bandwith,2);
  
  //dans le cercle de obs de rayon_bandwith bandwith
  for (i=0;i<  as<int>(selec["taille"]);i++)
  {
    d2=(as<std::vector<double> >(selec["dist2"]))[i];
    if (d2 <= band2)
    {
      dist2.push_back(d2);
      pt_xi=(as<std::vector<double> >(selec["pt_ex_x"]))[i];
      pt_yi=(as<std::vector<double> >(selec["pt_ex_y"]))[i];
      pt_ex_x.push_back(pt_xi);
      pt_ex_y.push_back(pt_yi);
      taille++;
    }
  }
  
  //on sample les points
  IntegerVector index=samplepoint(taille,taille_sample);
  
  taille_sp=index.length();
  
  for (i=0;i<taille_sp;i++)
  {
    dist2_sp.push_back(dist2[index[i]]);
    pt_ex_sp_x.push_back(pt_ex_x[index[i]]);
    pt_ex_sp_y.push_back(pt_ex_y[index[i]]);
  }
  
  
  return List::create(
    _["dist2"] = dist2_sp,
    _["pt_ex_x"] = pt_ex_sp_x,
    _["pt_ex_y"] = pt_ex_sp_y,
    _["taille"] = taille_sp
  );
}



NumericMatrix dist_sp_pt(List sp, List pt,NumericVector rayons,double bandwith){
  int i,j,taille_sp,taille_pt,taille_rayons,tt,tt2;
  double d2;
  double band2=pow(bandwith,2);
  taille_rayons=rayons.length();
  NumericMatrix nb_pt_rayons(as<int>(sp["taille"]),taille_rayons+1);
  NumericVector compteur(rayons.length());
  taille_sp=as<int>(sp["taille"]);
  taille_pt=as<int>(pt["taille"]);
  NumericVector rayons2(taille_rayons); 
  
  
  for (tt=0;tt<taille_rayons;tt++){rayons2(tt)=pow(rayons(tt),2); }
  
  NumericVector spx(taille_sp);
  NumericVector spy(taille_sp);
  NumericVector ptx(taille_pt);
  NumericVector pty(taille_pt);
  
  for (i=0;i<taille_sp;i++)
  {
    spx(i)=(as<std::vector<double> >(sp["pt_ex_x"]))[i];      
    spy(i)=(as<std::vector<double> >(sp["pt_ex_y"]))[i];      
  }
  for (i=0;i<taille_pt;i++)
  {
    ptx(i)=(as<std::vector<double> >(pt["pt_ex_x"]))[i];      
    pty(i)=(as<std::vector<double> >(pt["pt_ex_y"]))[i];      
  }
  
  for (i=0;i<taille_sp;i++)
  {
    //for(tt=0;tt<taille_rayons;tt++){compteur[tt]=0;}
    nb_pt_rayons(i,0)= pow(1- ((as<std::vector<double> >(sp["dist2"]))[i]/band2),2)  ;      
    
    for (j=0;j<taille_pt;j++)
    {
      d2=pow(ptx(j)-spx(i),2)+pow(pty(j)-spy(i),2);
      
      for (tt=0;tt<taille_rayons;tt++)
      {
        //if (d2<=rayons2(tt)){compteur[tt]++;}
        if (d2<=rayons2(tt)){
          for (tt2=tt;tt2<taille_rayons;tt2++){nb_pt_rayons(i,tt2+1)++;}
          break;    
          //nb_pt_rayons(i,tt)++;
        }
      }
    }
    //for(tt=0;tt<taille_rayons;tt++){nb_pt_rayons(i,tt)=compteur[tt];}  
  }
  return(nb_pt_rayons);
}



NumericVector entropie(NumericMatrix nb, NumericVector qs){
  NumericVector ent(qs.length()*(nb.ncol()-1));
  int i,j,qi;
  double sump=0;
  for (i=0;i<nb.nrow();i++)
  {
    
    sump+=nb(i,0);
    for (j=0;j<(nb.ncol()-1);j++)
    {
      for (qi=0;qi<qs.length();qi++)
      {
        if (qs[qi]==1){ent[j +qi*(nb.ncol()-1)] += nb(i,0)* log2(nb(i,j+1));} 
        else
        {ent[j +qi*(nb.ncol()-1)] += nb(i,0)*pow(nb(i,j+1),qs[qi]-1);}
        
        
      }
    }
  }
  
  
  for (j=0;j<(nb.ncol()-1);j++)
  {
    for (qi=0;qi<qs.length();qi++)
    {
      if (qs[qi]==1){ent[j +qi*(nb.ncol()-1)]= ent[j +qi*(nb.ncol()-1)]/sump;}
      else  {ent[j +qi*(nb.ncol()-1)]= (1/(qs[qi]-1))*log2(ent[j +qi*(nb.ncol()-1)]);}
    }
  }
  
  return(ent);
}

// [[Rcpp::export]]
NumericMatrix gwfa_c(NumericVector pt_x, NumericVector pt_y, NumericVector coord_x,NumericVector coord_y,float bandwith,int taille_sample,NumericVector rayons,NumericVector qs  ) {
  int i,j,k,count,taille_sp,compteur; 
  float band2=pow(bandwith,2);
  NumericVector temp(qs.length()*rayons.length());
  List selec_pt;
  List sample_select_pt;
  NumericMatrix nb_pt_rayons;
  NumericMatrix entropie_sand_box(coord_x.length(),1+qs.length()*rayons.length()); 
  
  
  compteur=0;

  for (i=0;i<coord_x.length();i++)
  {
    compteur=ceil(((float)(i*100))/ coord_x.length());
    Rcpp::Rcout << "\rCalculation in progress : " << compteur<<"%";

    selec_pt=selection_pt(pt_x,pt_y,coord_x(i),coord_y(i),bandwith+max(rayons));//selection des points participants aux calculs pour une obs
    
    count=0;
    
    taille_sp=as<int>(selec_pt["taille"]);
    for (k=0;k<taille_sp;k++)
    {
      if  ((as<std::vector<float> >(selec_pt["dist2"]))[k]<=band2){count++;}      
    }
    entropie_sand_box(i,0)=count;
    
      
    sample_select_pt=sample_selection_pt(selec_pt,bandwith,taille_sample);//selection de echantillon 
    
    //start = std::clock();
    nb_pt_rayons=dist_sp_pt(sample_select_pt,selec_pt,rayons,bandwith);  
    //duration+= std::clock()-start;
    
    temp=entropie(nb_pt_rayons,qs);
    for (j=0;j<(qs.length()*rayons.length());j++){entropie_sand_box(i,j+1)=temp(j); }      
  }
  
  Rcpp::Rcout<<"\rCalculation completed            ";
  return (entropie_sand_box);
}



