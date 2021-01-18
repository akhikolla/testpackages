#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

static double theta(arma::vec p, int size, double type){
  double res;
  if (type==00)
  {
    res = pow(p(0),size);
  }
  else if (type==10)
  {
    res = pow(p(0)+p(1),size)-pow(p(0),size);
  }
  else if (type==01)
  {
    res = pow(p(0)+p(2),size)-pow(p(0),size);
  }
  else
  {
    res = 1-pow(p(0)+p(1),size)-pow(p(0)+p(2),size)+pow(p(0),size);
  }
  return res;
}

//Function to summarize all possible group sizes
arma::mat poolsize(int MAMPS,int S, std::string algorithm){
  arma::mat res(0,S,fill::ones);
  if (algorithm=="optimal"){
    if (S==2){
      res.resize(MAMPS-1,S);
      for (int i=2; i<=MAMPS;i++){
        res(i-2,0)=i; res(i-2,1)=1;
      }
    }

    if (S==3){
      for (int h=2; h<=floor(MAMPS/pow(2,S-2));h++){
        for(int h1=2; h1<=floor(MAMPS/pow(2,S-2));h1++){
          int length=res.n_rows;
          if (h*h1<=MAMPS){
            res.resize(length+1,S);
            res(length,0)=h*h1; res(length,1)=h1; res(length,2)=1;
          }
        }
      }
    }

    if (S==4){
      for (int h=2; h<=floor(MAMPS/pow(2,S-2));h++){
        for(int h1=2; h1<=floor(MAMPS/pow(2,S-2));h1++){
          for(int h2=2; h2<=floor(MAMPS/pow(2,S-2));h2++){
            int length=res.n_rows;
            if (h*h1*h2<=MAMPS){
              res.resize(length+1,S);
              res(length,0)=h*h1*h2; res(length,1)=h1*h2; res(length,2)=h2; res(length,3)=1;
            }
          }
        }
      }
    }

    if (S==5){
      for (int h=2; h<=floor(MAMPS/pow(2,S-2));h++){
        for(int h1=2; h1<=floor(MAMPS/pow(2,S-2));h1++){
          for(int h2=2; h2<=floor(MAMPS/pow(2,S-2));h2++){
            for(int h3=2; h3<=floor(MAMPS/pow(2,S-2));h3++){
              int length=res.n_rows;
              if (h*h1*h2*h3<=MAMPS){
                res.resize(length+1,S);
                res(length,0)=h*h1*h2*h3; res(length,1)=h1*h2*h3; res(length,2)=h2*h3; res(length,3)=h3; res(length,4)=1;
              }
            }
          }
        }
      }
    }

    if (S==6){
      for (int h=2; h<=floor(MAMPS/pow(2,S-2));h++){
        for(int h1=2; h1<=floor(MAMPS/pow(2,S-2));h1++){
          for(int h2=2; h2<=floor(MAMPS/pow(2,S-2));h2++){
            for(int h3=2; h3<=floor(MAMPS/pow(2,S-2));h3++){
              for(int h4=2; h4<=floor(MAMPS/pow(2,S-2));h4++){
                int length=res.n_rows;
                if (h*h1*h2*h3*h4<=MAMPS){
                  res.resize(length+1,S);
                  res(length,0)=h*h1*h2*h3*h4; res(length,1)=h1*h2*h3*h4; res(length,2)=h2*h3*h4; res(length,3)=h3*h4; res(length,4)=h4;res(length,5)=1;
                }
              }
            }
          }
        }
      }
    }} else if (algorithm=="halving"){
      double mid=pow(2,S-1);
      while(mid<MAMPS){
        int length=res.n_rows;
        res.resize(length+1,S);
        res(length,0)=mid;
        mid+=pow(2,S-2);
      }

      if (S>2){
        for (int i=1;i<=S-2;i++){
          res(span::all,i)=res(span::all,i-1)/2;
        }
      }
      res(span::all,S-1).fill(1);
    }
    return res;
}

//Function to calculate expected number of tests per individual for 2 infections
// [[Rcpp::export]]
double EFF_HIER(arma::vec p, int S, arma::mat SE, arma::mat SP, arma::vec ns){
  //double se1=SE(0), se2=SE(1), sp1=SP(0), sp2=SP(1);
  int n1=ns(0);
  arma::mat mid(4,4);
  arma::mat M(4,4,fill::zeros); M(0,0)=theta(p,n1,00); M(1,1)=theta(p,n1,10); M(2,2)=theta(p,n1,01); M(3,3)=theta(p,n1,11);
  arma::mat P_operator(4,4,fill::zeros);
  P_operator(0,0)=1-SP(0,0)*SP(0,1);
  P_operator(1,1)=1-(1-SE(0,0))*SP(0,1);
  P_operator(2,2)=1-SP(0,0)*(1-SE(0,1));
  P_operator(3,3)=1-(1-SE(0,0))*(1-SE(0,1));

  double ET=1+n1/ns(1)*accu(M*P_operator);
  arma::mat pi_transition(4,(S-2)*4,fill::zeros);

  if (S>2){

    for (int i=0;i<=S-3;i++){
      pi_transition(0,(i*4))=1;
      pi_transition(1,(i*4))=theta(p,ns(i+1),00)*theta(p,ns(i)-ns(i+1),10)/theta(p,ns(i),10);
      pi_transition(1,(i*4+1))=1-pi_transition(1,(i*4));
      pi_transition(2,(i*4))=theta(p,ns(i+1),00)*theta(p,ns(i)-ns(i+1),01)/theta(p,ns(i),01);
      pi_transition(2,(i*4+2))=1-pi_transition(2,(i*4));
      pi_transition(3,(i*4))=theta(p,ns(i+1),00)*theta(p,ns(i)-ns(i+1),11)/theta(p,ns(i),11);
      pi_transition(3,(i*4+1))=theta(p,ns(i+1),10)*(theta(p,ns(i)-ns(i+1),01)+theta(p,ns(i)-ns(i+1),11))/theta(p,ns(i),11);
      pi_transition(3,(i*4+2))=theta(p,ns(i+1),01)*(theta(p,ns(i)-ns(i+1),10)+theta(p,ns(i)-ns(i+1),11))/theta(p,ns(i),11);
      pi_transition(3,(i*4+3))=1-pi_transition(3,(i*4))-pi_transition(3,(i*4+1))-pi_transition(3,(i*4+2));
    }

    for (int s=2; s<=S-1;s++){
      mid= M *P_operator;
      for (int j=0;j<=s-2;j++){
        arma::mat P_operator_temp(4,4,fill::zeros);
        P_operator_temp(0,0)=1-SP(j+1,0)*SP(j+1,1);
        P_operator_temp(1,1)=1-(1-SE(j+1,0))*SP(j+1,1);
        P_operator_temp(2,2)=1-SP(j+1,0)*(1-SE(j+1,1));
        P_operator_temp(3,3)=1-(1-SE(j+1,0))*(1-SE(j+1,1));
        mid = mid * pi_transition(span(0,3), span(j*4,j*4+3)) * P_operator_temp;
      }
      ET=ET+ns(0)/ns(s)*accu(mid);
    }
  }


  double res=ET/n1;

  return res;
}

// [[Rcpp::export]]
arma::mat prob(arma::vec p, int S, arma::mat SE, arma::mat SP, arma::vec ns){
  //double se1=SE(0), se2=SE(1), sp1=SP(0), sp2=SP(1);
  int n1=ns(0);
  arma::mat mid(4,4);
  arma::mat M(4,4,fill::zeros); M(0,0)=theta(p,n1,00); M(1,1)=theta(p,n1,10); M(2,2)=theta(p,n1,01); M(3,3)=theta(p,n1,11);
  arma::mat P_operator(4,4,fill::zeros);
  P_operator(0,0)=1-SP(0,0)*SP(0,1);
  P_operator(1,1)=1-(1-SE(0,0))*SP(0,1);
  P_operator(2,2)=1-SP(0,0)*(1-SE(0,1));
  P_operator(3,3)=1-(1-SE(0,0))*(1-SE(0,1));

  double ET=1+n1/ns(1)*accu(M*P_operator);
  arma::mat pi_transition(4,(S-2)*4,fill::zeros);

  if (S>2){

    for (int i=0;i<=S-3;i++){
      pi_transition(0,(i*4))=1;
      pi_transition(1,(i*4))=theta(p,ns(i+1),00)*theta(p,ns(i)-ns(i+1),10)/theta(p,ns(i),10);
      pi_transition(1,(i*4+1))=1-pi_transition(1,(i*4));
      pi_transition(2,(i*4))=theta(p,ns(i+1),00)*theta(p,ns(i)-ns(i+1),01)/theta(p,ns(i),01);
      pi_transition(2,(i*4+2))=1-pi_transition(2,(i*4));
      pi_transition(3,(i*4))=theta(p,ns(i+1),00)*theta(p,ns(i)-ns(i+1),11)/theta(p,ns(i),11);
      pi_transition(3,(i*4+1))=theta(p,ns(i+1),10)*(theta(p,ns(i)-ns(i+1),01)+theta(p,ns(i)-ns(i+1),11))/theta(p,ns(i),11);
      pi_transition(3,(i*4+2))=theta(p,ns(i+1),01)*(theta(p,ns(i)-ns(i+1),10)+theta(p,ns(i)-ns(i+1),11))/theta(p,ns(i),11);
      pi_transition(3,(i*4+3))=1-pi_transition(3,(i*4))-pi_transition(3,(i*4+1))-pi_transition(3,(i*4+2));
    }

    for (int s=2; s<=S-1;s++){
      mid= M *P_operator;
      for (int j=0;j<=s-2;j++){
        arma::mat P_operator_temp(4,4,fill::zeros);
        P_operator_temp(0,0)=1-SP(j+1,0)*SP(j+1,1);
        P_operator_temp(1,1)=1-(1-SE(j+1,0))*SP(j+1,1);
        P_operator_temp(2,2)=1-SP(j+1,0)*(1-SE(j+1,1));
        P_operator_temp(3,3)=1-(1-SE(j+1,0))*(1-SE(j+1,1));
        mid = mid * pi_transition(span(0,3), span(j*4,j*4+3)) * P_operator_temp;
      }
      ET=ET+ns(0)/ns(s)*accu(mid);
    }
  }

  return mid;
}

//Function to calculate classification accuracy measures for 2 infections
// [[Rcpp::export]]
List ACCU_HIER(arma::vec p, int S, arma::mat SE, arma::mat SP, arma::vec ns){
  double se1=SE(0,0),se2=SE(0,1),sp1=SP(0,0),sp2=SP(0,1);
  int n1=ns(0);

  arma::mat M(4,4,fill::zeros); M(0,0)=theta(p,n1-1,00); M(1,1)=theta(p,n1-1,10); M(2,2)=theta(p,n1-1,01); M(3,3)=theta(p,n1-1,11);
  arma::mat P_operator(4,4,fill::zeros);P_operator(0,0)= 1-sp1*sp2; P_operator(1,1)= 1-(1-se1)*sp2;P_operator(2,2)= 1-sp1*(1-se2);P_operator(3,3)= 1-(1-se1)*(1-se2);
  arma::mat P_pos_neg(4,4,fill::zeros); P_pos_neg(0,0)=1-(1-se1)*sp2; P_pos_neg(1,1)=1-(1-se1)*sp2; P_pos_neg(2,2)=1-(1-se1)*(1-se2); P_pos_neg(3,3)=1-(1-se1)*(1-se2);
  arma::mat P_neg_pos(4,4,fill::zeros); P_neg_pos(0,0)=1-sp1*(1-se2); P_neg_pos(1,1)=1-(1-se1)*(1-se2); P_neg_pos(2,2)=1-sp1*(1-se2); P_neg_pos(3,3)=1-(1-se1)*(1-se2);
  arma::mat P_pos_pos(4,4,fill::zeros); P_pos_pos.diag().fill(1-(1-se1)*(1-se2));

  arma::mat PSe1_mid1(4,4,fill::zeros);  PSe1_mid1=M*P_pos_neg;
  arma::mat PSe1_mid2(4,4,fill::zeros);  PSe1_mid2=M*P_pos_pos;
  arma::mat PSe2_mid1(4,4,fill::zeros);  PSe2_mid1=M*P_neg_pos;
  arma::mat PSe2_mid2(4,4,fill::zeros);  PSe2_mid2=M*P_pos_pos;
  arma::mat PSp1_mid1(4,4,fill::zeros);  PSp1_mid1=M*P_operator;
  arma::mat PSp1_mid2(4,4,fill::zeros);  PSp1_mid2=M*P_neg_pos;
  arma::mat PSp2_mid1(4,4,fill::zeros);  PSp2_mid1=M*P_operator;
  arma::mat PSp2_mid2(4,4,fill::zeros);  PSp2_mid2=M*P_pos_neg;

  arma::mat pi_transition(4,(S-2)*4,fill::zeros);

  if (S>2){

    for (int i=0;i<=S-3;i++){
      pi_transition(0,(i*4))=1;
      pi_transition(1,(i*4))=theta(p,ns(i+1)-1,00)*theta(p,ns(i)-ns(i+1),10)/theta(p,ns(i)-1,10);
      pi_transition(1,(i*4+1))=1-pi_transition(1,(i*4));
      pi_transition(2,(i*4))=theta(p,ns(i+1)-1,00)*theta(p,ns(i)-ns(i+1),01)/theta(p,ns(i)-1,01);
      pi_transition(2,(i*4+2))=1-pi_transition(2,(i*4));
      pi_transition(3,(i*4))=theta(p,ns(i+1)-1,00)*theta(p,ns(i)-ns(i+1),11)/theta(p,ns(i)-1,11);
      pi_transition(3,(i*4+1))=theta(p,ns(i+1)-1,10)*(theta(p,ns(i)-ns(i+1),01)+theta(p,ns(i)-ns(i+1),11))/theta(p,ns(i)-1,11);
      pi_transition(3,(i*4+2))=theta(p,ns(i+1)-1,01)*(theta(p,ns(i)-ns(i+1),10)+theta(p,ns(i)-ns(i+1),11))/theta(p,ns(i)-1,11);
      pi_transition(3,(i*4+3))=1-pi_transition(3,(i*4))-pi_transition(3,(i*4+1))-pi_transition(3,(i*4+2));
    }

    for (int s=0; s<=S-3;s++){
      se1=SE(s+1,0); se2=SE(s+1,1); sp1=SP(s+1,0); sp2=SP(s+1,1);
      arma::mat P_operator(4,4,fill::zeros);P_operator(0,0)= 1-sp1*sp2; P_operator(1,1)= 1-(1-se1)*sp2;P_operator(2,2)= 1-sp1*(1-se2);P_operator(3,3)= 1-(1-se1)*(1-se2);
      arma::mat P_pos_neg_temp(4,4,fill::zeros); P_pos_neg(0,0)=1-(1-se1)*sp2; P_pos_neg(1,1)=1-(1-se1)*sp2; P_pos_neg(2,2)=1-(1-se1)*(1-se2); P_pos_neg(3,3)=1-(1-se1)*(1-se2);
      arma::mat P_neg_pos_temp(4,4,fill::zeros); P_neg_pos(0,0)=1-sp1*(1-se2); P_neg_pos(1,1)=1-(1-se1)*(1-se2); P_neg_pos(2,2)=1-sp1*(1-se2); P_neg_pos(3,3)=1-(1-se1)*(1-se2);
      arma::mat P_pos_pos_temp(4,4,fill::zeros); P_pos_pos.diag().fill(1-(1-se1)*(1-se2));

      PSe1_mid1=PSe1_mid1*pi_transition(span(0,3),span(s*4,s*4+3))*P_pos_neg;
      PSe1_mid2=PSe1_mid2*pi_transition(span(0,3),span(s*4,s*4+3))*P_pos_pos;
      PSe2_mid1=PSe2_mid1*pi_transition(span(0,3),span(s*4,s*4+3))*P_neg_pos;
      PSe2_mid2=PSe2_mid2*pi_transition(span(0,3),span(s*4,s*4+3))*P_pos_pos;
      PSp1_mid1=PSp1_mid1*pi_transition(span(0,3),span(s*4,s*4+3))*P_operator;
      PSp1_mid2=PSp1_mid2*pi_transition(span(0,3),span(s*4,s*4+3))*P_neg_pos;
      PSp2_mid1=PSp2_mid1*pi_transition(span(0,3),span(s*4,s*4+3))*P_operator;
      PSp2_mid2=PSp2_mid2*pi_transition(span(0,3),span(s*4,s*4+3))*P_pos_neg;
    }
  }

  se1=SE(S-1,0); se2=SE(S-1,1); sp1=SP(S-1,0); sp2=SP(S-1,1);
  double PSE1=p(1)/(p(1)+p(3))*accu(PSe1_mid1)*se1+p(3)/(p(1)+p(3))*accu(PSe1_mid2)*se1;
  double PSE2=p(2)/(p(2)+p(3))*accu(PSe2_mid1)*se2+p(3)/(p(2)+p(3))*accu(PSe2_mid2)*se2;
  double PSP1=1-p(0)/(p(2)+p(0))*accu(PSp1_mid1)*(1-sp1)-p(2)/(p(2)+p(0))*accu(PSp1_mid2)*(1-sp1);
  double PSP2=1-p(0)/(p(1)+p(0))*accu(PSp2_mid1)*(1-sp2)-p(1)/(p(1)+p(0))*accu(PSp2_mid2)*(1-sp2);
  double PPV1=(p(1)+p(3))*PSE1/(PSE1*(p(1)+p(3))+(p(0)+p(2))*(1-PSP1));
  double PPV2=(p(2)+p(3))*PSE2/(PSE2*(p(2)+p(3))+(p(0)+p(1))*(1-PSP2));
  double NPV1=(p(2)+p(0))*PSP1/(PSP1*(p(2)+p(0))+(p(1)+p(3))*(1-PSE1));
  double NPV2=(p(1)+p(0))*PSP2/(PSP2*(p(1)+p(0))+(p(2)+p(3))*(1-PSE2));

  double c00=p(0)*(1-accu(PSp1_mid1)* (1-sp1*sp2));
  double c10=p(1)*accu(PSp2_mid2)*se1*sp2;
  double c01=p(2)*accu(PSe2_mid1)*sp1*se2;
  double c11=p(3)*accu(PSe1_mid2)*se1*se2;

  double res = c00+c10+c01+c11;

  return Rcpp::List::create(Rcpp::Named("PSE1")=PSE1, Rcpp::Named("PSE2")=PSE2,
                            Rcpp::Named("PSP1")=PSP1, Rcpp::Named("PSP2")=PSP2,
                            Rcpp::Named("PPV1")=PPV1, Rcpp::Named("PPV2")=PPV2,
                            Rcpp::Named("NPV1")=NPV1, Rcpp::Named("NPV2")=NPV2,
                            Rcpp::Named("ACCU")=res);
}

//Function to calculate optimal configuration for S=2,3,4,5,6 stages algorithm
// [[Rcpp::export]]
List OPT(arma::vec p, arma::mat SE, arma::mat SP,int MAMPS,std::string obj,std::string algorithm){
  int S_temp = 0;
  double eps=0.00001;
  if (MAMPS<4-eps) {S_temp=2;}
  else if (MAMPS<8-eps) {S_temp=3;}
  else if (MAMPS<16-eps) {S_temp=4;}
  else if (MAMPS<32-eps) {S_temp=5;}
  else S_temp=6;

  arma::mat efficiency(S_temp-1,10,fill::zeros); arma::mat ps(S_temp-1,6,fill::zeros);

  for (int S=2; S<=S_temp; S++){
    arma::mat pools = poolsize(MAMPS,S,algorithm);
    int num_pools = pools.n_rows;
    arma::mat temp(num_pools,2,fill::zeros);
    for (int i=0; i<num_pools; i++){
      temp(i,0)=EFF_HIER(p,S,SE(span((S+1)*(S-2)/2,(S+1)*(S-2)/2+S-1),span::all),SP(span((S+1)*(S-2)/2,(S+1)*(S-2)/2+S-1),span::all),trans(pools(i,span::all)));
      temp(i,1)=ACCU_HIER(p,S,SE(span((S+1)*(S-2)/2,(S+1)*(S-2)/2+S-1),span::all),SP(span((S+1)*(S-2)/2,(S+1)*(S-2)/2+S-1),span::all),trans(pools(i,span::all)))["ACCU"];
    }

    if (obj=="minimize"){
      double mid=min(temp(span::all,0));
      arma::uvec q1=find(temp(span::all,0)==mid);
      double q1d=as_scalar(q1);
      List mid1=ACCU_HIER(p,S,SE(span((S+1)*(S-2)/2,(S+1)*(S-2)/2+S-1),span::all),SP(span((S+1)*(S-2)/2,(S+1)*(S-2)/2+S-1),span::all),trans(pools(q1d,span::all)));
      efficiency(S-2,0)=mid; efficiency(S-2,1)=mid1["PSE1"]; efficiency(S-2,2)=mid1["PSE2"];  efficiency(S-2,3)=mid1["PSP1"]; efficiency(S-2,4)=mid1["PSP2"];
      efficiency(S-2,5)=mid1["PPV1"]; efficiency(S-2,6)=mid1["PPV2"];efficiency(S-2,7)=mid1["NPV1"]; efficiency(S-2,8)=mid1["NPV2"];efficiency(S-2,9)=mid1["ACCU"];
      ps(S-2,span(0,S-1))=pools(q1d,span::all);
    } else if (obj=="maximize"){
      arma::vec MAR=temp(span::all,1)/temp(span::all,0);
      double mid=max(MAR);
      arma::uvec q1=find(MAR(span::all)==mid);
      double q1d=as_scalar(q1);
      List mid1=ACCU_HIER(p,S,SE(span((S+1)*(S-2)/2,(S+1)*(S-2)/2+S-1),span::all),SP(span((S+1)*(S-2)/2,(S+1)*(S-2)/2+S-1),span::all),trans(pools(q1d,span::all)));
      efficiency(S-2,0)=temp(q1d,0); efficiency(S-2,1)=mid1["PSE1"]; efficiency(S-2,2)=mid1["PSE2"];  efficiency(S-2,3)=mid1["PSP1"]; efficiency(S-2,4)=mid1["PSP2"];
      efficiency(S-2,5)=mid1["PPV1"]; efficiency(S-2,6)=mid1["PPV2"];efficiency(S-2,7)=mid1["NPV1"]; efficiency(S-2,8)=mid1["NPV2"];efficiency(S-2,9)=mid1["ACCU"];
      ps(S-2,span(0,S-1))=pools(q1d,span::all);
    }
  }

  return Rcpp::List::create(Rcpp::Named("OperatingCharacteristics")=efficiency,
                            Rcpp::Named("PoolSize")=ps);
}

//Function to find optimal stage
// [[Rcpp::export]]
List optimal_stage(double rho, arma::mat SE, arma::mat SP, arma::vec pi1,arma::vec pi2, int MAMPS, std::string obj, std::string algorithm){
  int size = pi1.n_elem;
  arma::mat stage(size,size,fill::zeros);
  arma::mat pool2=poolsize(MAMPS,2,algorithm), pool3=poolsize(MAMPS,3,algorithm),pool4=poolsize(MAMPS,4,algorithm),pool5=poolsize(MAMPS,5,algorithm),pool6=poolsize(MAMPS,6,algorithm);
  arma::vec two(pool2.n_rows,fill::ones), three(pool3.n_rows,fill::ones),  four(pool4.n_rows,fill::ones), five(pool5.n_rows,fill::ones), six(pool6.n_rows,fill::ones);
  arma::vec two1(pool2.n_rows,fill::ones), three1(pool3.n_rows,fill::ones),  four1(pool4.n_rows,fill::ones), five1(pool5.n_rows,fill::ones), six1(pool6.n_rows,fill::ones);

  int num_config2=pool2.n_rows, num_config3=pool3.n_rows, num_config4=pool4.n_rows, num_config5=pool5.n_rows, num_config6=pool6.n_rows;

  for (int i=0;i<size;i++)
  {
    for (int j=i;j<size;j++)
    {
      double p11=rho*sqrt(pi1(j)*(1-pi1(j))*pi2(i)*(1-pi2(i)))+pi1(j)*pi2(i);
      double p10=pi1(j)-p11;
      double p01=pi2(i)-p11;
      double p00=1-p11-p10-p01;
      arma::vec p(4,fill::zeros); p(0)=p00; p(1)=p10; p(2)=p01; p(3)=p11;

      arma::vec max_stage(5,fill::zeros);
      arma::vec min_stage(5,fill::zeros);

      if (min(p)>=0 && max(p)<=1 )
      {

        for (int k=0;k<num_config2;k++)
        {
          two(k)=EFF_HIER(p,2,SE(span(0,1),span::all),SP(span(0,1),span::all),trans(pool2(k,span::all)));
          two1(k)=ACCU_HIER(p,2,SE(span(0,1),span::all),SP(span(0,1),span::all),trans(pool2(k,span::all)))["ACCU"];
        }

        min_stage(0)=min(two);
        max_stage(0)=max(two1/two);

        for (int k=0;k<num_config3;k++)
        {
          three(k)=EFF_HIER(p,3,SE(span(2,4),span::all),SP(span(2,4),span::all),trans(pool3(k,span::all)));
          three1(k)=ACCU_HIER(p,3,SE(span(2,4),span::all),SP(span(2,4),span::all),trans(pool3(k,span::all)))["ACCU"];
        }

        max_stage(1)=max(three1/three);
        min_stage(1)=min(three);

        for (int k=0;k<num_config4;k++)
        {
          four(k)=EFF_HIER(p,4,SE(span(5,8),span::all),SP(span(5,8),span::all),trans(pool4(k,span::all)));
          four1(k)=ACCU_HIER(p,4,SE(span(5,8),span::all),SP(span(5,8),span::all),trans(pool4(k,span::all)))["ACCU"];
        }
        max_stage(2)=max(four1/four);
        min_stage(2)=min(four);

        for (int k=0;k<num_config5;k++)
        {
          five(k)=EFF_HIER(p,5,SE(span(9,13),span::all),SP(span(9,13),span::all),trans(pool5(k,span::all)));
          five1(k)=ACCU_HIER(p,5,SE(span(9,13),span::all),SP(span(9,13),span::all),trans(pool5(k,span::all)))["ACCU"];
        }
        max_stage(3)=max(five1/five);
        min_stage(3)=min(five);

        for (int k=0;k<num_config6;k++)
        {
          six(k)=EFF_HIER(p,6,SE(span(14,19),span::all),SP(span(14,19),span::all),trans(pool6(k,span::all)));
          six1(k)=ACCU_HIER(p,6,SE(span(14,19),span::all),SP(span(14,19),span::all),trans(pool6(k,span::all)))["ACCU"];
        }
        max_stage(4)=max(six1/six);
        min_stage(4)=min(six);

        if (obj=="minimize"){
          double min_all=min(min_stage);
          if (min_all>1){
            stage(i,j)=1; stage(j,i)=1;
          } else {
            arma::uvec q1=find(min_stage(span::all)==min_all);
            double q1d=as_scalar(q1);
            stage(i,j)=q1d+2; stage(j,i)=stage(i,j);
          }
        } else if (obj=="maximize"){
          double max_all=max(max_stage);
          if (max_all<p00*SP(1,0)*SP(1,1)+p10*SE(1,0)*SP(1,1)+p01*SP(1,0)*SE(1,1)+p11*SE(1,0)*SE(1,1)){
            stage(i,j)=1; stage(j,i)=stage(i,j);
          } else{
            arma::uvec q1=find(max_stage(span::all)==max_all);
            double q1d=as_scalar(q1);
            stage(i,j)=q1d+2; stage(j,i)=stage(i,j);
          }
        }
      }
      else {
        stage(i,j)=0;
        stage(j,i)=0;
      }

    }
  }
  return Rcpp::List::create(Rcpp::Named("stage")=stage);
}
