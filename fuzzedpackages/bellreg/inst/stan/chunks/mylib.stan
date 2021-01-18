
functions{

    real lambertW(real x){
      real y;
      real w;
      y= sqrt( 1 + exp(1) *x);
      w=-1 + 2.036 * log( (1 + 1.14956131 * y)/(1 + 0.45495740*log(1+y)) );
      w = (w/(1+w)) * ( 1 + log(x/w) );
      w = (w/(1+w)) * ( 1 + log(x/w) );
      w = (w/(1+w)) * ( 1 + log(x/w) );
      return(w);
    }

   real bellnumber(int n){
    if(n < 2){
      return(1);
    }else{
      int k;
      vector[n] B;
      vector[n] Bneu;
      B[1] = 1;
      for (i in 1:(n - 1)){
        k = i;
        Bneu[1] = B[i];
        for (j in 2:(i + 1)){
          Bneu[j] = B[j - 1] + Bneu[j - 1];
        }
        for(j in 1:n){
          B[j] = Bneu[j];
        }
      }
      return(Bneu[k + 1]);
    }
  }


  real bell_lpmf(int x, real theta){
      real Bx;
      real lprob;
      Bx = bellnumber(x);
      lprob = x*log(theta) - exp(theta) + 1 + log(Bx) - lgamma(x+1);
      return lprob;
    }
  
  // real bell_lpmf(int[] x, real theta){
  //     int n = num_elements(x);
  //     real Bx[n];
  //     real lprob[n];
  //     for(i in 1:n){
  //       Bx[i] = bellnumber(x[i]);
  //       lprob[i] = x[i]*log(theta) - exp(theta) + 1 + log(Bx[i]) - lgamma(x[i]+1);  
  //     }
  //     return sum(lprob);
  // }  

  real loglik_bell(int[] x, real[] theta){
      real lprob[num_elements(x)];
      for(i in 1:num_elements(x)){
        lprob[i] = x[i]*log(theta[i]) - exp(theta[i]);
      }
      return sum(lprob);
    }

}

