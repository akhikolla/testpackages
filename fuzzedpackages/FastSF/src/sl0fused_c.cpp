#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List sl0fused_c(arma::vec y, int T0, int T02, int max_steps)
{
  int n=y.n_elem;
  double rho = n*n;
  int m = n-1;
  
  arma::vec beta = y;
  
  arma::vec d=y;
  arma::vec z=diff(y);

  arma::vec u(m);
  u.zeros();
  
  


  arma::uvec A(T0);
  A=A.fill(0);
  arma::uvec A2(T0);
  A2=A2.fill(0);

  int k=0;

  arma::vec w(n);
  w.fill(rho*2);
  w(0)=rho;
  w(n-1)=rho;
  w = 1+w;




  while(k < max_steps)
  {
    // cout<<k<<endl;
    arma::vec bd=abs(z+u/rho);
    arma::uvec bd_order=sort_index(bd,"descend");
    arma::uvec A_new=bd_order.head(T0);
    arma::uvec I_new=bd_order.tail(n-1-T0);
    A_new = sort(A_new);
    I_new = sort(I_new);

    arma::vec bd2=sqrt(w)%abs(beta+d);
    arma::uvec bd2_order=sort_index(bd2,"descend");
    arma::uvec A2_new=bd2_order.head(T02);
    arma::uvec I2_new=bd2_order.tail(n-T02);

    A2_new = sort(A2_new);
    I2_new = sort(I2_new);

    beta=beta.fill(0);
    d=d.fill(0);

    arma::vec dtu(u.n_elem+1);
    dtu(0)=-u(0);
    dtu(u.n_elem)=u(u.n_elem-1);
    arma::uvec u_index=arma::linspace<arma::uvec>(1,u.n_elem-1,u.n_elem-1);
    dtu(u_index) = -diff(u);

    beta(A2_new)=y(A2_new)-dtu(A2_new);


    arma::vec dtz(z.n_elem+1);
    dtz(0)=-z(0);
    dtz(z.n_elem)=z(z.n_elem-1);
    arma::uvec z_index=arma::linspace<arma::uvec>(1,z.n_elem-1,z.n_elem-1);
    dtz(z_index) = -diff(z);

    arma::vec dbeta = diff(beta);

    arma::vec dtdbeta(dbeta.n_elem+1);
    dtdbeta(0)=-dbeta(0);
    dtdbeta(dbeta.n_elem)=dbeta(dbeta.n_elem-1);
    arma::uvec dbeta_index=arma::linspace<arma::uvec>(1,dbeta.n_elem-1,dbeta.n_elem-1);
    dtdbeta(dbeta_index) = -diff(dbeta);

    arma::vec dtemp=((y-dtu+rho*dtz-beta-rho*dtdbeta)/w);

    d(I2_new)=dtemp(I2_new);


    z=z.fill(0);
    u=u.fill(0);

    arma::vec dy = diff(y);


    arma::uvec dI_new = diff(I_new);
    arma::uvec cp = find(dI_new!=1);






    if (cp.n_elem == 0)
    {
      arma::vec y_temp = diff(y);
      arma::vec x_col = arma::linspace<arma::vec>(1,I_new.n_elem,I_new.n_elem);

      arma::uvec I_new_index(1);
      I_new_index = 1;
      while(I_new_index(0)<I_new.n_elem+1)
      {
        arma::vec x_row(I_new.n_elem);
        x_row.fill(I_new_index(0));
        arma::vec x_min;
        arma::vec x_max;
        x_min = arma::min(x_row, x_col);
        x_max = arma::max(x_row, x_col);
        arma::vec di_inv = x_min % (I_new.n_elem + 1 - x_max)/(I_new.n_elem + 1);
        arma::vec y_I_temp;
        y_I_temp = y_temp(I_new);
        u(I_new(I_new_index(0)-1)) = dot(di_inv, y_I_temp);
        I_new_index = I_new_index + 1;
      }
    }
    else
    {
      arma::uvec tail_temp(1);
      tail_temp = I_new.n_elem-1;
      arma::uvec index_temp = arma::linspace<arma::uvec>(1,cp.n_elem, cp.n_elem);

      arma::uvec cp_new;
      cp_new = cp_new.zeros(cp.n_elem+1);

      cp_new(index_temp-1)=cp;

      arma::uvec index_temp2;
      index_temp2  = cp.n_elem;
      cp_new(index_temp2) = tail_temp;
      cp=cp_new;

      arma::uvec block_index;
      block_index=1;
      while(block_index(0)<cp.n_elem+1)
      {
        if (block_index(0)==1){
          arma::uvec I_temp;
          arma::uvec index_temp3;
          index_temp3 = arma::linspace<arma::uvec>(1,cp(block_index(0)-1)+1,cp(block_index(0)-1)+1);
          index_temp3 = index_temp3-1;

          I_temp = I_new(index_temp3);

          arma::vec y_temp = diff(y);
          arma::vec x_col = arma::linspace<arma::vec>(1,I_temp.n_elem,I_temp.n_elem);

          arma::uvec I_temp_index(1);
          I_temp_index = 1;
          while(I_temp_index(0)<I_temp.n_elem+1)
          {
            arma::vec x_row(I_temp.n_elem);
            x_row.fill(I_temp_index(0));
            arma::vec x_min;
            arma::vec x_max;
            x_min = arma::min(x_row, x_col);
            x_max = arma::max(x_row, x_col);
            arma::vec di_inv = x_min % (I_temp.n_elem + 1 - x_max)/(I_temp.n_elem + 1);
            arma::vec y_I_temp;
            y_I_temp = y_temp(I_temp);
            u(I_temp(I_temp_index(0)-1)) = dot(di_inv, y_I_temp);
            I_temp_index = I_temp_index + 1;
          }



        }


        else{
          arma::uvec I_temp;
          arma::uvec index_temp3;
          index_temp3 = arma::linspace<arma::uvec>(cp(block_index(0)-2)+2,cp(block_index(0)-1)+1,cp(block_index(0)-1)-cp(block_index(0)-2));
          index_temp3 = index_temp3-1;

          I_temp = I_new(index_temp3);
          arma::vec y_temp = diff(y);
          arma::vec x_col = arma::linspace<arma::vec>(1,I_temp.n_elem,I_temp.n_elem);

          arma::uvec I_temp_index(1);
          I_temp_index = 1;
          while(I_temp_index(0)<I_temp.n_elem+1)
          {
            arma::vec x_row(I_temp.n_elem);
            x_row.fill(I_temp_index(0));
            arma::vec x_min;
            arma::vec x_max;
            x_min = arma::min(x_row, x_col);
            x_max = arma::max(x_row, x_col);
            arma::vec di_inv = x_min % (I_temp.n_elem + 1 - x_max)/(I_temp.n_elem + 1);
            arma::vec y_I_temp;
            y_I_temp = y_temp(I_temp);
            u(I_temp(I_temp_index(0)-1)) = dot(di_inv, y_I_temp);
            I_temp_index = I_temp_index + 1;
          }

        }


        block_index = block_index+1;

      }


    }






    dtu(0)=-u(0);
    dtu(u.n_elem)=u(u.n_elem-1);
    dtu(u_index) = -diff(u);




    arma::vec z_temp = diff(y-dtu);
    z(A_new) = z_temp(A_new);


    A_new=sort(A_new);
    A=sort(A);

    arma::uvec if_true=A_new==A;
    if(sum(if_true)==A.n_elem) break;
    else{
      k=k+1;
      A=A_new;
      A2=A2_new;
    }
  }
  return List::create(Named("beta")=beta,Named("d")=d,Named("z")=z,Named("u")=u);
}
