#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List l0fused_c(arma::vec y, int T0, int max_steps)
{
  arma::vec beta = y;
  int n = y.n_elem;
  int m = n-1;
  
  double rho = n*n;
  
  arma::vec z = diff(y);
  arma::vec u(m);
  u.zeros();
  
  
  arma::uvec A(T0);
  A = A.fill(0);
 
  
  int k = 0;
  
  
  while(k < max_steps)
  {
    arma::vec bd = abs(z+u/rho);
    arma::uvec bd_order = sort_index(bd,"descend");
    arma::uvec A_new = bd_order.head(T0);
    arma::uvec I_new = bd_order.tail(m-T0);
    A_new = sort(A_new);
    I_new = sort(I_new);
    
    
    z=z.zeros(m);
    u=z.zeros(m);
    
    arma::uvec index = diff(I_new);
    arma::uvec cp = find(index!=1);
    
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
    
    
    arma::vec u_diff=-diff(u);
    arma::vec u_head;
    u_head = -u(0);
    arma::vec u_tail;
    u_tail = u(u.n_elem-1);
    arma::vec u_temp(u_diff.n_elem+2);
    arma::uvec u_index;
    arma::uvec u_head_index;
    u_head_index = 1;
    
    u_temp(u_head_index-1) = u_head;
    
    u_index = arma::linspace<arma::uvec>(1,u_diff.n_elem,u_diff.n_elem);
    u_temp(u_index) = u_diff;
    arma::uvec u_tail_index;
    u_tail_index = u_diff.n_elem+1;
    
    u_temp(u_tail_index) = u_tail;
    u_temp = y-u_temp;
    beta  = u_temp;
    arma::vec z_temp = diff(u_temp);
    z(A_new) = z_temp(A_new);
    
    A_new=sort(A_new);
    A=sort(A);
    arma::uvec if_true=A_new==A;
    if(sum(if_true)==A.n_elem) break;
    else{
      k=k+1;
      A=A_new;
    }
    
  }
  
  
  return List::create(Named("beta")=beta,Named("z")=z,Named("u")=u);
  
}
