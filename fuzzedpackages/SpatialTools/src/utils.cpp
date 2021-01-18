#include "utils.h"

using namespace Rcpp;
using namespace arma;

arma::mat rcondsim(int nsim, arma::vec y, arma::mat w, arma::mat Vediag, arma::mat dV, int method, double m)
{
    int n = y.n_elem;
	int na = dV.n_rows;
	arma::mat condsim = arma::mat(na, nsim);

	RNGScope scope;
	
	NumericVector Zr = NumericVector(rnorm((na + n) * nsim, 0, 1));
	arma::mat Z(Zr.begin(), (na + n), nsim);
    
	arma::mat newsim = dV * Z.rows(0, na - 1);
  
    // Conditional simulation algorithm changes slightly for simple kriging
    // when m != 0.  
    if(m != 0)
    {    
        arma::mat newZ = repmat(y - m, 1, nsim) - newsim.rows(0, n - 1) + diagmat(sqrt(Vediag)) * Z.rows(na, na + n - 1);
    
        //accounts for the fact that indexing starts at zero
        condsim = m + newsim.rows(n, na - 1) + trans(w) * newZ;
    }
    else
    {    
        arma::mat newZ = repmat(y, 1, nsim) - newsim.rows(0, n - 1) + diagmat(sqrt(Vediag)) * Z.rows(na, na + n - 1);
        
        //accounts for the fact that indexing starts at zero
        condsim = newsim.rows(n, na - 1) + trans(w) * newZ;
    }

    return condsim;
}

arma::mat decomp_V(const arma::mat &V, int method)
{
	int n = V.n_rows;
    
	arma::mat dV = arma::mat(n, n);
	
	if(method == 1)
    {
        arma::vec eigval = arma::vec(n);
        arma::mat eigvec = arma::mat(n, n);
        
        //compute eigen values and vectors of V
        eig_sym(eigval, eigvec, V);
        
        for(int i = 0; i < eigval.n_rows; i++)
        {
            if(eigval(i) < 0)
            {
                eigval(i) = 0;
            }
        }
        
        dV = eigvec * diagmat(sqrt(eigval));
    }
    else if(method == 2)
    {
        dV = trans(arma::chol(V));
    }
    else
    {
        arma::mat U = arma::mat(n, n);
        arma::mat U2 = arma::mat(n, n);
        arma::vec sv = arma::vec(n);
        
        svd_econ(U, sv, U2, V);
        
        dV = U * diagmat(sqrt(sv)) * trans(U2);
    }
    
    return dV;
}

arma::mat rmvnorm(int nsim, const arma::mat &mu, const arma::mat &V, int method){
    
	arma::mat dV = decomp_V(V, method);
	
	RNGScope scope;
	NumericVector Zr = NumericVector(rnorm(V.n_rows * nsim, 0, 1));
	arma::mat Z(Zr.begin(), V.n_rows, nsim);
	
	return repmat(mu, 1, nsim) + dV * Z;
}

arma::mat rcondnorm(int nsim, const arma::mat &y,
                    const arma::mat &mu, const arma::mat &mup,
                    const arma::mat &V, const arma::mat &Vp, const arma::mat &Vop, int method)
{
	arma::mat ViVop = solve(V, Vop);
    
	arma::mat mc = mup + trans(ViVop) * (y - mu);
	arma::mat Vc = Vp - trans(Vop) * ViVop;
    
	arma::mat sim = rmvnorm(nsim, mc, Vc, method);
	
	return sim;
}

arma::mat dist1(const arma::mat &coords)
{
	int n = coords.n_rows;
	arma::mat D = arma::zeros(n, n);
    
	for(int j = 1; j < n; j++)
	{
		for(int i = 0; i <=j ; i++)
		{
			D(i, j)=sqrt(pow(coords(i,0)-coords(j,0),2)+pow(coords(i,1)-coords(j,1),2));
			D(j, i) = D(i, j);
		}
	}
	return D;
}

arma::mat dist2(const arma::mat &coords, const arma::mat &pcoords)
{
	int n = coords.n_rows;
	int np = pcoords.n_rows;
    
	arma::mat D = arma::zeros(n, np);
    
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < np; j++)
		{
			D(i, j) = sqrt(pow(coords(i,0)-pcoords(j,0),2) +
                           pow(coords(i,1)-pcoords(j,1),2));
		}
	}
	return D;
}

arma::mat cov_spBayes(const arma::mat &D, int sp_type, double sigmasq,
                      double phi, double nu, double ev, double fv)
{
	int nr = D.n_rows;
	int nc = D.n_cols;
    
	arma::mat V = arma::zeros(nr, nc);
    
	for(int i = 0; i < nr; i++)
	{
		for(int j = 0; j < nc; j++)
		{
			if(D(i, j) == 0)
			{
				V(i, j) = sigmasq + ev + fv;
			}
			else
			{
				if(sp_type == 1)
				{
					V(i, j) = sigmasq * exp(-D(i, j) * phi);
				}
				else if(sp_type == 2)
				{
					V(i, j) = sigmasq * exp(-pow(D(i, j) * phi, 2));
				}
				else if(sp_type == 3)
				{
					V(i, j) = sigmasq * pow(D(i,j)*phi, nu)/(pow(2, nu-1)*R::gammafn(nu))*R::bessel_k(D(i,j)*phi, nu, 1.0);
				}
				else
				{
					if(D(i, j) <= 1.0/phi)
					{
						V(i, j) = sigmasq * (1.0 - 1.5*phi*D(i,j) + 0.5*pow(phi*D(i,j),3));
					}
					
				}
			}
		}
	}	
	return V;
}
