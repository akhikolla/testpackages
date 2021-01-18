
#include <RcppArmadillo.h>

#define ALL 0
#define LAST 1

using namespace arma;


template< class T>
void print_vec(const std::vector<T>& v)
{
	const int n=v.size();
	int r=0;
	printf("[1]   ");
	for(int i=0;i<n;i++)
	{
		printf("%012.3f     ", (double)v[i]);
		r++;
		if(r==7)
		{
			r=0;
			printf("\n[%i]  ",i);
		}
	}
	printf("\n");
	fflush(stdout);
}



int elsym(const int routing, const vec& b, const ivec& a, 
			int* first_ptr, int* last_ptr, const int nI,
			const int* mod_min_ptr, const int* mod_max_ptr, const int nmod,
			int* nit_ptr, 
			std::vector<long double>& g_out, const int item1_first=-1, const int aij=1, 
			const int item2_first=-1, const int akl=1)
{
	// trick to get lvalued vector views
	const ivec first(first_ptr, nI, false, true);
	const ivec last(last_ptr, nI, false, true);
	const ivec nit(nit_ptr, nmod, false, true);	
	
	// these do not use aux mem as they can be changed
	ivec mod_min(mod_min_ptr, nmod);
	ivec mod_max(mod_max_ptr, nmod);
	
	ivec cnit(nmod+1);
	cnit[0] = 0;
	for(int m=0; m<nmod; m++)
		cnit[m+1] = cnit[m] + nit[m];	
	
	std::vector<long double> g(g_out.size(), 0); 
	  
	std::fill(g_out.begin(), g_out.end(), 0); 


	if(routing == LAST)
	{
		int cummin = 0, cummax = 0;
		
		for(int m=0; m < nmod; m++)
		{
			std::fill(g.begin(), g.end(), 0); 
			g[0] = 1;
			int Msc = 0;
			
			for(int i = cnit[m]; i<cnit[m+1]; i++)
			{
				if(first[i] == item1_first)
				{
					mod_max[m] -= aij; //a[last[i]];
					mod_min[m] = std::max(mod_min[m] - aij, 0);
				}
				if(first[i] == item2_first)
				{
				    mod_max[m] -= akl; //a[last[i]];
				    mod_min[m] = std::max(mod_min[m] - akl, 0);
				}
			}
			if(mod_min[m] >= mod_max[m])
				return -1;
			
			
			for (int i=cnit[m]; i < cnit[m+1]; i++)
			{
				if(first[i] != item1_first && first[i] != item2_first)
				{
					for (int s=Msc; s>=0;s--)
						for (int j=last[i];j>=first[i];j--)
							g[s+a[j]] += g[s]*b[j];

					Msc+=a[last[i]];
					Msc = std::min(Msc, mod_max[m]); // a little extra speed
				} 			
			}
			// merge
			if(m==0)
			{
				for(int s=mod_min[0]; s<=mod_max[0]; s++)
					g_out[s] = g[s];			
			}
			else
			{
				for(int s1 = cummax; s1>=cummin; s1--)
					for(int s2 = mod_min[m]; s2<=mod_max[m]; s2++)
						if(s2>0)
							g_out[s1+s2] += g_out[s1] * g[s2];			
			}		
			cummin += mod_min[m];
			cummax += mod_max[m];
		}		
		for(int s=0;s<cummin;s++)
			g_out[s] = 0;
	}
	else //all
	{
		for(int m=0; m < nmod; m++)
		{
			std::fill(g.begin(), g.end(), 0); 
			g[0] = 1;
			int Msc = 0;
			
			for(int i = cnit[m]; i<cnit[m+1]; i++)
			{
				if(first[i] == item1_first)
				{
					for(int nm = m; nm < nmod; nm++)
					{
						mod_max[nm] -= aij;
						mod_min[nm] = std::max(mod_min[nm] - aij, 0); 
					}
				}
				if(first[i] == item2_first)
				{
					for(int nm = m; nm < nmod; nm++)
				    {
						mod_max[nm] -= akl;
						mod_min[nm] = std::max(mod_min[nm] - akl, 0); 
				    }
				}
			}
			if(mod_min[m] >= mod_max[m])
				return -1;
			
			for (int i=cnit[m]; i < cnit[m+1]; i++)
			{
				if(first[i] != item1_first && first[i] != item2_first)
				{
					for (int s=Msc; s>=0;s--)
						for (int j=last[i];j>=first[i];j--)
							g[s+a[j]] += g[s]*b[j];

					Msc+=a[last[i]];
					Msc = std::min(Msc, mod_max[m]); // a little extra speed
				} 			
			}
				
			// merge
			if(m==0)
			{
				for(int s=mod_min[0]; s<=mod_max[0]; s++)
					g_out[s] = g[s];			
			}
			else
			{
				for(int s1 = mod_max[m-1]; s1>=mod_min[m-1]; s1--)
					for(int s2 = mod_max[m]-s1; s2>=mod_min[m]-s1; s2--) 
						if(s2 > 0 )
							g_out[s1+s2] += g_out[s1] * g[s2];
			}		
		}
		for(int s=0;s<mod_min[nmod-1];s++)
			g_out[s] = 0;
	}
	return 0;
}


// [[Rcpp::export]]
void Expect( const arma::vec& b, const arma::ivec& a, 
             arma::ivec& bfirst, arma::ivec& blast, const arma::ivec& bmax,
			 const arma::ivec& nmod,  const arma::ivec& brouting,
			 arma::ivec& mnit, const arma::ivec& mod_min, const arma::ivec& mod_max,
			 const arma::ivec& scoretab,  /* out */ arma::vec& E)
{
	const int nb = nmod.n_elem;
	const int len_g = max(bmax) + 1;
	
	// bookkeeping
	ivec cbmax(nb+1), bnit(nb, fill::zeros), cbnit(nb+1), cbmod(nb+1);
	cbmax[0] = 0, cbnit[0] = 0, cbmod[0] = 0;
	{
		int i=0;
		for(int bi=0; bi<nb; bi++)
		{			
			for(int m=0; m<nmod[bi]; m++)
				bnit[bi] += mnit[i++];
			cbnit[bi+1] = cbnit[bi] + bnit[bi];
			cbmax[bi+1] = cbmax[bi] + bmax[bi]+1;
			cbmod[bi+1] = cbmod[bi] + nmod[bi];

		}	
	}
	

	E.zeros();
	std::vector<long double> g(len_g), gi(len_g);

	for(int bi=0; bi<nb; bi++)
	{
		elsym(brouting[bi], b, a,  bfirst.memptr() + cbnit[bi], blast.memptr() + cbnit[bi], bnit[bi], 
				mod_min.memptr() + cbmod[bi], mod_max.memptr() + cbmod[bi], nmod[bi],
				mnit.memptr() + cbmod[bi], g);
		
		for(int i=cbnit[bi]; i<cbnit[bi+1]; i++)
		{
			for(int j=bfirst[i]; j<= blast[i]; j++)
			{
			  elsym(brouting[bi], b, a, bfirst.memptr() + cbnit[bi], blast.memptr() + cbnit[bi], bnit[bi], 
           mod_min.memptr() + cbmod[bi], mod_max.memptr() + cbmod[bi], nmod[bi],
            mnit.memptr() + cbmod[bi],  gi, bfirst[i], a[j]);
			  
				for (int s = a[j]; s <= bmax[bi]; s++) 
				{
				  if(g[s]>0)
				  {
					  E[j] += scoretab[cbmax[bi]+s] * b[j] * (gi[s-a[j]]/g[s]); 
				  }
				}
			}		
		}
	}

}


// [[Rcpp::export]]
void NR( const arma::vec& b, const arma::ivec& a, 
             arma::ivec& bfirst, arma::ivec& blast, const arma::ivec& bmax,
			 const arma::ivec& nmod,  const arma::ivec& brouting,
			 arma::ivec& mnit, const arma::ivec& mod_min, const arma::ivec& mod_max,
			 const arma::ivec& scoretab,  /* out */ arma::vec& E, arma::mat& H)
{
	const int nb = nmod.n_elem;
	const int len_g = max(bmax) + 1;// if crashes + maxa
	double tmp;
	
	// bookkeeping
	ivec cbmax(nb+1), bnit(nb, fill::zeros), cbnit(nb+1), cbmod(nb+1);
	cbmax[0] = 0, cbnit[0] = 0, cbmod[0] = 0;
	{
		int i=0;
		for(int bi=0; bi<nb; bi++)
		{
			for(int m=0; m<nmod[bi]; m++)
				bnit[bi] += mnit[i++];
			cbnit[bi+1] = cbnit[bi] + bnit[bi];
			cbmax[bi+1] = cbmax[bi] + bmax[bi]+1;
			cbmod[bi+1] = cbmod[bi] + nmod[bi];
		}	
	}	
	
	E.zeros();
	H.zeros();
	std::vector<long double> g(len_g), gi(len_g), gk(len_g), gik(len_g);
	vec cc(len_g, fill::zeros);
		
	for(int bi=0; bi<nb; bi++)
	{
		elsym(brouting[bi],b, a, bfirst.memptr() + cbnit[bi], blast.memptr() + cbnit[bi], bnit[bi], 
				mod_min.memptr() + cbmod[bi], mod_max.memptr() + cbmod[bi], nmod[bi],
				mnit.memptr() + cbmod[bi], g);

		for(int i=cbnit[bi]; i<cbnit[bi+1]; i++)
		{
			for (int j=bfirst[i]; j<=blast[i]; j++)
			{
				elsym(brouting[bi], b, a,  bfirst.memptr() + cbnit[bi], blast.memptr() + cbnit[bi], bnit[bi], 
						mod_min.memptr() + cbmod[bi], mod_max.memptr() + cbmod[bi], nmod[bi],
						mnit.memptr() + cbmod[bi], gi, bfirst[i], a[j]);
				
				for (int s = a[j]; s <= bmax[bi]; s++)
				{
					if(g[s]>0)
					{
						tmp = b[j] * (gi[s-a[j]]/g[s]);
						cc[s] = scoretab[cbmax[bi]+s] * tmp;
						E[j] += cc[s]; 
						H.at(j,j) += cc[s] * (1-tmp);
					}
				}
				
				// between categories of the same item
				for (int k = (j+1); k <= blast[i]; k++)
				{
					for (int s = a[k]; s <= bmax[bi]; s++)
					{
						if(g[s]>0)
							H.at(k,j) -= cc[s] * b[k] * (gi[s-a[k]]/g[s]);						
					}
					H.at(j,k) = H.at(k,j);
				}
				
				for (int k=i+1; k < cbnit[bi+1]; k++)
				{
					for (int l=bfirst[k]; l<=blast[k];l++)
					{
						elsym(brouting[bi], b, a,  bfirst.memptr() + cbnit[bi], blast.memptr() + cbnit[bi], bnit[bi], 
								mod_min.memptr() + cbmod[bi], mod_max.memptr() + cbmod[bi], nmod[bi],
								mnit.memptr() + cbmod[bi], gk, bfirst[k],a[l]);
					  
						int success = elsym(brouting[bi], b, a,  bfirst.memptr() + cbnit[bi], blast.memptr() + cbnit[bi], bnit[bi], 
											mod_min.memptr() + cbmod[bi], mod_max.memptr() + cbmod[bi], nmod[bi],
											mnit.memptr() + cbmod[bi], gik, bfirst[i], a[j], bfirst[k], a[l]);
						if(success>=0)
						for (int s=0;s<=bmax[bi];s++)
						{
							if (g[s]>0)
							{
								if (s >= a[j]+a[l]) 
								{
									H.at(l,j) +=  scoretab[cbmax[bi]+s] * (gik[s-a[j]-a[l]]/g[s]) * b[j]*b[l];
								}
								if (s >= a[j] && s >= a[l]) // hier gaat em mis waarschijnlijk 
								{
									H.at(l,j) -= cc[s] * b[l] * (gk[s-a[l]] / g[s]);//scoretab[cbmax[bi]+s] * b[j] * b[l] * ((gi[s-a[j]]/g[s])*(gk[s-a[l]]/g[s]));
								}
							}
						}
						H.at(j,l) = H.at(l,j);
					}
				}
			}		
		}
	}
	// there seems to be no sum of H with it's transpose like in dexter, is this caused by omission of 0 cat?
}

// [[Rcpp::export]]
void dirichlet(const arma::vec& alpha, arma::vec& out)
{
	const int n = alpha.n_elem;
	
	out.zeros();
	for(int i=0; i<n; i++)
		out[i] = R::rgamma(alpha[i], 1);	

	out /= accu(out);
}

// to do min scores bmin kunnen er uit
// [[Rcpp::export]]
arma::mat calibrate_Bayes(const arma::ivec& a, const arma::ivec& first, const arma::ivec& last, 
							arma::ivec& bfirst, arma::ivec& blast, const arma::ivec& bmax, const arma::ivec& bmin,
							const arma::ivec& nmod,  const arma::ivec& brouting,
							arma::ivec& mnit, const arma::ivec& mod_min, const arma::ivec& mod_max,
						    const arma::ivec& itb, const arma::ivec& itnb,
							const arma::ivec& sufI, const arma::ivec& scoretab,
							arma::vec& b,  const arma::vec& fixed_b, 
							const int from, const int step, const int ndraws,
							const double prior_eta=0.5, const double prior_rho=0.5,
							const double prior_nu=0.1, const double prior_sigma=0)
{

	const bool free_calibration = all(fixed_b != fixed_b); // NA != NA

	const int nIter = ndraws * step + from;

	const int nb = nmod.n_elem;
	const int len_g = max(bmax) + 1;
	const int nit = first.n_elem;
	const int max_cat = max(last-first)+1;
	
	// bookkeeping
	ivec cbmax(nb+1), bnit(nb, fill::zeros), cbnit(nb+1), cbmod(nb+1), citnb(itnb.n_elem+1);
	cbmax[0] = 0, cbnit[0] = 0, cbmod[0] = 0, citnb[0] = 0;
	{
		int i=0;
		for(int k=0; k<nb; k++)
		{
			for(int m=0; m<nmod[k]; m++)
				bnit[k] += mnit[i++];
			cbnit[k+1] = cbnit[k] + bnit[k];
			cbmax[k+1] = cbmax[k] + bmax[k]+1;
			cbmod[k+1] = cbmod[k] + nmod[k];
		}	
	}	
	for(int i=0; i<nit; i++)
		citnb[i+1] = itnb[i] + citnb[i];
		
	ivec m(nb);
	for (int k=0; k<nb; k++)
		m[k] = accu(scoretab.subvec(cbmax[k],cbmax[k+1]-1));

	// working variables
	vec y(max_cat), z(nb, fill::zeros);
	vec bklambda(scoretab.n_elem, fill::zeros);	
	vec pi_k(len_g, fill::zeros);
	
	// set lambda 1 for existing scores
	for(int k=0; k<nb; k++)
		for(int s=bmin[k]; s<=bmax[k]; s++)
			bklambda[cbmax[k]+s] = 1;
	
	std::vector<long double> g(len_g);
	
	vec fpwr(len_g);
	fpwr[0] = 1;
	
	//output
	mat bx(b.n_elem, ndraws);
	
	int col_index = 0;
	for (int iter=0; iter<nIter; iter++)
	{
		for (int k=0; k<nb; k++)
		{
			
			// data augmentation
			elsym(brouting[k], b, a,  bfirst.memptr() + cbnit[k], blast.memptr() + cbnit[k], bnit[k], 
				mod_min.memptr() + cbmod[k], mod_max.memptr() + cbmod[k], nmod[k],
				mnit.memptr() + cbmod[k], g);
			
			long double sm = 0;
			for (int s=0; s<=bmax[k];s++)
			{
			  if (g[s]>0){
			    pi_k[s] = R::rgamma(scoretab[cbmax[k]+s] + prior_nu, 1.0);
			    sm += pi_k[s];
			  }
			}
			
			for (int s=0; s<=bmax[k];s++)
			{
			  if(g[s]>0)
			  {
			    pi_k[s] = pi_k[s]/sm;
			    bklambda[cbmax[k]+s] = (pi_k[s]*m[k])/g[s];
			  }
			}
			z[k] = R::rgamma(m[k], 1.0/m[k]);
		}

		// bayes_items
		for(int i=0; i<nit;i++)
		{
			y.zeros();
			for(int ibnr = citnb[i]; ibnr<citnb[i+1]; ibnr++)
			{
				int k = itb[ibnr];
			  for(int j=first[i], c=0; j<=last[i]; j++,c++)//for(int s=0; s<=bmax[k]-a[last[i]]; s++)
				{
			    elsym(brouting[k], b, a,  bfirst.memptr() + cbnit[k], blast.memptr() + cbnit[k], bnit[k], 
             mod_min.memptr() + cbmod[k], mod_max.memptr() + cbmod[k], nmod[k],
             mnit.memptr() + cbmod[k], g, first[i],a[j]);
					for(int s=0; s<=bmax[k]-a[last[i]]; s++)
						y[c] += z[k] * g[s] * bklambda[cbmax[k]+s+a[j]];	
				}
			}

			for(int j=first[i], c=0; j<=last[i]; j++,c++)
				b[j] = R::rgamma(sufI[j] + prior_eta, 1/(y[c]+prior_rho));
		}

		if (free_calibration)
		{
			/*
			double f = b[0]; 
			if(a[0] != 1)
				for(int i=0; i<nit; i++)
					for(int j = first[i]; j<= last[i]; j++)
						b[j] /= std::pow(f, ((double)(a[j]))/a[0]);
			
			// Lambda
			fpwr[1] = f;
			for(int s=2; s < len_g; s++)
				fpwr[s] = fpwr[s-1] * f;
			
			for (int k=0; k<nb; k++)
			{	
				bklambda.subvec(cbmax[k], cbmax[k+1]-1) %= fpwr.head(bmax[k]+1);

				// 0 score or min score???
				if (bklambda[cbmax[k]  ] > 0)
					bklambda.subvec(cbmax[k], cbmax[k+1]-1) /= bklambda[cbmax[k]];
			}
			*/
		}
		else
		{
			b.elem(find_finite(fixed_b)) = fixed_b.elem(find_finite(fixed_b));
		}
		
		if(iter >= from && (iter-from) % step == 0)
			bx.col(col_index++) = b;
	}
	return bx;
}


// to do: possible scores (in range)
// to do: void may be slightly faster
//[[Rcpp::export]]
arma::mat  ittotmat_mst( const arma::vec& b, const arma::ivec& a, const arma::vec& c, 
             arma::ivec& first, arma::ivec& last, 
			 const int bmin, const int bmax, const int nmod, const int brouting,
			 arma::ivec& mnit, const arma::ivec& mod_min, const arma::ivec& mod_max)
{
	
	const int nI = last.n_elem;
	const int npar = accu(last-first) + nI;
	const int nscores = bmax+1;
	const vec logb = log(b);
	const vec alogc = a % log(c);
	int indx;
	  
	mat pi(npar, nscores, fill::zeros);	  
 
		std::vector<long double> g(nscores), gi(nscores);
		vec eta(npar);

		for (int s = bmin; s <= bmax; s++)
		{
			//if(ps[s] == 1)
			//{
				int k = 0; 
				eta = exp(logb + s * alogc);

				elsym(brouting, eta, a,  first.memptr(), last.memptr(), nI, 
						mod_min.memptr(), mod_max.memptr(), nmod,
						mnit.memptr(), g);

				for (int it = 0; it < nI; it++)
				{
					for (int j = first[it]; j <= last[it]; j++) 
					{
					  elsym(brouting, eta, a,  first.memptr(), last.memptr(), nI, 
             mod_min.memptr(), mod_max.memptr(), nmod,mnit.memptr(), gi, first[it],a[j]);
						indx = s-a[j];
						if ( indx >= 0 && indx < (nscores - a[last[it]])) 
							pi.at(k,s) = exp(log(eta[j]) + log(gi[indx]) - log(g[s]));
						k++;
					}
				}
			//}
		}
	return pi;
}



// might change to log scale

// for one booklet
//[[Rcpp::export]]
arma::vec elsym_C(const int routing, const arma::vec& b, const arma::ivec& a, arma::ivec& first, arma::ivec& last, 
					arma::ivec& mod_min, arma::ivec& mod_max, arma::ivec& mnit, const int max_score,
					const int item1_first=-1, const int aij=1, const int item2_first=-1, const int akl=1)
{
	std::vector<long double> g(max_score+1);
	vec out(max_score+1);
	elsym(routing, b, a, first.memptr(), last.memptr(), first.n_elem, mod_min.memptr(), mod_max.memptr(), 
       mnit.n_elem, mnit.memptr(), g, item1_first, aij, item2_first, akl);
	for(int s=0; s<=max_score; s++)
		out[s] = g[s];
	return out;	
}  

// to normalize, divide by g
//[[Rcpp::export]]
arma::vec prof_enorm(const arma::vec& b, const arma::ivec& a, arma::ivec& first, arma::ivec& last,
					const int routing, arma::ivec& mnit,
					arma::ivec& mod_min, arma::ivec& mod_max, const int max_score,
					const arma::ivec& AB) 
{

	const int nmod = mnit.n_elem;
	int idx=0, mA=0, mB=0;
	
	std::vector<long double> gA(max_score+1), gB(max_score+1);
	
	ivec cnit(nmod+1);
	cnit[0] = 0;
	for(int m=0; m<nmod; m++)
		cnit[m+1] = cnit[m] + mnit[m];	
	
	const int dimA = 1 + accu(a(conv_to<uvec>::from((1-AB) % last)));
	const int dimB = 1 + accu(a(conv_to<uvec>::from(AB % last)));
	
	cube ps(dimA, dimB, 2, fill::zeros); // contains non normalized probability of score combination, alternating slices
	
	ps(0,0,1-idx) = 1;
	
	vec out(max_score+1,fill::zeros);
	
	for(int m=0; m<nmod; m++)
	{		
		ps.slice(idx).zeros();
		std::fill(gA.begin(), gA.end(), 0);
		std::fill(gB.begin(), gB.end(), 0);		
		gA[0] = 1;
		gB[0] = 1;
		
		int MscA=0, MscB=0;

		for (int i=cnit[m]; i < cnit[m+1]; i++)
		{
			if(AB[i]==0)
			{
				for (int s=MscA; s>=0;s--)
					for (int j=last[i];j>=first[i];j--)
						gA[s+a[j]] += gA[s]*b[j];
		
				MscA += a[last[i]];
			}
			else
			{
				for (int s=MscB; s>=0; s--)
					for (int j=last[i]; j>=first[i]; j--)
						gB[s+a[j]] += gB[s]*b[j];
			
				MscB += a[last[i]];					
			}
		}
		if(routing == LAST || m==0)
		{
			MscA = std::min(MscA, mod_max[m]);
			MscB = std::min(MscB, mod_max[m]);
		}
		else
		{
			MscA = std::min(MscA, mod_max[m] - mod_max[m-1]);
			MscB = std::min(MscB, mod_max[m] - mod_max[m-1]);		
		}
		
		if(routing == LAST)
		{
			for (int sA=0; sA <= MscA; sA++)
				for (int sB=0; sB <= MscB; sB++)				
					if(sA+sB >= mod_min[m] && sA+sB <= mod_max[m])
						for(int sAp=mA;sAp>=0; sAp--)
							for(int sBp=mB;sBp>=0; sBp--)
								if(	ps(sAp, sBp,1-idx) > 0)
									ps(sAp+sA,sBp+sB,idx) += ps(sAp, sBp,1-idx) * gA[sA] * gB[sB];
									
		}
		else
		{
			for (int sA=0; sA <= MscA; sA++)
				for (int sB=0; sB <= MscB; sB++)				
					for(int sAp=mA;sAp>=0; sAp--)
						for(int sBp=mB;sBp>=0; sBp--)		
							if( ps(sAp, sBp,1-idx) > 0 && sAp + sA + sBp + sB >= mod_min[m] && sAp + sA + sBp + sB <= mod_max[m])
								ps(sAp+sA,sBp+sB,idx) += ps(sAp, sBp,1-idx) * gA[sA] * gB[sB];
				
		}

		mA += MscA;
		mB += MscB;
		idx = 1-idx;
	}	

	for(int sA=0; sA<=mA; sA++)
		for(int sB=0; sB<=mB; sB++)
			out(sA+sB) += sA * ps(sA,sB,1-idx);
		
	return out;
}

