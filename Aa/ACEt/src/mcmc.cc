#include "mcmc.h"


typedef boost::minstd_rand base_generator_type;


void ci_mh_atctet(double *result,int * num_p_mz, int * num_p_dz, 
				  int * num_col_a, int * num_col_c, int * num_col_e, 
				  double *ph_m, double *ph_d, 
				  double *B_des_a_m, double *B_des_a_d, double *B_des_c_m, double *B_des_c_d, double *B_des_e_m, double *B_des_e_d, 
				  double *G_a, double *G_c, double *G_e, double *ei_a, double *ei_c, double *ei_e, double *var_b_a, double *var_b_c, double *var_b_e, 
				  double *beta_a, double *beta_c, double *beta_e, int *D_a, int *D_c, int *D_e, 
				  int *iter_n, int *burn, double *sd_mcmc)
{

	int ITER_NUM = (*iter_n);
	int burnin = (*burn);

	base_generator_type generator(24);

	typedef std::vector<double> Vec;
	typedef std::vector<short> Vec_i;
	typedef std::vector<Vec> Mat;
	typedef std::vector<Vec_i> Mat_i;

	typedef boost::bernoulli_distribution<> distribution_type;
	typedef boost::normal_distribution<> distribution_type2;
	typedef boost::gamma_distribution<> distribution_type3;
	typedef boost::uniform_01<> distribution_type4;

	typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type;
	typedef boost::variate_generator<base_generator_type&, distribution_type2> gen_type2;
	typedef boost::variate_generator<base_generator_type&, distribution_type3> gen_type3;
	typedef boost::variate_generator<base_generator_type&, distribution_type4> gen_type4;

	// read input files

	Vec pheno_m(0);
	Vec pheno_d(0);
	int NUM_SUB_M = (*num_p_mz);
	int NUM_SUB_D = (*num_p_dz);
	int COL_A = (*num_col_a);
	int COL_C = (*num_col_c);
	int COL_E = (*num_col_e);
	Mat b_a_m;
	Mat b_a_d;
	Mat b_c_m;
	Mat b_c_d;
	Mat b_e_m;
	Mat b_e_d;
	Mat g_a;
	Mat g_c;
	Mat g_e;
	double VAR_E = (*var_b_e);
	double VAR_A = (*var_b_a);
	double VAR_C = (*var_b_c);
	Mat_i D_C;
	Mat_i D_A;
	Mat_i D_E;
	Vec e_a;
	Vec e_c;
	Vec e_e;
	Vec b_a;
	Vec b_c;
	Vec b_e;
	
	double * p = ph_m;
	for(int i = 0; i < NUM_SUB_M; i++)
	{
		double temp = *p++;
		pheno_m.push_back(temp);
	}

	p = ph_d;
	for(int i = 0; i < NUM_SUB_D; i++)
	{
		double temp = *p++;
		pheno_d.push_back(temp);
	}

	double * p2 = B_des_a_m;
	for(int i = 0; i < NUM_SUB_M; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_A; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		b_a_m.push_back(row_ge);
	}

	p2 = B_des_a_d;
	for(int i = 0; i < NUM_SUB_D; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_A; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		b_a_d.push_back(row_ge);
	}	

	p2 = B_des_c_m;
	for(int i = 0; i < NUM_SUB_M; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_C; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		b_c_m.push_back(row_ge);
	}

	p2 = B_des_c_d;
	for(int i = 0; i < NUM_SUB_D; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_C; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		b_c_d.push_back(row_ge);
	}

	p2 = B_des_e_m;
	for(int i = 0; i < NUM_SUB_M; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_E; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		b_e_m.push_back(row_ge);
	}

	p2 = B_des_e_d;
	for(int i = 0; i < NUM_SUB_D; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_E; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		b_e_d.push_back(row_ge);
	}

	p2 = G_a;
	for(int i = 0; i < COL_A; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_A; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		g_a.push_back(row_ge);
	}

	p2 = G_c;
	for(int i = 0; i < COL_C; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_C; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		g_c.push_back(row_ge);
	}

	p2 = G_e;
	for(int i = 0; i < COL_E; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_E; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		g_e.push_back(row_ge);
	}
	
	p2 = ei_a;
	double * p_b = beta_a;
	for(int i = 0; i < COL_A; i++)
	{
		double temp = *p2++;
		e_a.push_back(temp);
		temp = *p_b++;
		b_a.push_back(temp);
	}

	p2 = ei_c;
	p_b = beta_c;
	for(int i = 0; i < COL_C; i++)
	{
		double temp = *p2++;
		e_c.push_back(temp);
		temp = *p_b++;
		b_c.push_back(temp);
	}

	p2 = ei_e;
	p_b = beta_e;
	for(int i = 0; i < COL_E; i++)
	{
		double temp = *p2++;
		e_e.push_back(temp);
		temp = *p_b++;
		b_e.push_back(temp);
	}
	

	int *p3 = D_a;
	for(int i = 0; i < COL_A; i++)
	{
		Vec_i row_ge;
		for(int j = 0; j < COL_A; j++)
		{
			int temp = *p3++;
			row_ge.push_back(temp);
		}
		D_A.push_back(row_ge);
	}

	p3 = D_c;
	for(int i = 0; i < COL_C; i++)
	{
		Vec_i row_ge;
		for(int j = 0; j < COL_C; j++)
		{
			int temp = *p3++;
			row_ge.push_back(temp);
		}
		D_C.push_back(row_ge);
	}

	p3 = D_e;
	for(int i = 0; i < COL_E; i++)
	{
		Vec_i row_ge;
		for(int j = 0; j < COL_E; j++)
		{
			int temp = *p3++;
			row_ge.push_back(temp);
		}
		D_E.push_back(row_ge);
	}

	int penal_a = 2;
	if(COL_A==1)
	{penal_a = 1;}
	int penal_c = 2;
	if(COL_C==1)
	{penal_c = 1;}
	int penal_e = 2;
	if(COL_E==1)
	{penal_e = 1;}

	Vec a_t(COL_A);
	Vec c_t(COL_C);
	Vec e_t(COL_E);
	Vec tr_a_t(COL_A);
	Vec tr_c_t(COL_C);
	Vec tr_e_t(COL_E);
	//Mat mcmc_a;
	//Mat mcmc_c;
	//Mat mcmc_e;
	Vec tr_t_t(COL_A+COL_C+COL_E,0);
	Mat mcmc_t;
	for(int i = 0; i < COL_A; i++)
	{
		a_t[i] = b_a[i];
		//tr_a_t[i] = 0;
		double prod_a_t = 0;
		for(int j = 0; j < COL_A; j++)
		{
			prod_a_t += g_a[i][j]*a_t[j];
		}
		tr_a_t[i] = prod_a_t;
	}
	for(int i = 0; i < COL_C; i++)
	{
		c_t[i] = b_c[i];
		//tr_c_t[i] = 0;
		double prod_c_t = 0;
		for (int j = 0; j < COL_C; j++)
		{
			prod_c_t += g_c[i][j]*c_t[j];
		}
		tr_c_t[i] = prod_c_t;
	}
	for(int i = 0; i < COL_E; i++)
	{
		e_t[i] = b_e[i];
		//tr_e_t[i] = 0;
		double prod_e_t = 0;
		for (int j = 0; j < COL_E; j++)
		{
			prod_e_t += g_e[i][j]*e_t[j];
		}
		tr_e_t[i] = prod_e_t;
	}
	//mcmc_a.push_back(a_t);
	//mcmc_c.push_back(c_t);
	//mcmc_e.push_back(e_t);
	mcmc_t.push_back(tr_t_t);
	double lik = 0;

	double YSY_m = 0;
		double YSY_d = 0;
		double D_m = 0;
		double D_d = 0;

		for(int i = 0; i < NUM_SUB_M; )
		{
			double temp_a = 0;
			for(int j = 0; j < COL_A; j++)
			{
				//temp_a += a_t[j]*b_a_m[i][j];
				temp_a += tr_a_t[j]*b_a_m[i][j];
			}
			temp_a = exp(temp_a);
			double temp_c = 0;
			for(int j = 0; j < COL_C; j++)
			{
				//temp_c += c_t[j]*b_c_m[i][j];
				temp_c += tr_c_t[j]*b_c_m[i][j];
			}
			temp_c = exp(temp_c);
			double temp_e = 0;
			for(int j = 0; j < COL_E; j++)
			{
				//temp_e += e_t[j]*b_e_m[i][j];
				temp_e += tr_e_t[j]*b_e_m[i][j];
			}
			temp_e = exp(temp_e);

			double a11 = temp_e + temp_a + temp_c;	
			double a12 = a11 - temp_e;
			
			double a22 = a11;
			double a21 = a12;

			double denom = a11*a22-a12*a21;
			double i_a11 = a22/denom;
			double i_a22 = i_a11;
			double i_a12 = (-1)*a12/denom;
			//double i_a21 = i_a12;
			YSY_m += pheno_m[i]*pheno_m[i]*i_a11 + pheno_m[i]*pheno_m[i+1]*i_a12*2 + pheno_m[i+1]*pheno_m[i+1]*i_a22;
			D_m += log(a11*a11-a12*a12);		
			
			i += 2;
		}

		for(int i = 0; i < NUM_SUB_D; )
		{
			double temp_a = 0;
			for(int j = 0; j < COL_A; j++)
			{
				//temp_a += a_t[j]*b_a_d[i][j];
				temp_a += tr_a_t[j]*b_a_d[i][j];
			}
			temp_a = exp(temp_a);
			double temp_c = 0;
			for(int j = 0; j < COL_C; j++)
			{
				//temp_c += c_t[j]*b_c_d[i][j];
				temp_c += tr_c_t[j]*b_c_d[i][j];
			}
			temp_c = exp(temp_c);

			double temp_e = 0;
			for(int j = 0; j < COL_E; j++)
			{
				//temp_e += e_t[j]*b_e_d[i][j];
				temp_e += tr_e_t[j]*b_e_d[i][j];
			}
			temp_e = exp(temp_e);

			double a11 = temp_e + temp_a + temp_c;
			double a12 = a11 - temp_e - temp_a/2;
			
			double a22 = a11;
			double a21 = a12;

			double denom = a11*a22-a12*a21;
			double i_a11 = a22/denom;
			double i_a22 = i_a11;
			double i_a12 = (-1)*a12/denom;
			//double i_a21 = i_a12;
			YSY_d += pheno_d[i]*pheno_d[i]*i_a11 + pheno_d[i]*pheno_d[i+1]*i_a12*2 + pheno_d[i+1]*pheno_d[i+1]*i_a22;
			D_d += log(a11*a11-a12*a12);
			i += 2;
		}
		
		double temp_d_a = 0;
		if ((COL_A > 2) && (VAR_A > 0))
		{
			for (int i = 0; i < COL_A; i++)
				for (int j = 0; j < COL_A; j++)
				{
					temp_d_a += tr_a_t[i] * D_A[i][j] * tr_a_t[j];
				}
			temp_d_a /= VAR_A;
		}
		
		double temp_d_c = 0;
		if ((COL_C > 2) && (VAR_C > 0))
		{
			for (int i = 0; i < COL_C; i++)
				for (int j = 0; j < COL_C; j++)
				{
					temp_d_c += tr_c_t[i] * D_C[i][j] * tr_c_t[j];
				}
			temp_d_c /= VAR_C;
		}

		double temp_d_e = 0;
		if ((COL_E > 2) && (VAR_E > 0))
		{
			for (int i = 0; i < COL_E; i++)
				for (int j = 0; j < COL_E; j++)
				{
					temp_d_e += tr_e_t[i] * D_E[i][j] * tr_e_t[j];
				}
			temp_d_e /= VAR_E;
		}

		lik = YSY_m + YSY_d + D_m + D_d + temp_d_a + temp_d_c + temp_d_e;
		/*
		Vec map_a(COL_A);
		Vec map_c(COL_C);
		Vec map_e(COL_E);
		*/

	for(int iter = 1; iter < ITER_NUM; iter++)
	{
		double step = 0;
		Vec a_n(COL_A);
		Vec tr_a(COL_A);
		if(COL_A>2)
		{
			step = VAR_A;
			if(step>0.05)
			{step = (*sd_mcmc);}
		
			if((penal_a==2)&&(VAR_A>0))
			{
				for(int i = 0; i < (COL_A-penal_a); i++)
				{
					//gen_type2 die_gen_a(generator, distribution_type2(a_t[i],step));
					gen_type2 die_gen_a(generator, distribution_type2(a_t[i],0.1*sqrt(step/e_a[i])));
					boost::generator_iterator<gen_type2> die_a(&die_gen_a);
					a_n[i] = *die_a++;
				}
			}
			else
			{
				for(int i = 0; i < (COL_A-penal_a); i++)
				{	
					a_n[i] = 0;
				}
			}

			for(int i = (COL_A-penal_a); i < COL_A; i++)
			{
				gen_type2 die_gen_a(generator, distribution_type2(a_t[i],(*sd_mcmc)));
				boost::generator_iterator<gen_type2> die_a(&die_gen_a);
				a_n[i] = *die_a++;
			}

			for(int i = 0; i < COL_A; i++)
			{
				double prod_a_t = 0;
				for(int j = 0; j < COL_A; j++)
				{
					prod_a_t += g_a[i][j]*a_n[j];
				}
				tr_a[i] = prod_a_t;
			}
		}
		else
		{
			if(COL_A==2)
			{
				for(int i = 0; i < 2; i++)
				{
					gen_type2 die_gen_a(generator, distribution_type2(a_t[i],(*sd_mcmc)));
					boost::generator_iterator<gen_type2> die_a(&die_gen_a);
					a_n[i] = *die_a++;
					tr_a[i] = a_n[i];
				}
			}
			else
			{
				gen_type2 die_gen_a(generator, distribution_type2(a_t[0],(*sd_mcmc)));
				boost::generator_iterator<gen_type2> die_a(&die_gen_a);
				a_n[0] = *die_a++;
				tr_a[0] = a_n[0];
			}
		}

		Vec c_n(COL_C);
		Vec tr_c(COL_C);
		if(COL_C>2)
		{
			step = VAR_C;
			if(step>0.05)
			{step = (*sd_mcmc);}

			if((penal_c==2)&&(VAR_C>0))
			{
				for(int i = 0; i < (COL_C-penal_c); i++)
				{
					gen_type2 die_gen_c(generator, distribution_type2(c_t[i],0.1*sqrt(step/e_c[i])));
					boost::generator_iterator<gen_type2> die_c(&die_gen_c);
					c_n[i] = *die_c++;
				}
			}
			else
			{
				for(int i = 0; i < (COL_C-penal_c); i++)
				{	
					c_n[i] = 0;
				}
			}
			for(int i = (COL_C-penal_c); i < COL_C; i++)
			{
				gen_type2 die_gen_c(generator, distribution_type2(c_t[i],(*sd_mcmc)));
				boost::generator_iterator<gen_type2> die_c(&die_gen_c);
				c_n[i] = *die_c++;
			}

			for(int i = 0; i < COL_C; i++)
			{
				double prod_c_t = 0;
				for(int j = 0; j < COL_C; j++)
				{
					prod_c_t += g_c[i][j]*c_n[j];
				}
				tr_c[i] = prod_c_t;
			}
		}
		else
		{
			if(COL_C==2)
			{
				for(int i = 0; i < 2; i++)
				{
					gen_type2 die_gen_c(generator, distribution_type2(c_t[i],(*sd_mcmc)));
					boost::generator_iterator<gen_type2> die_c(&die_gen_c);
					c_n[i] = *die_c++;
					tr_c[i] = c_n[i];
				}
			}
			else
			{
				gen_type2 die_gen_c(generator, distribution_type2(c_t[0],(*sd_mcmc)));
				boost::generator_iterator<gen_type2> die_c(&die_gen_c);
				c_n[0] = *die_c++;
				tr_c[0] = c_n[0];
			}
		}

		Vec e_n(COL_E);
		Vec tr_e(COL_E);
		
		if(COL_E>2)
		{
			step = VAR_E;
			if(step>0.05)
			{step = (*sd_mcmc);}

		
			if((penal_e==2)&&(VAR_E>0))
			{
				for(int i = 0; i < (COL_E-penal_e); i++)
				{
					gen_type2 die_gen_e(generator, distribution_type2(e_t[i],0.1*sqrt(step/e_e[i])));
					boost::generator_iterator<gen_type2> die_e(&die_gen_e);
					e_n[i] = *die_e++;
				}
			}
			else
			{
				for(int i = 0; i < (COL_E-penal_e); i++)
				{	
					e_n[i] = 0;
				}
			}
			for(int i = (COL_E-penal_e); i < COL_E; i++)
			{
				gen_type2 die_gen_e(generator, distribution_type2(e_t[i],(*sd_mcmc)));
				boost::generator_iterator<gen_type2> die_e(&die_gen_e);
				e_n[i] = *die_e++;
			}

			for(int i = 0; i < COL_E; i++)
			{
				double prod_e_t = 0;
				for(int j = 0; j < COL_E; j++)
				{
					prod_e_t += g_e[i][j]*e_n[j];
				}
				tr_e[i] = prod_e_t;
			}
		}
		else
		{
			if(COL_E==2)
			{
				for(int i = 0; i < 2; i++)
				{
					gen_type2 die_gen_e(generator, distribution_type2(e_t[i],(*sd_mcmc)));
					boost::generator_iterator<gen_type2> die_e(&die_gen_e);
					e_n[i] = *die_e++;
					tr_e[i] = e_n[i];
				}
			}
			else
			{
				gen_type2 die_gen_e(generator, distribution_type2(e_t[0],(*sd_mcmc)));
				boost::generator_iterator<gen_type2> die_e(&die_gen_e);
				e_n[0] = *die_e++;
				tr_e[0] = e_n[0];
			}
		}
		
		double new_lik = 0;
		YSY_m = 0;
		YSY_d = 0;
		D_m = 0;
		D_d = 0;

		for(int i = 0; i < NUM_SUB_M; )
		{
			double temp_a = 0;
			for(int j = 0; j < COL_A; j++)
			{
				temp_a += tr_a[j]*b_a_m[i][j];
			}
			temp_a = exp(temp_a);
			double temp_c = 0;
			for(int j = 0; j < COL_C; j++)
			{
				temp_c += tr_c[j]*b_c_m[i][j];
			}
			temp_c = exp(temp_c);
			double temp_e = 0;
			for(int j = 0; j < COL_E; j++)
			{
				temp_e += tr_e[j]*b_e_m[i][j];
			}
			temp_e = exp(temp_e);

			double a11 = temp_e + temp_a + temp_c;	
			double a12 = a11 - temp_e;
			
			double a22 = a11;
			double a21 = a12;

			double denom = a11*a22-a12*a21;
			double i_a11 = a22/denom;
			double i_a22 = i_a11;
			double i_a12 = (-1)*a12/denom;
			//double i_a21 = i_a12;
			YSY_m += pheno_m[i]*pheno_m[i]*i_a11 + pheno_m[i]*pheno_m[i+1]*i_a12*2 + pheno_m[i+1]*pheno_m[i+1]*i_a22;
			D_m += log(a11*a11-a12*a12);
			i += 2;
		}

		for(int i = 0; i < NUM_SUB_D; )
		{
			double temp_a = 0;
			for(int j = 0; j < COL_A; j++)
			{
				temp_a += tr_a[j]*b_a_d[i][j];
			}
			temp_a = exp(temp_a);
			double temp_c = 0;
			for(int j = 0; j < COL_C; j++)
			{
				temp_c += tr_c[j]*b_c_d[i][j];
			}
			temp_c = exp(temp_c);
			double temp_e = 0;
			for(int j = 0; j < COL_E; j++)
			{
				temp_e += tr_e[j]*b_e_d[i][j];
			}
			temp_e = exp(temp_e);

			double a11 = temp_e + temp_a + temp_c;
			double a12 = a11 - temp_e - temp_a/2;
			
			double a22 = a11;
			double a21 = a12;

			double denom = a11*a22-a12*a21;
			double i_a11 = a22/denom;
			double i_a22 = i_a11;
			double i_a12 = (-1)*a12/denom;
			//double i_a21 = i_a12;
			YSY_d += pheno_d[i]*pheno_d[i]*i_a11 + pheno_d[i]*pheno_d[i+1]*i_a12*2 + pheno_d[i+1]*pheno_d[i+1]*i_a22;
			D_d += log(a11*a11-a12*a12);
			i += 2;
		}

		temp_d_a = 0;
		if((COL_A>2)&&(VAR_A>0))
		{
			for(int i = 0; i < COL_A; i++)
				for(int j = 0; j < COL_A; j++)
				{
					temp_d_a += tr_a[i]*D_A[i][j]*tr_a[j];
				}
			temp_d_a /= VAR_A;
		}

		temp_d_c = 0;
		if((COL_C>2)&&(VAR_C>0))
		{
			for(int i = 0; i < COL_C; i++)
				for(int j = 0; j < COL_C; j++)
				{
					temp_d_c += tr_c[i]*D_C[i][j]*tr_c[j];
				}
			temp_d_c /= VAR_C;
		}

		temp_d_e = 0;
		if((COL_E>2)&&(VAR_E>0))
		{
			for(int i = 0; i < COL_E; i++)
				for(int j = 0; j < COL_E; j++)
				{
					temp_d_e += tr_e[i]*D_E[i][j]*tr_e[j];
				}
			temp_d_e /= VAR_E;
		}

		new_lik = YSY_m + YSY_d + D_m + D_d + temp_d_a + temp_d_c + temp_d_e;

		gen_type4 die_gen_u1(generator, distribution_type4());
		boost::generator_iterator<gen_type4> die_u1(&die_gen_u1);
		double u = *die_u1++;
		double r = exp((-0.5)*(new_lik-lik));
		/*if(new_lik<lik)
		{
			map_a = tr_a;
			map_c = tr_c;
			map_e = tr_e;
		}*/
		if(u<r)
		{
			a_t = a_n;
			c_t = c_n;
			e_t = e_n;
			tr_a_t = tr_a;
			tr_c_t = tr_c;
			tr_e_t = tr_e;
			lik = new_lik;
		}
		//mcmc_a.push_back(tr_a_t);
		//mcmc_c.push_back(tr_c_t);
		//mcmc_e.push_back(tr_e_t);
		tr_t_t = tr_a_t;
		tr_t_t.insert(tr_t_t.end(),tr_c_t.begin(),tr_c_t.end());
		tr_t_t.insert(tr_t_t.end(),tr_e_t.begin(),tr_e_t.end());
		mcmc_t.push_back(tr_t_t);
		// std::cout<<iter;
	}
	
	Vec post_mean_a(COL_A);
	Mat post_cov_a;
	Vec post_mean_c(COL_C);
	Mat post_cov_c;
	Vec post_mean_e(COL_E);
	Mat post_cov_e;
	for(int i = 0; i < COL_A; i++){post_mean_a[i] = 0;}
	for(int i = 0; i < COL_C; i++){post_mean_c[i] = 0;}
	for(int i = 0; i < COL_E; i++){post_mean_e[i] = 0;}

	for(int i = burnin; i < ITER_NUM; i++)
	{
		for(int j = 0; j < COL_A; j++)
		{
			//post_mean_a[j] += mcmc_a[i][j];
			post_mean_a[j] += mcmc_t[i][j];
		}
		for(int j = 0; j < COL_C; j++)
		{
			//post_mean_c[j] += mcmc_c[i][j];
			post_mean_c[j] += mcmc_t[i][j+COL_A];
		}
		for(int j = 0; j < COL_E; j++)
		{
			//post_mean_e[j] += mcmc_e[i][j];
			post_mean_e[j] += mcmc_t[i][j+COL_A+COL_C];
		}
	}
	int iter_count = ITER_NUM - burnin;
	for(int i = 0; i < COL_A; i++)
	{
		post_mean_a[i] /= iter_count;
		//if((penal_a==2)&&(VAR_A>0))
		//{
			result[i] = post_mean_a[i];
		/*}
		else
		{
			result[i] = map_a[i];
		}*/
	}
	
	for(int i = 0; i < COL_C; i++)
	{
		post_mean_c[i] /= iter_count;
		//if((penal_c==2)&&(VAR_C>0))
		//{
			result[i+COL_A] = post_mean_c[i];
		/*}
		else
		{
			result[i+COL_A] = map_c[i];
		}*/
	}

	for(int i = 0; i < COL_E; i++)
	{
		post_mean_e[i] /= iter_count;
		//if((penal_e==2)&&(VAR_E>0))
		//{
			result[i+COL_A+COL_C] = post_mean_e[i];
		/*}
		else
		{
			result[i+COL_A+COL_C] = map_e[i];
		}*/
	}
	
	int COL_T = COL_A + COL_C + COL_E;
	
	int cov_i = 0;
	for(int i = 0; i < COL_T; i++)
	{
		for(int j = i; j < COL_T; j++)
		{
			double temp_cov = 0;
			for(int k = burnin; k < ITER_NUM; k++)
			{
				// temp_cov += (mcmc_a[k][i] - post_mean_a[i])*(mcmc_a[k][j] - post_mean_a[j]);
				temp_cov += (mcmc_t[k][i] - result[i])*(mcmc_t[k][j] - result[j]);
			}
			temp_cov /= iter_count;
			//post_cov_a[i][j] = temp_cov;
			//post_cov_a[j][i] = temp_cov;
			result[cov_i+COL_A+COL_C+COL_E] = temp_cov;
			cov_i++;
		}
	}
	/*
	int cov_j = 0;
	for(int i = 0; i < COL_C; i++)
	{
		for(int j = i; j < COL_C; j++)
		{
			double temp_cov = 0;
			for(int k = burnin; k < ITER_NUM; k++)
			{
				temp_cov += (mcmc_c[k][i] - post_mean_c[i])*(mcmc_c[k][j] - post_mean_c[j]);
			}
			temp_cov /= iter_count;
			//post_cov_c[i][j] = temp_cov;
			//post_cov_c[j][i] = temp_cov;
			result[cov_j+cov_i+COL_A+COL_C+COL_E] = temp_cov;
			cov_j++;
		}
	}

	int cov_k = 0;
	for(int i = 0; i < COL_E; i++)
	{
		for(int j = i; j < COL_E; j++)
		{
			double temp_cov = 0;
			for(int k = burnin; k < ITER_NUM; k++)
			{
				temp_cov += (mcmc_e[k][i] - post_mean_e[i])*(mcmc_e[k][j] - post_mean_e[j]);
			}
			temp_cov /= iter_count;
			//post_cov_c[i][j] = temp_cov;
			//post_cov_c[j][i] = temp_cov;
			result[cov_k+cov_j+cov_i+COL_A+COL_C+COL_E] = temp_cov;
			cov_k++;
		}
	}
	*/
}

void ci_mh_atctet_2(double *result,int * num_p_mz, int * num_p_dz, 
  			  int * num_col_a, int * num_col_c, int * num_col_e, 
				  double *ph_m, double *ph_d, 
				  double *B_des_a_m, double *B_des_a_d, double *B_des_c_m, double *B_des_c_d, double *B_des_e_m, double *B_des_e_d, 
				  double *var_b_a, double *var_b_c, double *var_b_e, 
				  double *D_a, double *D_c, double *D_e, 
				  int *iter_n, int *burn, double *sd_mcmc)
{

	int ITER_NUM = (*iter_n);
	int burnin = (*burn);

	base_generator_type generator(24);

	typedef std::vector<double> Vec;
	typedef std::vector<short> Vec_i;
	typedef std::vector<Vec> Mat;
	typedef std::vector<Vec_i> Mat_i;

	typedef boost::bernoulli_distribution<> distribution_type;
	typedef boost::normal_distribution<> distribution_type2;
	typedef boost::gamma_distribution<> distribution_type3;
	typedef boost::uniform_01<> distribution_type4;

	typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type;
	typedef boost::variate_generator<base_generator_type&, distribution_type2> gen_type2;
	typedef boost::variate_generator<base_generator_type&, distribution_type3> gen_type3;
	typedef boost::variate_generator<base_generator_type&, distribution_type4> gen_type4;

	// read input files

	Vec pheno_m(0);
	Vec pheno_d(0);
	int NUM_SUB_M = (*num_p_mz);
	int NUM_SUB_D = (*num_p_dz);
	int COL_A = (*num_col_a);
	int COL_C = (*num_col_c);
	int COL_E = (*num_col_e);
	Mat b_a_m;
	Mat b_a_d;
	Mat b_c_m;
	Mat b_c_d;
	Mat b_e_m;
	Mat b_e_d;
	double VAR_E = (*var_b_e);
	double VAR_A = (*var_b_a);
	double VAR_C = (*var_b_c);
	Mat D_C;
	Mat D_A;
	Mat D_E;
	
	double * p = ph_m;
	for(int i = 0; i < NUM_SUB_M; i++)
	{
		double temp = *p++;
		pheno_m.push_back(temp);
	}

	p = ph_d;
	for(int i = 0; i < NUM_SUB_D; i++)
	{
		double temp = *p++;
		pheno_d.push_back(temp);
	}

	double * p2 = B_des_a_m;
	for(int i = 0; i < NUM_SUB_M; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_A; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		b_a_m.push_back(row_ge);
	}

	p2 = B_des_a_d;
	for(int i = 0; i < NUM_SUB_D; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_A; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		b_a_d.push_back(row_ge);
	}	

	p2 = B_des_c_m;
	for(int i = 0; i < NUM_SUB_M; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_C; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		b_c_m.push_back(row_ge);
	}

	p2 = B_des_c_d;
	for(int i = 0; i < NUM_SUB_D; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_C; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		b_c_d.push_back(row_ge);
	}

	p2 = B_des_e_m;
	for(int i = 0; i < NUM_SUB_M; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_E; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		b_e_m.push_back(row_ge);
	}

	p2 = B_des_e_d;
	for(int i = 0; i < NUM_SUB_D; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_E; j++)
		{
			double temp = *p2++;
			row_ge.push_back(temp);
		}
		b_e_d.push_back(row_ge);
	}

	double *p3 = D_a;
	for(int i = 0; i < COL_A; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_A; j++)
		{
			double temp = *p3++;
			row_ge.push_back(temp);
		}
		D_A.push_back(row_ge);
	}

	p3 = D_c;
	for(int i = 0; i < COL_C; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_C; j++)
		{
			double temp = *p3++;
			row_ge.push_back(temp);
		}
		D_C.push_back(row_ge);
	}

	p3 = D_e;
	for(int i = 0; i < COL_E; i++)
	{
		Vec row_ge;
		for(int j = 0; j < COL_E; j++)
		{
			double temp = *p3++;
			row_ge.push_back(temp);
		}
		D_E.push_back(row_ge);
	}

	int penal_a = 2;
	if(COL_A==1)
	{penal_a = 1;}
	int penal_c = 2;
	if(COL_C==1)
	{penal_c = 1;}
	int penal_e = 2;
	if(COL_E==1)
	{penal_e = 1;}

	Vec a_t(COL_A);
	Vec c_t(COL_C);
	Vec e_t(COL_E);
	//Mat mcmc_a;
	//Mat mcmc_c;
	//Mat mcmc_e;
	Vec t_t(COL_A+COL_C+COL_E,0);
	Mat mcmc_t;
	for(int i = 0; i < COL_A; i++)
	{
		a_t[i] = 0;
	}
	for(int i = 0; i < COL_C; i++)
	{
		c_t[i] = 0;
	}
	for(int i = 0; i < COL_E; i++)
	{
		e_t[i] = 0;
	}
	//mcmc_a.push_back(a_t);
	//mcmc_c.push_back(c_t);
	//mcmc_e.push_back(e_t);
	mcmc_t.push_back(t_t);
	double lik = 0;

	double YSY_m = 0;
		double YSY_d = 0;
		double D_m = 0;
		double D_d = 0;

		for(int i = 0; i < NUM_SUB_M; )
		{
			double temp_a = 0;
			for(int j = 0; j < COL_A; j++)
			{
				temp_a += a_t[j]*b_a_m[i][j];
			}
			temp_a = exp(temp_a);
			double temp_c = 0;
			for(int j = 0; j < COL_C; j++)
			{
				temp_c += c_t[j]*b_c_m[i][j];
			}
			temp_c = exp(temp_c);
			double temp_e = 0;
			for(int j = 0; j < COL_E; j++)
			{
				temp_e += e_t[j]*b_e_m[i][j];
			}
			temp_e = exp(temp_e);

			double a11 = temp_e + temp_a + temp_c;	
			double a12 = a11 - temp_e;
			
			double a22 = a11;
			double a21 = a12;

			double denom = a11*a22-a12*a21;
			double i_a11 = a22/denom;
			double i_a22 = i_a11;
			double i_a12 = (-1)*a12/denom;
			//double i_a21 = i_a12;
			YSY_m += pheno_m[i]*pheno_m[i]*i_a11 + pheno_m[i]*pheno_m[i+1]*i_a12*2 + pheno_m[i+1]*pheno_m[i+1]*i_a22;
			D_m += log(a11*a11-a12*a12);		
			
			i += 2;
		}

		for(int i = 0; i < NUM_SUB_D; )
		{
			double temp_a = 0;
			for(int j = 0; j < COL_A; j++)
			{
				temp_a += a_t[j]*b_a_d[i][j];
			}
			temp_a = exp(temp_a);
			double temp_c = 0;
			for(int j = 0; j < COL_C; j++)
			{
				temp_c += c_t[j]*b_c_d[i][j];
			}
			temp_c = exp(temp_c);

			double temp_e = 0;
			for(int j = 0; j < COL_E; j++)
			{
				temp_e += e_t[j]*b_e_d[i][j];
			}
			temp_e = exp(temp_e);

			double a11 = temp_e + temp_a + temp_c;
			double a12 = a11 - temp_e - temp_a/2;
			
			double a22 = a11;
			double a21 = a12;

			double denom = a11*a22-a12*a21;
			double i_a11 = a22/denom;
			double i_a22 = i_a11;
			double i_a12 = (-1)*a12/denom;
			//double i_a21 = i_a12;
			YSY_d += pheno_d[i]*pheno_d[i]*i_a11 + pheno_d[i]*pheno_d[i+1]*i_a12*2 + pheno_d[i+1]*pheno_d[i+1]*i_a22;
			D_d += log(a11*a11-a12*a12);
			i += 2;
		}
		
		double temp_d_a = 0;
		
		double temp_d_c = 0;

		double temp_d_e = 0;

		lik = YSY_m + YSY_d + D_m + D_d + temp_d_a + temp_d_c + temp_d_e;
		Vec map_a(COL_A);
		Vec map_c(COL_C);
		Vec map_e(COL_E);

	for(int iter = 1; iter < ITER_NUM; iter++)
	{
		double step = 0;
		Vec a_n(COL_A);
		if(COL_A>2)
		{
			step = VAR_A;
			if(step>0.05)
			{step = (*sd_mcmc);}
		
			if((penal_a==2)&&(VAR_A>0))
			{
				for(int i = 0; i < (COL_A-penal_a); i++)
				{
					gen_type2 die_gen_a(generator, distribution_type2(a_t[i],step));
					boost::generator_iterator<gen_type2> die_a(&die_gen_a);
					a_n[i] = *die_a++;
				}
			}
			else
			{
				for(int i = 0; i < (COL_A-penal_a); i++)
				{	
					a_n[i] = 0;
				}
			}

			for(int i = (COL_A-penal_a); i < COL_A; i++)
			{
				gen_type2 die_gen_a(generator, distribution_type2(a_t[i],(*sd_mcmc)));
				boost::generator_iterator<gen_type2> die_a(&die_gen_a);
				a_n[i] = *die_a++;
			}

		}
		else
		{
			if(COL_A==2)
			{
				for(int i = 0; i < 2; i++)
				{
					gen_type2 die_gen_a(generator, distribution_type2(a_t[i],(*sd_mcmc)));
					boost::generator_iterator<gen_type2> die_a(&die_gen_a);
					a_n[i] = *die_a++;
				}
			}
			else
			{
				gen_type2 die_gen_a(generator, distribution_type2(a_t[0],(*sd_mcmc)));
				boost::generator_iterator<gen_type2> die_a(&die_gen_a);
				a_n[0] = *die_a++;
			}
		}

		Vec c_n(COL_C);
		if(COL_C>2)
		{
			step = VAR_C;
			if(step>0.05)
			{step = (*sd_mcmc);}

			if((penal_c==2)&&(VAR_C>0))
			{
				for(int i = 0; i < (COL_C-penal_c); i++)
				{
					gen_type2 die_gen_c(generator, distribution_type2(c_t[i],step));
					boost::generator_iterator<gen_type2> die_c(&die_gen_c);
					c_n[i] = *die_c++;
				}
			}
			else
			{
				for(int i = 0; i < (COL_C-penal_c); i++)
				{	
					c_n[i] = 0;
				}
			}
			for(int i = (COL_C-penal_c); i < COL_C; i++)
			{
				gen_type2 die_gen_c(generator, distribution_type2(c_t[i],(*sd_mcmc)));
				boost::generator_iterator<gen_type2> die_c(&die_gen_c);
				c_n[i] = *die_c++;
			}

		}
		else
		{
			if(COL_C==2)
			{
				for(int i = 0; i < 2; i++)
				{
					gen_type2 die_gen_c(generator, distribution_type2(c_t[i],(*sd_mcmc)));
					boost::generator_iterator<gen_type2> die_c(&die_gen_c);
					c_n[i] = *die_c++;
				}
			}
			else
			{
				gen_type2 die_gen_c(generator, distribution_type2(c_t[0],(*sd_mcmc)));
				boost::generator_iterator<gen_type2> die_c(&die_gen_c);
				c_n[0] = *die_c++;
			}
		}

		Vec e_n(COL_E);
		
		if(COL_E>2)
		{
			step = VAR_E;
			if(step>0.05)
			{step = (*sd_mcmc);}

		
			if((penal_e==2)&&(VAR_E>0))
			{
				for(int i = 0; i < (COL_E-penal_e); i++)
				{
					gen_type2 die_gen_e(generator, distribution_type2(e_t[i],step));
					boost::generator_iterator<gen_type2> die_e(&die_gen_e);
					e_n[i] = *die_e++;
				}
			}
			else
			{
				for(int i = 0; i < (COL_E-penal_e); i++)
				{	
					e_n[i] = 0;
				}
			}
			for(int i = (COL_E-penal_e); i < COL_E; i++)
			{
				gen_type2 die_gen_e(generator, distribution_type2(e_t[i],(*sd_mcmc)));
				boost::generator_iterator<gen_type2> die_e(&die_gen_e);
				e_n[i] = *die_e++;
			}

		}
		else
		{
			if(COL_E==2)
			{
				for(int i = 0; i < 2; i++)
				{
					gen_type2 die_gen_e(generator, distribution_type2(e_t[i],(*sd_mcmc)));
					boost::generator_iterator<gen_type2> die_e(&die_gen_e);
					e_n[i] = *die_e++;
				}
			}
			else
			{
				gen_type2 die_gen_e(generator, distribution_type2(e_t[0],(*sd_mcmc)));
				boost::generator_iterator<gen_type2> die_e(&die_gen_e);
				e_n[0] = *die_e++;
			}
		}
		
		double new_lik = 0;
		YSY_m = 0;
		YSY_d = 0;
		D_m = 0;
		D_d = 0;

		for(int i = 0; i < NUM_SUB_M; )
		{
			double temp_a = 0;
			for(int j = 0; j < COL_A; j++)
			{
				temp_a += a_n[j]*b_a_m[i][j];
			}
			temp_a = exp(temp_a);
			double temp_c = 0;
			for(int j = 0; j < COL_C; j++)
			{
				temp_c += c_n[j]*b_c_m[i][j];
			}
			temp_c = exp(temp_c);
			double temp_e = 0;
			for(int j = 0; j < COL_E; j++)
			{
				temp_e += e_n[j]*b_e_m[i][j];
			}
			temp_e = exp(temp_e);

			double a11 = temp_e + temp_a + temp_c;	
			double a12 = a11 - temp_e;
			
			double a22 = a11;
			double a21 = a12;

			double denom = a11*a22-a12*a21;
			double i_a11 = a22/denom;
			double i_a22 = i_a11;
			double i_a12 = (-1)*a12/denom;
			//double i_a21 = i_a12;
			YSY_m += pheno_m[i]*pheno_m[i]*i_a11 + pheno_m[i]*pheno_m[i+1]*i_a12*2 + pheno_m[i+1]*pheno_m[i+1]*i_a22;
			D_m += log(a11*a11-a12*a12);
			i += 2;
		}

		for(int i = 0; i < NUM_SUB_D; )
		{
			double temp_a = 0;
			for(int j = 0; j < COL_A; j++)
			{
				temp_a += a_n[j]*b_a_d[i][j];
			}
			temp_a = exp(temp_a);
			double temp_c = 0;
			for(int j = 0; j < COL_C; j++)
			{
				temp_c += c_n[j]*b_c_d[i][j];
			}
			temp_c = exp(temp_c);
			double temp_e = 0;
			for(int j = 0; j < COL_E; j++)
			{
				temp_e += e_n[j]*b_e_d[i][j];
			}
			temp_e = exp(temp_e);

			double a11 = temp_e + temp_a + temp_c;
			double a12 = a11 - temp_e - temp_a/2;
			
			double a22 = a11;
			double a21 = a12;

			double denom = a11*a22-a12*a21;
			double i_a11 = a22/denom;
			double i_a22 = i_a11;
			double i_a12 = (-1)*a12/denom;
			//double i_a21 = i_a12;
			YSY_d += pheno_d[i]*pheno_d[i]*i_a11 + pheno_d[i]*pheno_d[i+1]*i_a12*2 + pheno_d[i+1]*pheno_d[i+1]*i_a22;
			D_d += log(a11*a11-a12*a12);
			i += 2;
		}

		temp_d_a = 0;
		if((COL_A>2)&&(VAR_A>0))
		{
			for(int i = 0; i < COL_A; i++)
				for(int j = 0; j < COL_A; j++)
				{
					temp_d_a += a_n[i]*D_A[i][j]*a_n[j];
				}
			temp_d_a /= VAR_A;
		}

		temp_d_c = 0;
		if((COL_C>2)&&(VAR_C>0))
		{
			for(int i = 0; i < COL_C; i++)
				for(int j = 0; j < COL_C; j++)
				{
					temp_d_c += c_n[i]*D_C[i][j]*c_n[j];
				}
			temp_d_c /= VAR_C;
		}

		temp_d_e = 0;
		if((COL_E>2)&&(VAR_E>0))
		{
			for(int i = 0; i < COL_E; i++)
				for(int j = 0; j < COL_E; j++)
				{
					temp_d_e += e_n[i]*D_E[i][j]*e_n[j];
				}
			temp_d_e /= VAR_E;
		}

		new_lik = YSY_m + YSY_d + D_m + D_d + temp_d_a + temp_d_c + temp_d_e;

		gen_type4 die_gen_u1(generator, distribution_type4());
		boost::generator_iterator<gen_type4> die_u1(&die_gen_u1);
		double u = *die_u1++;
		double r = exp((-0.5)*(new_lik-lik));
		if(u<r)
		{
			a_t = a_n;
			c_t = c_n;
			e_t = e_n;
			lik = new_lik;
		}
		//mcmc_a.push_back(tr_a_t);
		//mcmc_c.push_back(tr_c_t);
		//mcmc_e.push_back(tr_e_t);
		t_t = a_t;
		t_t.insert(t_t.end(),c_t.begin(),c_t.end());
		t_t.insert(t_t.end(),e_t.begin(),e_t.end());
		mcmc_t.push_back(t_t);
		// std::cout<<iter;
	}
	
	Vec post_mean_a(COL_A);
	Mat post_cov_a;
	Vec post_mean_c(COL_C);
	Mat post_cov_c;
	Vec post_mean_e(COL_E);
	Mat post_cov_e;
	for(int i = 0; i < COL_A; i++){post_mean_a[i] = 0;}
	for(int i = 0; i < COL_C; i++){post_mean_c[i] = 0;}
	for(int i = 0; i < COL_E; i++){post_mean_e[i] = 0;}

	for(int i = burnin; i < ITER_NUM; i++)
	{
		for(int j = 0; j < COL_A; j++)
		{
			//post_mean_a[j] += mcmc_a[i][j];
			post_mean_a[j] += mcmc_t[i][j];
		}
		for(int j = 0; j < COL_C; j++)
		{
			//post_mean_c[j] += mcmc_c[i][j];
			post_mean_c[j] += mcmc_t[i][j+COL_A];
		}
		for(int j = 0; j < COL_E; j++)
		{
			//post_mean_e[j] += mcmc_e[i][j];
			post_mean_e[j] += mcmc_t[i][j+COL_A+COL_C];
		}
	}
	int iter_count = ITER_NUM - burnin;
	for(int i = 0; i < COL_A; i++)
	{
		post_mean_a[i] /= iter_count;
		//if((penal_a==2)&&(VAR_A>0))
		//{
			result[i] = post_mean_a[i];
		/*}
		else
		{
			result[i] = map_a[i];
		}*/
	}
	
	for(int i = 0; i < COL_C; i++)
	{
		post_mean_c[i] /= iter_count;
		//if((penal_c==2)&&(VAR_C>0))
		//{
			result[i+COL_A] = post_mean_c[i];
		/*}
		else
		{
			result[i+COL_A] = map_c[i];
		}*/
	}

	for(int i = 0; i < COL_E; i++)
	{
		post_mean_e[i] /= iter_count;
		//if((penal_e==2)&&(VAR_E>0))
		//{
			result[i+COL_A+COL_C] = post_mean_e[i];
		/*}
		else
		{
			result[i+COL_A+COL_C] = map_e[i];
		}*/
	}
	
	int COL_T = COL_A + COL_C + COL_E;
	
	int cov_i = 0;
	for(int i = 0; i < COL_T; i++)
	{
		for(int j = i; j < COL_T; j++)
		{
			double temp_cov = 0;
			for(int k = burnin; k < ITER_NUM; k++)
			{
				// temp_cov += (mcmc_a[k][i] - post_mean_a[i])*(mcmc_a[k][j] - post_mean_a[j]);
				temp_cov += (mcmc_t[k][i] - result[i])*(mcmc_t[k][j] - result[j]);
			}
			temp_cov /= iter_count;
			//post_cov_a[i][j] = temp_cov;
			//post_cov_a[j][i] = temp_cov;
			result[cov_i+COL_A+COL_C+COL_E] = temp_cov;
			cov_i++;
		}
	}
	/*
	int cov_j = 0;
	for(int i = 0; i < COL_C; i++)
	{
		for(int j = i; j < COL_C; j++)
		{
			double temp_cov = 0;
			for(int k = burnin; k < ITER_NUM; k++)
			{
				temp_cov += (mcmc_c[k][i] - post_mean_c[i])*(mcmc_c[k][j] - post_mean_c[j]);
			}
			temp_cov /= iter_count;
			//post_cov_c[i][j] = temp_cov;
			//post_cov_c[j][i] = temp_cov;
			result[cov_j+cov_i+COL_A+COL_C+COL_E] = temp_cov;
			cov_j++;
		}
	}

	int cov_k = 0;
	for(int i = 0; i < COL_E; i++)
	{
		for(int j = i; j < COL_E; j++)
		{
			double temp_cov = 0;
			for(int k = burnin; k < ITER_NUM; k++)
			{
				temp_cov += (mcmc_e[k][i] - post_mean_e[i])*(mcmc_e[k][j] - post_mean_e[j]);
			}
			temp_cov /= iter_count;
			//post_cov_c[i][j] = temp_cov;
			//post_cov_c[j][i] = temp_cov;
			result[cov_k+cov_j+cov_i+COL_A+COL_C+COL_E] = temp_cov;
			cov_k++;
		}
	}
	*/
}
