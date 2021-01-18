
/* This model use a approximation in [-T,T] of the liklihood
 * Let g: x-> 1/2*x + log(1-1/(1+exp(-x)))
 * Let h: x-> g(T*x)
 * the approximation is made for h in [-1,1], definied by the order 2*K
 * polynomial P minimizing the norm of the difference with h in [-1,1] and
 * respecting the constraint P''<max(-f'')*tol, the have quasi concavity.
 * tol is choose to be 1e-3, and K=7 (degree 14 of the polynomial)
 * Coefs (only even), starting for the constant:
 *  -6.9972974002661037e-01
 *  -1.1648509443490106e+01
 *  +3.0762276187045025e+01
 *  -7.8199766154064548e+01
 *  +1.3263488729186136e+02
 *  -1.3718047517968654e+02
 *  +7.7835307129731049e+01
 *  -1.8514863011774406e+01
 *
 *
 * bcf as defined as follow:
 * bcf[k][t] = coefs(2*k)*bincoef(2*k,t)
 * and therefore
 * P(x+y) = sum_{k=0..K} sum_{t=0..2*k} bcf[k][t] * x^t y^(2*k-t)
 * and therefore g is approximated by:
 * g2(x+y) = sum_{k=0..K} sum_{t=0..2*k} bcf[k][t] * (x/T)^t (y/T)^(2*k-t)
 */

#define INIT_bcf_000 { -7.0296697793969232e-01}
#define INIT_bcf_001 { -2.5747451970441507e+01, -5.1494903940883013e+01, -2.5747451970441507e+01}
#define INIT_bcf_002 { +1.4653501920584080e+02, +5.8614007682336319e+02, +8.7921011523504478e+02, +5.8614007682336319e+02, +1.4653501920584080e+02}
#define INIT_bcf_003 { -8.4109650268336179e+02, -5.0465790161001705e+03, -1.2616447540250427e+04, -1.6821930053667234e+04, -1.2616447540250427e+04, -5.0465790161001705e+03, -8.4109650268336179e+02}
#define INIT_bcf_004 { +3.5293685520550280e+03, +2.8234948416440224e+04, +9.8822319457540783e+04, +1.9764463891508157e+05, +2.4705579864385194e+05, +1.9764463891508157e+05, +9.8822319457540783e+04, +2.8234948416440224e+04, +3.5293685520550280e+03}
#define INIT_bcf_005 { -1.0415093220417288e+04, -1.0415093220417287e+05, -4.6867919491877797e+05, -1.2498111864500744e+06, -2.1871695762876305e+06, -2.6246034915451566e+06, -2.1871695762876305e+06, -1.2498111864500744e+06, -4.6867919491877797e+05, -1.0415093220417287e+05, -1.0415093220417288e+04}
#define INIT_bcf_006 { +2.1440933275506373e+04, +2.5729119930607648e+05, +1.4151015961834206e+06, +4.7170053206114024e+06, +1.0613261971375655e+07, +1.6981219154201049e+07, +1.9811422346567888e+07, +1.6981219154201049e+07, +1.0613261971375655e+07, +4.7170053206114024e+06, +1.4151015961834206e+06, +2.5729119930607648e+05, +2.1440933275506373e+04}
#define INIT_bcf_007 { -3.0502522204942401e+04, -4.2703531086919364e+05, -2.7757295206497586e+06, -1.1102918082599035e+07, -3.0533024727147344e+07, -6.1066049454294689e+07, -9.1599074181442037e+07, -1.0468465620736232e+08, -9.1599074181442037e+07, -6.1066049454294689e+07, -3.0533024727147344e+07, -1.1102918082599035e+07, -2.7757295206497586e+06, -4.2703531086919364e+05, -3.0502522204942401e+04}
#define INIT_bcf_008 { +2.9322607715088521e+04, +4.6916172344141634e+05, +3.5187129258106225e+06, +1.6420660320449572e+07, +5.3367146041461110e+07, +1.2808115049950667e+08, +2.3481544258242887e+08, +3.3545063226061267e+08, +3.7738196129318929e+08, +3.3545063226061267e+08, +2.3481544258242887e+08, +1.2808115049950667e+08, +5.3367146041461110e+07, +1.6420660320449572e+07, +3.5187129258106225e+06, +4.6916172344141634e+05, +2.9322607715088521e+04}
#define INIT_bcf_009 { -1.8154258054405320e+04, -3.2677664497929579e+05, -2.7776014823240140e+06, -1.4813874572394742e+07, -5.5552029646480277e+07, -1.5554568301014477e+08, -3.3701564652198035e+08, -5.7774110832339489e+08, -7.9439402394466794e+08, -8.8266002660518670e+08, -7.9439402394466794e+08, -5.7774110832339489e+08, -3.3701564652198035e+08, -1.5554568301014477e+08, -5.5552029646480277e+07, -1.4813874572394742e+07, -2.7776014823240140e+06, -3.2677664497929579e+05, -1.8154258054405320e+04}
#define INIT_bcf_010 { +6.5274731672945736e+03, +1.3054946334589148e+05, +1.2402199017859690e+06, +7.4413194107158137e+06, +3.1625607495542210e+07, +1.0120194398573507e+08, +2.5300485996433768e+08, +5.0600971992867535e+08, +8.2226579488409746e+08, +1.0963543931787965e+09, +1.2059898324966762e+09, +1.0963543931787965e+09, +8.2226579488409746e+08, +5.0600971992867535e+08, +2.5300485996433768e+08, +1.0120194398573507e+08, +3.1625607495542210e+07, +7.4413194107158137e+06, +1.2402199017859690e+06, +1.3054946334589148e+05, +6.5274731672945736e+03}
#define INIT_bcf_011 { -1.0350067135477861e+03, -2.2770147698051296e+04, -2.3908655082953861e+05, -1.5939103388635907e+06, -7.5710741096020555e+06, -2.7255866794567399e+07, -7.7224955917940974e+07, -1.7651418495529366e+08, -3.3096409679117560e+08, -5.1483303945293981e+08, -6.6928295128882170e+08, -7.3012685595144188e+08, -6.6928295128882170e+08, -5.1483303945293981e+08, -3.3096409679117560e+08, -1.7651418495529366e+08, -7.7224955917940974e+07, -2.7255866794567399e+07, -7.5710741096020555e+06, -1.5939103388635907e+06, -2.3908655082953861e+05, -2.2770147698051296e+04, -1.0350067135477861e+03}
#define INIT_bcf {INIT_bcf_000, INIT_bcf_001, INIT_bcf_002, INIT_bcf_003, INIT_bcf_004, INIT_bcf_005, INIT_bcf_006, INIT_bcf_007, INIT_bcf_008, INIT_bcf_009, INIT_bcf_010, INIT_bcf_011}

#define bcf_K 11
#define bcf_T 15.0

const double bcf[500][500] = INIT_bcf;


class bernoulli_covariates_fast
{
    public:

    class network
    {
        public:
        /* Here you should put all the variable which define the network
         * for scalar model, a mat field is good choice.
         * See bernoulli, poisson model for example.
         *
         * Covariates must also stocked here.
         * See poisson_covariates model for example
         */

        mat adj;
        cube covariates;


        /* Here you should add all precomputed values which depends only on the
         * network usefull in various functions, to avoid computing these value
         * many time.
         */

        mat adjD;
        mat adjZ;
        mat adjDZ;
        mat Mones;
        mat MonesZ;
        
        network(Rcpp::List & network_from_R)
        {
            /* Here you must define how initialize the fields describing the
             * network. Precomputed value must be computed here.
             * See bernoulli, poisson or poisson_covariates for example.
             *
             * For scalar network, the provided list have a adjacency field
             * which contains the adjacency matrix.
             * For scalar network with covariates vectors on edges (existent or
             * not), the provided list have a covariates field which contain a
             * list of matrices. This matrices have the same size than the
             * adjacency matrix, and the i-th matrix is the matrix of the i-th
             * covariate on all edges.
             */
            adj = Rcpp::as<mat>(network_from_R["adjacency"]);

            Rcpp::List covariates_list = network_from_R["covariates"];
            covariates.set_size(adj.n_rows,adj.n_cols,covariates_list.size());
            for(int k=0; k<covariates_list.size(); k++)
                covariates.slice(k) = Rcpp::as<mat>(covariates_list[k]);

            adjD = adj - .5;
            adjZ = fill_diag(adj,0);
            adjDZ = fill_diag(adjD,0);
            Mones = ones<mat>(adj.n_rows,adj.n_cols);
            MonesZ = fill_diag(Mones,0);
        }
    };

    // parameters
    unsigned int n_parameters;
    bool symmetric;

    /* Here you must put the bernoulli_covariates_fast model parameters */
    mat m;
    vec beta;

    bernoulli_covariates_fast(SBM & membership, bernoulli_covariates_fast::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM case
         */

        m = (membership.Z.t() * net.adjZ * membership.Z)
            /
            (membership.Z.t() * net.MonesZ * membership.Z);

        m=logit(m);

        beta.set_size(net.covariates.n_slices);
        beta.fill(0);

        n_parameters = m.n_elem + beta.n_elem;
        symmetric=false;
    }

    
    bernoulli_covariates_fast(SBM_sym & membership, bernoulli_covariates_fast::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM_sym case
         */

        m = (membership.Z.t() * net.adjZ * membership.Z)
            /
            (membership.Z.t() * net.MonesZ * membership.Z);

        m=logit(m);

        beta.set_size(net.covariates.n_slices);
        beta.fill(0);

        n_parameters = m.n_rows*(m.n_rows+1)/2 + beta.n_elem;
        symmetric=true;
    }

    
    bernoulli_covariates_fast(LBM & membership, bernoulli_covariates_fast::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * LBM case 
         */
        m = (membership.Z1.t() * net.adj * membership.Z2)
            /
            (membership.Z1.t() * net.Mones * membership.Z2);

        m=logit(m);

        beta.set_size(net.covariates.n_slices);
        beta.fill(0);

        n_parameters = m.n_rows * m.n_cols + net.covariates.n_slices;
        symmetric=false;
    }


    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["n_parameters"] = n_parameters;

        /* Here you must define the export way of your model parameters in the
         * list values. This is this list which is returned to R. */

        values["m"] = m;
        values["beta"] = beta;


        return values;
    }

    /* Keep this. When this is not usefull, this is not used. */
    inline vec to_vector();
    bernoulli_covariates_fast(SBM &, const vec &);
    bernoulli_covariates_fast(SBM_sym &, const vec &);
    bernoulli_covariates_fast(LBM &, const vec &);
};





/* Function for the model.
 */

/******************************************************************************
 *
 * Notations :
 *
 * i a row node
 * j a col node
 * q a row class
 * l a col class
 *
 * (for SBM there are no distinctions between row and col nodes and classes)
 *
 * f(i,j,q,l) the prob or density of the variable
 * X_{ij} conditionnally to i belong to class q and j to class l
 *
 * logf(i,j,q,l) is a quicker way to write log(f(i,j,q,l))
 *
 ******************************************************************************/





/* Usefull, Optional
 *
 * If you have direct formula to compute the following form, use it, else skip
 * this definition and delete this definition.
 *
 * forall i,q ; lZ_{iq} += \sum_{jl} Z_{jl} (logf(i,j,q,l) + logf(j,i,l,q)
 *
 * This is usefull for SBM. See bernoulli and poisson for example.
 */

template<>
inline
void e_fixed_step(SBM & membership,
                  bernoulli_covariates_fast & model,
                  bernoulli_covariates_fast::network & net,
                  mat & lZ)
{
    std::vector<mat> power_BoT(2*bcf_K+1);
    std::vector<mat> power_moT(2*bcf_K+1);

    power_BoT[0] = net.MonesZ;
    power_moT[0] = ones<mat>(lZ.n_cols,lZ.n_cols);
    power_BoT[1] = fill_diag(compute_B(model.beta, net.covariates)/bcf_T,0);
    power_moT[1] = model.m/bcf_T;
    for(unsigned int t=2;t<=2*bcf_K;t++)
    {
        power_BoT[t]=power_BoT[t-1] % power_BoT[1];
        power_moT[t]=power_moT[t-1] % power_moT[1];
    }
    
    lZ += net.adjDZ * membership.Z * model.m.t();
        + net.adjDZ.t() * membership.Z * model.m;

    for(unsigned int k=0; k<=bcf_K; k++)
    {
        for(unsigned int t=0; t<=2*k; t++)
        {
            lZ += bcf[k][t] * power_BoT[2*k-t] * membership.Z * power_moT[t].t();
                + bcf[k][t] * power_BoT[2*k-t].t() * membership.Z * power_moT[t];
        }
    }
}





/* Usefull, Optional
 *
 * If you have direct formula to compute the following form, use it, else skip
 * this definition and delete this definition.
 *
 * forall i,q ; lZ_{iq} += \sum_{jl} Z_{jl} logf(i,j,q,l)
 *
 * This is usefull for SBM_sym. See bernoulli and poisson for example.
 */

template<>
inline
void e_fixed_step(SBM_sym & membership,
                  bernoulli_covariates_fast & model,
                  bernoulli_covariates_fast::network & net,
                  mat & lZ)
{
    std::vector<mat> power_BoT(2*bcf_K+1);
    std::vector<mat> power_moT(2*bcf_K+1);

    power_BoT[0] = net.MonesZ;
    power_moT[0] = ones<mat>(lZ.n_cols,lZ.n_cols);
    power_BoT[1] = fill_diag(compute_B(model.beta, net.covariates)/bcf_T,0);
    power_moT[1] = model.m/bcf_T;
    for(unsigned int t=2;t<=2*bcf_K;t++)
    {
        power_BoT[t]=power_BoT[t-1] % power_BoT[1];
        power_moT[t]=power_moT[t-1] % power_moT[1];
    }


    lZ += net.adjDZ * membership.Z * model.m;

    for(unsigned int k=0; k<=bcf_K; k++)
    {
        for(unsigned int t=0; t<=2*k; t++)
        {
            lZ += bcf[k][t] * power_BoT[2*k-t] * power_moT[t];
        }
    }
}





/* Usefull, Optional
 *
 * If you have direct formula to compute the following form, use it, else skip
 * this definition and delete this definition.
 *
 * forall i,q ; lZ1_{iq} += \sum_{jl} Z2_{jl} logf(i,j,q,l)
 * forall j,l ; lZ2_{jl} += \sum_{iq} Z1_{iq} logf(i,j,q,l)
 *
 * This is usefull for LBM. See bernoulli and poisson for example.
 */

template<>
inline
void e_fixed_step(LBM & membership,
                  bernoulli_covariates_fast & model,
                  bernoulli_covariates_fast::network & net,
                  mat & lZ1,
                  mat & lZ2)
{
    std::vector<mat> power_BoT(2*bcf_K+1);
    std::vector<mat> power_moT(2*bcf_K+1);

    power_BoT[0] = net.Mones;
    power_moT[0] = ones<mat>(lZ1.n_cols,lZ2.n_cols);
    power_BoT[1] = compute_B(model.beta, net.covariates)/bcf_T;
    power_moT[1] = model.m/bcf_T;
    for(unsigned int t=2;t<=2*bcf_K;t++)
    {
        power_BoT[t]=power_BoT[t-1] % power_BoT[1];
        power_moT[t]=power_moT[t-1] % power_moT[1];
    }


    lZ1 += net.adjD * membership.Z2 * model.m.t();
    lZ2 += net.adjD.t() * membership.Z1 * model.m;

    for(unsigned int k=0; k<=bcf_K; k++)
    {
        for(unsigned int t=0; t<=2*k; t++)
        {
            lZ1 += bcf[k][t] * power_BoT[2*k-t] * membership.Z2 * power_moT[t].t();
            lZ2 += bcf[k][t] * power_BoT[2*k-t].t() * membership.Z1 * power_moT[t];
        }
    }
}





/* Usefull, Optional
 *
 * If you have a explicit formula to compute the maximum of PL in respect to
 * model parameters, knowing the membership, use it, else skip and delete this
 * definition. In the case of SBM.
 *
 * This function must return PL in the maximum.
 * 
 * This is usefull for SBM. See bernoulli and poisson for example.
 */

// template<>
// inline
// double m_step(SBM & membership,
//               bernoulli_covariates_fast & model,
//               bernoulli_covariates_fast::network & net)
// {
// }





/* Usefull, Optional
 *
 * If you have a explicit formula to compute the maximum of PL in respect to
 * model parameters, knowing the membership, use it, else skip and delete this
 * definition. In the case of SBM_sym.
 *
 * This function must return PL in the maximum.
 * 
 * This is usefull for SBM_sym. See bernoulli and poisson for example.
 */

// template<>
// inline
// double m_step(SBM_sym & membership,
//               bernoulli_covariates_fast & model,
//               bernoulli_covariates_fast::network & net)
// {
// }





/* Usefull, Optional
 *
 * If you have a explicit formula to compute the maximum of PL in respect to
 * model parameters, knowing the membership, use it, else skip and delete this
 * definition. In the case of LBM.
 *
 * This function must return PL in the maximum.
 * 
 * This is usefull for LBM. See bernoulli and poisson for example.
 */

// template<>
// inline
// double m_step(LBM & membership,
//               bernoulli_covariates_fast & model,
//               bernoulli_covariates_fast::network & net)
// {
// }





/* If you have defined m_step(SBM &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have direct formula to compute PL, without loops over (i,j,q,l), in
 * the case of SBM, use it.
 *
 * See poisson_covariate for example.
 */

template<>
inline
double PL(bernoulli_covariates_fast & model,
          SBM & membership,
          bernoulli_covariates_fast::network & net)
{
    std::vector<mat> power_BoT(2*bcf_K+1);
    std::vector<mat> power_moT(2*bcf_K+1);

    power_BoT[0] = net.MonesZ;
    power_moT[0] = ones<mat>(membership.Z.n_cols,membership.Z.n_cols);
    power_BoT[1] = fill_diag(compute_B(model.beta, net.covariates)/bcf_T,0);
    power_moT[1] = model.m/bcf_T;
    for(unsigned int t=2;t<=2*bcf_K;t++)
    {
        power_BoT[t]=power_BoT[t-1] % power_BoT[1];
        power_moT[t]=power_moT[t-1] % power_moT[1];
    }

    std::vector<mat> Z_power_BoT_Z(2*bcf_K+1);
    for(unsigned int t=0;t<=2*bcf_K;t++)
    {
        Z_power_BoT_Z[t] = membership.Z.t() * power_BoT[t] * membership.Z;
    }
    
    double value = 0;
    value += accu( (membership.Z.t() * net.adjDZ * membership.Z) % model.m );
    value += bcf_T * accu( 
                membership.Z.t() * (net.adjD % power_BoT[1]) * membership.Z
            );

    for(unsigned int k=0; k<=bcf_K; k++)
    {
        for(unsigned int t=0; t<=2*k; t++)
        {
            value += bcf[k][t] * accu(
                    ( Z_power_BoT_Z[2*k-t] )
                    %
                    ( power_moT[t] )
                    );
        }
    }
    return(value);
}






/* If you have defined m_step(SBM_sym &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have direct formula to compute PL, without loops over (i,j,q,l), in
 * the case of SBM_sym, use it.
 *
 * See poisson_covariate for example.
 */

template<>
inline
double PL(bernoulli_covariates_fast & model,
          SBM_sym & membership,
          bernoulli_covariates_fast::network & net)
{
    std::vector<mat> power_BoT(2*bcf_K+1);
    std::vector<mat> power_moT(2*bcf_K+1);

    power_BoT[0] = net.MonesZ;
    power_moT[0] = ones<mat>(membership.Z.n_cols,membership.Z.n_cols);
    power_BoT[1] = fill_diag(compute_B(model.beta, net.covariates)/bcf_T,0);
    power_moT[1] = model.m/bcf_T;
    for(unsigned int t=2;t<=2*bcf_K;t++)
    {
        power_BoT[t]=power_BoT[t-1] % power_BoT[1];
        power_moT[t]=power_moT[t-1] % power_moT[1];
    }
    
    std::vector<mat> Z_power_BoT_Z(2*bcf_K+1);
    for(unsigned int t=0;t<=2*bcf_K;t++)
    {
        Z_power_BoT_Z[t] = membership.Z.t() * power_BoT[t] * membership.Z;
    }
    
    double value = 0;
    value += accu( (membership.Z.t() * net.adjDZ * membership.Z) % model.m );
    value += bcf_T * accu( 
                membership.Z.t() * (net.adjD % power_BoT[1]) * membership.Z
            );

    for(unsigned int k=0; k<=bcf_K; k++)
    {
        for(unsigned int t=0; t<=2*k; t++)
        {
            value += bcf[k][t] * accu(
                    ( Z_power_BoT_Z[2*k-t] )
                    %
                    ( power_moT[t] )
                    );
        }
    }
    return(.5 * value);
}




/* If you have defined m_step(LBM &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have direct formula to compute PL, without loops over (i,j,q,l), in
 * the case of LBM, use it.
 *
 * See poisson_covariate for example.
 */

template<>
inline
double PL(bernoulli_covariates_fast & model,
          LBM & membership,
          bernoulli_covariates_fast::network & net)
{
    std::vector<mat> power_BoT(2*bcf_K+1);
    std::vector<mat> power_moT(2*bcf_K+1);

    power_BoT[0] = net.Mones;
    power_moT[0] = ones<mat>(membership.Z1.n_cols,membership.Z2.n_cols);
    power_BoT[1] = compute_B(model.beta, net.covariates)/bcf_T;
    power_moT[1] = model.m/bcf_T;
    for(unsigned int t=2;t<=2*bcf_K;t++)
    {
        power_BoT[t]=power_BoT[t-1] % power_BoT[1];
        power_moT[t]=power_moT[t-1] % power_moT[1];
    }
    
    std::vector<mat> Z1_power_BoT_Z2(2*bcf_K+1);
    for(unsigned int t=0;t<=2*bcf_K;t++)
    {
        Z1_power_BoT_Z2[t] = membership.Z1.t() * power_BoT[t] * membership.Z2;
    }
    

    double value = 0;
    value += accu( (membership.Z1.t() * net.adjD * membership.Z2) % model.m );
    value += bcf_T * accu( 
                membership.Z1.t() * (net.adjD % power_BoT[1]) * membership.Z2
            );

    for(unsigned int k=0; k<=bcf_K; k++)
    {
        for(unsigned int t=0; t<=2*k; t++)
        {
            value += bcf[k][t] * accu(
                    ( Z1_power_BoT_Z2[2*k-t] )
                    %
                    ( power_moT[t] )
                    );
        }
    }
    return(value);
}
            






/* If one of the following specialization
 *     - m_step(LBM &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network)
 *     - m_step(SBM &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network)
 *     - m_step(SBM_sym &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network)
 *   is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * Export the model parameters as a vector.
 *
 * See naive_bernoulli and poisson_covariates for example.
 */

inline
vec bernoulli_covariates_fast::to_vector()
{
    vec out(n_parameters);
    vec vm = (symmetric) ? vech(m) : reshape(m,m.n_elem,1);

    out.subvec(0,vm.n_elem-1) = vm;
    out.subvec(vm.n_elem,n_parameters-1) = beta;
    return(out);
}





/* If m_step(SBM &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network) is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 * 
 * Constructor of the model, from a vector of parameter named vectorized with
 * the given membership of type SBM.
 *
 * See naive_bernoulli for example.
 */

bernoulli_covariates_fast::bernoulli_covariates_fast(SBM & membership, const vec & vectorized)
{
    unsigned int Q=membership.Z.n_cols;
    m = reshape(vectorized.subvec(0,Q*Q-1),Q,Q);
    beta = vectorized.subvec(Q*Q,vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric = false;
}





/* If m_step(SBM_sym &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network) is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 * 
 * Constructor of the model, from a vector of parameter named vectorized with
 * the given membership of type SBM_sym.
 *
 * See naive_bernoulli for example.
 */

bernoulli_covariates_fast::bernoulli_covariates_fast(SBM_sym & membership, const vec & vectorized)
{
    unsigned int Q=membership.Z.n_cols;
    m = unvech(vectorized.subvec(0,Q*(Q+1)/2-1));
    beta = vectorized.subvec(Q*(Q+1)/2,vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric = true;
}





/* If m_step(LBM &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network) is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 * 
 * Constructor of the model, from a vector of parameter named vectorized with
 * the given membership of type LBM.
 *
 * See naive_bernoulli for example.
 */

bernoulli_covariates_fast::bernoulli_covariates_fast(LBM & membership, const vec & vectorized)
{
    unsigned int Q1=membership.Z1.n_cols;
    unsigned int Q2=membership.Z2.n_cols;
    m = reshape(vectorized.subvec(0,Q1*Q2-1),Q1,Q2);
    beta = vectorized.subvec(Q1*Q2,vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric=false;
}





/* If one of the following specialization
 *     - m_step(LBM &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network)
 *     - m_step(SBM &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network)
 *     - m_step(SBM_sym &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network)
 *   is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * From the given model, in the given direction for model parameters, return the
 * maximum betwenn 1 and the maximum step the descend algorithm can do.
 *
 * See naive_bernoulli for example.
 *
 * If you don't have any constraint on you parameter vector, return 1, else, use
 * your brain.
 *
 * You can specialize this template in three function for SBM, LBM, and SBM_sym
 * if you need
 */

template<class membership_type>
inline
double maximum_step_in_direction(membership_type & membership,
                                 bernoulli_covariates_fast & model,
                                 bernoulli_covariates_fast::network & net,
                                 vec & direction)
{
    return 1;
}





/* If you have defined m_step(SBM &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have a direct formula to compute the gradient of PL in respect to your
 * parameters vectors without loops over (i,j,q,l), in the case of SBM, use it.
 *
 * See poisson_covariates for example.
 */

template<>
inline
vec grad(bernoulli_covariates_fast & model,
         SBM & membership,
         bernoulli_covariates_fast::network & net)
{
    std::vector<mat> power_BoT(2*bcf_K+1);
    std::vector<mat> power_moT(2*bcf_K+1);

    power_BoT[0] = net.MonesZ;
    power_moT[0] = ones<mat>(membership.Z.n_cols,membership.Z.n_cols);
    power_BoT[1] = fill_diag(compute_B(model.beta, net.covariates)/bcf_T,0);
    power_moT[1] = model.m/bcf_T;
    for(unsigned int t=2;t<=2*bcf_K;t++)
    {
        power_BoT[t]=power_BoT[t-1] % power_BoT[1];
        power_moT[t]=power_moT[t-1] % power_moT[1];
    }

    std::vector<mat> Z_power_BoT_Z(2*bcf_K+1);
    std::vector<mat> Z_power_moT_Z(2*bcf_K+1);

    for(unsigned int t=0; t<2*bcf_K;t++) // the last is never used
    {
        Z_power_BoT_Z[t]=membership.Z.t() * power_BoT[t] * membership.Z;
        Z_power_moT_Z[t]=membership.Z * power_moT[t] * membership.Z.t();
    }

    mat dm = membership.Z.t() * net.adjDZ * membership.Z;

    mat mask = (membership.Z * power_moT[0] * membership.Z.t() ) % net.adjDZ;

    for(unsigned int k=0;k<=bcf_K;k++)
    {
        for(unsigned int t=1;t<=2*k;t++) // t=1..
        {
            dm += bcf[k][t] * t / bcf_T * 
                (
                 power_moT[t-1] 
                 % 
                 (Z_power_BoT_Z[2*k-t])
                );
            mask += bcf[k][t] * t / bcf_T *
                (
                 power_BoT[t-1] 
                 % 
                 (Z_power_moT_Z[2*k-t])
                );
        }
    }

    vec dbeta(model.beta.n_elem);
    for(unsigned int k=0;k<dbeta.n_elem;k++)
    {
        dbeta(k) = accu( mask % net.covariates.slice(k));
    }

    vec out(model.n_parameters);
    vec vdm = reshape(dm,dm.n_elem,1);

    out.subvec(0,vdm.n_elem-1) = vdm;
    out.subvec(vdm.n_elem,model.n_parameters-1) = dbeta;

    return(out);
}





/* If you have defined m_step(SBM_sym &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have a direct formula to compute the gradient of PL in respect to your
 * parameters vectors without loops over (i,j,q,l), in the case of SBM, use it.
 *
 * See poisson_covariates for example.
 */

template<>
inline
vec grad(bernoulli_covariates_fast & model,
         SBM_sym & membership,
         bernoulli_covariates_fast::network & net)
{
    std::vector<mat> power_BoT(2*bcf_K+1);
    std::vector<mat> power_moT(2*bcf_K+1);

    power_BoT[0] = net.MonesZ;
    power_moT[0] = ones<mat>(membership.Z.n_cols,membership.Z.n_cols);
    power_BoT[1] = fill_diag(compute_B(model.beta, net.covariates)/bcf_T,0);
    power_moT[1] = model.m/bcf_T;
    for(unsigned int t=2;t<=2*bcf_K;t++)
    {
        power_BoT[t]=power_BoT[t-1] % power_BoT[1];
        power_moT[t]=power_moT[t-1] % power_moT[1];
    }

    std::vector<mat> Z_power_BoT_Z(2*bcf_K+1);
    std::vector<mat> Z_power_moT_Z(2*bcf_K+1);

    for(unsigned int t=0; t<2*bcf_K;t++) // the last is never used
    {
        Z_power_BoT_Z[t]=membership.Z.t() * power_BoT[t] * membership.Z;
        Z_power_moT_Z[t]=membership.Z * power_moT[t] * membership.Z.t();
    }
    
    mat dm = membership.Z.t() * net.adjDZ * membership.Z;

    mat mask = (membership.Z * power_moT[0] * membership.Z.t() ) % net.adjDZ;

    for(unsigned int k=0;k<=bcf_K;k++)
    {
        for(unsigned int t=1;t<=2*k;t++) // t=1..
        {
            dm += bcf[k][t] * t / bcf_T * 
                (
                 power_moT[t-1] 
                 % 
                 ( Z_power_BoT_Z[2*k-t] )
                );
            mask += bcf[k][t] * t / bcf_T *
                (
                 power_BoT[t-1] 
                 % 
                 ( Z_power_moT_Z[2*k-t] )
                );
        }
    }

    vec dbeta(model.beta.n_elem);
    for(unsigned int k=0;k<dbeta.n_elem;k++)
    {
        dbeta(k) = accu( mask % net.covariates.slice(k));
    }

    vec out(model.n_parameters);
    vec vdm = vech(dm);

    out.subvec(0,vdm.n_elem-1) = vdm;
    out.subvec(vdm.n_elem,model.n_parameters-1) = dbeta;

    return(.5*out);
}





/* If you have defined m_step(LBM &, bernoulli_covariates_fast &, bernoulli_covariates_fast::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have a direct formula to compute the gradient of PL in respect to your
 * parameters vectors without loops over (i,j,q,l), in the case of LBM, use it.
 *
 * See poisson_covariates for example.
 */

template<>
inline
vec grad(bernoulli_covariates_fast & model,
         LBM & membership,
         bernoulli_covariates_fast::network & net)
{
    std::vector<mat> power_BoT(2*bcf_K+1);
    std::vector<mat> power_moT(2*bcf_K+1);

    power_BoT[0] = net.Mones;
    power_moT[0] = ones<mat>(membership.Z1.n_cols,membership.Z2.n_cols);
    power_BoT[1] = compute_B(model.beta, net.covariates)/bcf_T;
    power_moT[1] = model.m/bcf_T;
    for(unsigned int t=2;t<=2*bcf_K;t++)
    {
        power_BoT[t]=power_BoT[t-1] % power_BoT[1];
        power_moT[t]=power_moT[t-1] % power_moT[1];
    }
    
    std::vector<mat> Z1_power_BoT_Z2(2*bcf_K+1);
    std::vector<mat> Z1_power_moT_Z2(2*bcf_K+1);

    for(unsigned int t=0; t<2*bcf_K;t++) // the last is never used
    {
        Z1_power_BoT_Z2[t]=membership.Z1.t() * power_BoT[t] * membership.Z2;
        Z1_power_moT_Z2[t]=membership.Z1 * power_moT[t] * membership.Z2.t();
    }

    mat dm = membership.Z1.t() * net.adjD * membership.Z2;

    mat mask = (membership.Z1 * power_moT[0] * membership.Z2.t() ) % net.adjD;

    for(unsigned int k=0;k<=bcf_K;k++)
    {
        for(unsigned int t=1;t<=2*k;t++) // t=1..
        {
            dm += bcf[k][t] * t / bcf_T * 
                (
                 power_moT[t-1] 
                 % 
                 ( Z1_power_BoT_Z2[2*k-t] )
                );
            mask += bcf[k][t] * t / bcf_T *
                (
                 power_BoT[t-1] 
                 % 
                 ( Z1_power_moT_Z2[2*k-t] )
                );
        }
    }

    vec dbeta(model.beta.n_elem);
    for(unsigned int k=0;k<dbeta.n_elem;k++)
    {
        dbeta(k) = accu( mask % net.covariates.slice(k));
    }

    vec out(model.n_parameters);
    vec vdm = reshape(dm,dm.n_elem,1);

    out.subvec(0,vdm.n_elem-1) = vdm;
    out.subvec(vdm.n_elem,model.n_parameters-1) = dbeta;

    return(out);
}





/* If one of the following specialization
 *     - grad(bernoulli_covariates_fast &, LBM &, bernoulli_covariates_fast::network)
 *     - grad(bernoulli_covariates_fast &, SBM &, bernoulli_covariates_fast::network)
 *     - grad(bernoulli_covariates_fast &, SBM_sym &, bernoulli_covariates_fast::network)
 *   is Usefull and not defined
 * Then
 *     Usefull, Optional
 * Else
 *     Useless
 *
 * The gradient of logf(i,j,q,l) in respect to the parameter vector.
 */

// template<>
// inline
// vec grad_logf(bernoulli_covariates_fast & model,
//               bernoulli_covariates_fast::network & net,
//               unsigned int i,
//               unsigned int j,
//               unsigned int q,
//               unsigned int l)
// {
// }





/* If grad_logf(bernoulli_covariates_fast &,
 *         bernoulli_covariates_fast::network &,
 *         unsigned int,
 *         unsigned int,
 *         unsigned int,
 *         unsigned int)
 *    is Usefull and not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 *
 * The derivative of logf(i,j,q,l) in respect to the k-th parameters of the
 * parameter vector.
 *
 * See naive_bernoulli for example.
 */

// inline
// double grad_logf(bernoulli_covariates_fast & model,
//                  bernoulli_covariates_fast::network & net,
//                  unsigned int i,
//                  unsigned int j,
//                  unsigned int q,
//                  unsigned int l,
//                  unsigned int k)
// {
// }





/* If one of the following specialization
 *     - PL(bernoulli_covariates_fast &, SBM &, bernoulli_covariates_fast::network)
 *     - PL(bernoulli_covariates_fast &, SBM_sym &, bernoulli_covariates_fast::network)
 *     - PL(bernoulli_covariates_fast &, LBM &, bernoulli_covariates_fast::network)
 *   is Usefull and not definied
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * logf(i,j,q,l) see notations.
 */

// inline
// double logf(bernoulli_covariates_fast & model,
//             bernoulli_covariates_fast::network & net,
//             unsigned int i,
//             unsigned int j,
//             unsigned int q,
//             unsigned int l)
// {
// }
