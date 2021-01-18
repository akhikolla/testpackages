
/*
 * This is package skeleton for bernoulli_multiplex model
 *
 * You should delete this header when the bernoulli_multiplex model is ready for use to avoid
 * any confusion
 */

class bernoulli_multiplex
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
         *cube n*n*d where d is the number of multiplex compounant
         */
         cube adj;
         


        /* Here you should add all precomputed values which depends only on the
         * network usefull in various functions, to avoid computing these value
         * many time.
         */


         mat MonesZD;
         mat Mones;
         field<mat> adj_indicator;
         field<mat> adj_indicatorZD;


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
            Rcpp::List adj_list = network_from_R["adjacency"];
            mat first_mat = Rcpp::as<mat>( adj_list[0] );
            adj.set_size(first_mat.n_rows, first_mat.n_cols, adj_list.size());
            for(int k=0; k<adj_list.size(); k++)
                adj.slice(k) = Rcpp::as<mat>(adj_list[k]);

           // precomputation
           Mones = ones<mat>(first_mat.n_rows, first_mat.n_cols);
           MonesZD = fill_diag(Mones,0);
           
           unsigned int K=(1U << adj.n_slices);
           
           adj_indicator.set_size(K);
           adj_indicatorZD.set_size(K);
           

           for (unsigned int k=0;k<K;k++)
           {
                adj_indicator(k)=Mones;
                for (unsigned int s=0;s<adj.n_slices;s++)
                {
                    if  ((k>> (adj.n_slices-1-s)) & 1U )
                        adj_indicator(k)%= adj.slice(s); 
                    else 
                        adj_indicator(k)%= 1 - adj.slice(s);    

                }
                
                adj_indicatorZD(k)=fill_diag(adj_indicator(k),0);  

           }   
          
        }
    };

    // parameters
    unsigned int n_parameters;

    /* Here you must put the bernoulli_multiplex model parameters */
    field<mat> pi;
    
    unsigned int d; //pour avoir le nombre de reseaux

    bernoulli_multiplex(SBM & membership, bernoulli_multiplex::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM case
         */
         n_parameters =((1U << net.adj.n_slices)-1) * membership.Z.n_cols * membership.Z.n_cols; 
         pi.set_size(1U << net.adj.n_slices);
         d=net.adj.n_slices;
         for (unsigned int k=0; k<pi.n_elem;k++)
             pi(k).set_size(membership.Z.n_cols,membership.Z.n_cols);
    }
    
    bernoulli_multiplex(SBM_sym & membership, bernoulli_multiplex::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM_sym case
         */
         n_parameters = ((1U << net.adj.n_slices)-1) * membership.Z.n_cols*(membership.Z.n_cols+1)/2;
         pi.set_size(1U << net.adj.n_slices);
         d=net.adj.n_slices;
         for (unsigned int k=0; k<pi.n_elem;k++)
             pi(k).set_size(membership.Z.n_cols,membership.Z.n_cols);
    }
    
    bernoulli_multiplex(LBM & membership, bernoulli_multiplex::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * LBM case 
         */
         n_parameters = ((1U << net.adj.n_slices)-1) * membership.Z1.n_cols*membership.Z2.n_cols;
         pi.set_size(1U << net.adj.n_slices);
         d=net.adj.n_slices;
         for (unsigned int k=0; k<pi.n_elem;k++)
             pi(k).set_size(membership.Z1.n_cols,membership.Z2.n_cols);
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["n_parameters"] = n_parameters;

        /* Here you must define the export way of your model parameters in the
         * list values. This is this list which is returned to R. */
        Rcpp::List pi_as_list;
 

        for(unsigned int k=0;k<pi.n_elem;k++)
        {         
		
                char name[64];
                for(unsigned int t=0; t<64; name[t++]=0);
                for (unsigned int s=0;s<d;s++)
                {
                    if  ((k>> (d-1-s)) & 1U )
                        name[s]='1';
                    else
                        name[s]='0';
                }
	   
            pi_as_list[name]=pi(k);
        }
        values["pi"] = pi_as_list;
        
  
        return values;
    }

    /* Keep this. When this is not usefull, this is not used. */
    inline vec to_vector();
    bernoulli_multiplex(SBM &, const vec &);
    bernoulli_multiplex(SBM_sym &, const vec &);
    bernoulli_multiplex(LBM &, const vec &);
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
                  bernoulli_multiplex & model,
                   bernoulli_multiplex::network & net,
                   mat & lZ)
{
        for(unsigned int k=0;k<net.adj_indicator.n_elem;k++)
        {
           lZ+=net.adj_indicatorZD(k) * membership.Z * log(model.pi(k).t())
              + net.adj_indicatorZD(k).t() * membership.Z * log(model.pi(k));
  
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
                   bernoulli_multiplex & model,
                   bernoulli_multiplex::network & net,
                   mat & lZ)
{

        for(unsigned int k=0;k<net.adj_indicator.n_elem;k++)
        {
           lZ+=net.adj_indicatorZD(k) * membership.Z * log(model.pi(k).t());
              
  
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
                   bernoulli_multiplex & model,
                   bernoulli_multiplex::network & net,
                   mat & lZ1,
                   mat & lZ2)
{
        for(unsigned int k=0;k<net.adj_indicator.n_elem;k++)
        {
           lZ1+=net.adj_indicator(k) * membership.Z2 * log(model.pi(k).t());
           lZ2+=net.adj_indicator(k).t() * membership.Z1 * log(model.pi(k));
  
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

template<>
inline
double m_step(SBM & membership,
               bernoulli_multiplex & model,
               bernoulli_multiplex::network & net)
{
        
        double PL=0;

        for(unsigned int k=0;k<net.adj_indicator.n_elem;k++)
        {
               mat machin = membership.Z.t() * net.adj_indicatorZD(k) * membership.Z;
               model.pi(k) = (machin) / (membership.Z.t() * net.MonesZD * membership.Z ); 

               PL+=accu(machin % log(model.pi(k)));

        }

       return(PL);

}





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

template<>
inline
double m_step(SBM_sym & membership,
               bernoulli_multiplex & model,
               bernoulli_multiplex::network & net)
{
      
       return m_step<SBM>(membership,model,net)/2;


}





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

template<>
inline
double m_step(LBM & membership,
               bernoulli_multiplex & model,
               bernoulli_multiplex::network & net)
{

        double PL=0;

        for(unsigned int k=0;k<net.adj_indicator.n_elem;k++)
        {
               mat machin = membership.Z1.t() * net.adj_indicator(k) * membership.Z2;
               model.pi(k) = (machin) / (membership.Z1.t() * net.Mones * membership.Z2 ); 

               PL+=accu(machin % log(model.pi(k)));

        }

       return(PL);


}





/* If you have defined m_step(SBM &, bernoulli_multiplex &, bernoulli_multiplex::network &)
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

// template<>
// inline
// double PL(bernoulli_multiplex & model,
//           SBM & membership,
//           bernoulli_multiplex::network & net)





/* If you have defined m_step(SBM_sym &, bernoulli_multiplex &, bernoulli_multiplex::network &)
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

// template<>
// inline
// double PL(bernoulli_multiplex & model,
//           SBM_sym & membership,
//           bernoulli_multiplex::network & net)





/* If you have defined m_step(LBM &, bernoulli_multiplex &, bernoulli_multiplex::network &)
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

// template<>
// inline
// double PL(bernoulli_multiplex & model,
//           LBM & membership,
//           bernoulli_multiplex::network & net)





/* If one of the following specialization
 *     - m_step(LBM &, bernoulli_multiplex &, bernoulli_multiplex::network)
 *     - m_step(SBM &, bernoulli_multiplex &, bernoulli_multiplex::network)
 *     - m_step(SBM_sym &, bernoulli_multiplex &, bernoulli_multiplex::network)
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

// inline
// vec bernoulli_multiplex::to_vector()
// {
// }





/* If m_step(SBM &, bernoulli_multiplex &, bernoulli_multiplex::network) is not defined
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

// bernoulli_multiplex::bernoulli_multiplex(SBM & membership, const vec & vectorized)
// {
// }





/* If m_step(SBM_sym &, bernoulli_multiplex &, bernoulli_multiplex::network) is not defined
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

// bernoulli_multiplex::bernoulli_multiplex(SBM_sym & membership, const vec & vectorized)
// {
// }





/* If m_step(LBM &, bernoulli_multiplex &, bernoulli_multiplex::network) is not defined
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

// bernoulli_multiplex::bernoulli_multiplex(LBM & membership, const vec & vectorized)
// {
// }





/* If one of the following specialization
 *     - m_step(LBM &, bernoulli_multiplex &, bernoulli_multiplex::network)
 *     - m_step(SBM &, bernoulli_multiplex &, bernoulli_multiplex::network)
 *     - m_step(SBM_sym &, bernoulli_multiplex &, bernoulli_multiplex::network)
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

// template<class membership_type>
// inline
// double maximum_step_in_direction(membership_type & membership,
//                                  bernoulli_multiplex & model,
//                                  bernoulli_multiplex::network & net,
//                                  vec & direction)
// {
// }





/* If you have defined m_step(SBM &, bernoulli_multiplex &, bernoulli_multiplex::network &)
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

// template<>
// inline
// vec grad(bernoulli_multiplex & model,
//          SBM & membership,
//          bernoulli_multiplex::network & net)
// {
// }





/* If you have defined m_step(SBM_sym &, bernoulli_multiplex &, bernoulli_multiplex::network &)
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

// template<>
// inline
// vec grad(bernoulli_multiplex & model,
//          SBM_sym & membership,
//          bernoulli_multiplex::network & net)
// {
// }





/* If you have defined m_step(LBM &, bernoulli_multiplex &, bernoulli_multiplex::network &)
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

// template<>
// inline
// vec grad(bernoulli_multiplex & model,
//          LBM & membership,
//          bernoulli_multiplex::network & net)
// {
// }





/* If one of the following specialization
 *     - grad(bernoulli_multiplex &, LBM &, bernoulli_multiplex::network)
 *     - grad(bernoulli_multiplex &, SBM &, bernoulli_multiplex::network)
 *     - grad(bernoulli_multiplex &, SBM_sym &, bernoulli_multiplex::network)
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
// vec grad_logf(bernoulli_multiplex & model,
//               bernoulli_multiplex::network & net,
//               unsigned int i,
//               unsigned int j,
//               unsigned int q,
//               unsigned int l)
// {
// }





/* If grad_logf(bernoulli_multiplex &,
 *         bernoulli_multiplex::network &,
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
// double grad_logf(bernoulli_multiplex & model,
//                  bernoulli_multiplex::network & net,
//                  unsigned int i,
//                  unsigned int j,
//                  unsigned int q,
//                  unsigned int l,
//                  unsigned int k)
// {
// }





/* If one of the following specialization
 *     - PL(bernoulli_multiplex &, SBM &, bernoulli_multiplex::network)
 *     - PL(bernoulli_multiplex &, SBM_sym &, bernoulli_multiplex::network)
 *     - PL(bernoulli_multiplex &, LBM &, bernoulli_multiplex::network)
 *   is Usefull and not definied
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * logf(i,j,q,l) see notations.
 */

// inline
// double logf(bernoulli_multiplex & model,
//             bernoulli_multiplex::network & net,
//             unsigned int i,
//             unsigned int j,
//             unsigned int q,
//             unsigned int l)
// {
// }
