
class poisson
{
    public:

    class network
    {
        public:
        mat adj; // adjacency matrix

        // precalculated matrices for SBM
        mat adjZD;
        mat adjZDt;
        mat MonesZD;

        // precalculated matrices for LBM
        mat Mones;
        mat Monest;
        mat adjt;

        // precalculated value for liklihood
        double accu_log_fact_XZD;
        double accu_log_fact_X;


        network(Rcpp::List & network_from_R)
        {
            mat adj_orig = network_from_R["adjacency"];

            adj = adj_orig;
            adjt = adj_orig.t();
            Mones = ones<mat>(adj.n_rows,adj.n_cols);
            Monest = Mones.t();
            adjZD = fill_diag(adj_orig,0);
            adjZDt = adjZD.t();
            MonesZD = fill_diag(Mones,0);

            accu_log_fact(adj,accu_log_fact_X,accu_log_fact_XZD);

        }
    };

    /* parameters */
    unsigned int n_parameters;
    mat lambda;

    poisson(SBM & membership, poisson::network & net)
    {
        n_parameters = membership.Z.n_cols * membership.Z.n_cols;
        lambda.set_size(membership.Z.n_cols,membership.Z.n_cols);
    }
    
    poisson(SBM_sym & membership, poisson::network & net)
    {
        n_parameters = membership.Z.n_cols * (membership.Z.n_cols+1)/2;
        lambda.set_size(membership.Z.n_cols,membership.Z.n_cols);
    }
    
    poisson(LBM & membership, poisson::network & net)
    {
        n_parameters = membership.Z1.n_cols * membership.Z2.n_cols;
        lambda.set_size(membership.Z1.n_cols,membership.Z2.n_cols);
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["lambda"] = lambda;
        values["n_parameters"] = n_parameters;

        return values;
    }
};

template<>
inline
void e_fixed_step(SBM & membership, poisson & model, poisson::network & net, mat & lZ)
{
    lZ+= net.adjZD * membership.Z * log(model.lambda).t()
       - net.MonesZD * membership.Z * model.lambda.t()
       + net.adjZDt * membership.Z * log(model.lambda)
       - net.MonesZD * membership.Z * model.lambda;
}

template<>
inline
void e_fixed_step(SBM_sym & membership, poisson & model, poisson::network & net, mat & lZ)
{
    lZ+= net.adjZD * membership.Z * log(model.lambda).t()
       - net.MonesZD * membership.Z * model.lambda.t();
}

template<>
inline
void e_fixed_step(LBM & membership, poisson & model, poisson::network & net, mat & lZ1, mat & lZ2)
{
    lZ1 += net.adj * membership.Z2 * log(model.lambda).t()
         - net.Mones * membership.Z2 * model.lambda.t();
    lZ2 += net.adjt * membership.Z1 * log(model.lambda)
         - net.Monest * membership.Z1 * model.lambda;
}

template<>
inline
double m_step(SBM & membership, poisson & model, poisson::network & net)
{
    model.lambda = (membership.Z.t() * net.adjZD * membership.Z)
                    /
                   (membership.Z.t() * net.MonesZD * membership.Z);

    return
        (
            accu(
                
                    -(
                        model.lambda 
                        % 
                        (membership.Z.t() * net.MonesZD * membership.Z)
                    )
                    +
                    (
                        log(model.lambda)
                        %
                        (membership.Z.t() * net.adjZD * membership.Z)
                    )
                )
            - net.accu_log_fact_X
        );
}

template<>
inline
double m_step(SBM_sym & membership, poisson & model, poisson::network & net)
{
    return m_step<SBM>(membership,model,net)/2;
}

template<>
inline
double m_step(LBM & membership, poisson & model, poisson::network & net)
{
    model.lambda = (membership.Z1.t() * net.adj * membership.Z2)
                    /
                   (membership.Z1.t() * net.Mones * membership.Z2);

    return
        (
            accu(
                
                    -(
                        model.lambda 
                        % 
                        (membership.Z1.t() * net.Mones * membership.Z2)
                    )
                    +
                    (
                        log(model.lambda)
                        %
                        (membership.Z1.t() * net.adj * membership.Z2)
                    )
                )
            - net.accu_log_fact_X
        );
}


