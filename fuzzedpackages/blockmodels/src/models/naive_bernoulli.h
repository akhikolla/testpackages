
class naive_bernoulli
{
    public:

    class network
    {
        public:
        mat adj;

        network(Rcpp::List & network_from_R)
        {
            mat adj_orig = network_from_R["adjacency"];
            adj = adj_orig;
        }
    };

    /* parameters */
    unsigned int n_parameters;
    bool symmetric;
    mat pi;

    naive_bernoulli(SBM & membership, naive_bernoulli::network & net)
    {
        n_parameters = membership.Z.n_cols * membership.Z.n_cols;
        pi.set_size(membership.Z.n_cols,membership.Z.n_cols);
        pi.fill(accu(net.adj)/(net.adj.n_rows*net.adj.n_cols));
        symmetric=false;
    }
    
    naive_bernoulli(SBM_sym & membership, naive_bernoulli::network & net)
    {
        n_parameters = membership.Z.n_cols * (membership.Z.n_cols + 1)/2;
        pi.set_size(membership.Z.n_cols,membership.Z.n_cols);
        pi.fill(accu(net.adj)/(net.adj.n_rows*net.adj.n_cols));
        symmetric=true;
    }
    
    naive_bernoulli(LBM & membership, naive_bernoulli::network & net)
    {
        n_parameters = membership.Z1.n_cols * membership.Z2.n_cols;
        pi.set_size(membership.Z1.n_cols,membership.Z2.n_cols);
        pi.fill(accu(net.adj)/(net.adj.n_rows*net.adj.n_cols));
        symmetric=false;
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["pi"] = pi;
        values["n_parameters"] = n_parameters;

        return values;
    }

    inline
    vec to_vector()
    {
        if(symmetric)
            return vech(pi);
        else
            return reshape(pi,n_parameters,1);
    }

    naive_bernoulli(SBM & membership, const vec & vectorized)
    {
        n_parameters = membership.Z.n_cols * membership.Z.n_cols;
        pi = reshape(vectorized, membership.Z.n_cols, membership.Z.n_cols);
    }
    
    naive_bernoulli(SBM_sym & membership, const vec & vectorized)
    {
        n_parameters = membership.Z.n_cols * (membership.Z.n_cols+1)/2;
        pi = unvech(vectorized);
    }
    
    naive_bernoulli(LBM & membership, const vec & vectorized)
    {
        n_parameters = membership.Z1.n_cols * membership.Z2.n_cols;
        pi = reshape(vectorized, membership.Z1.n_cols, membership.Z2.n_cols);
    }
};

inline
double logf(naive_bernoulli & model, naive_bernoulli::network & net, unsigned int i, unsigned int j, unsigned int q, unsigned int l)
{
    return net.adj(i,j)*log(model.pi(q,l)) + (1-net.adj(i,j))*log(1-model.pi(q,l));
}

inline
double grad_logf(naive_bernoulli & model, naive_bernoulli::network & net, unsigned int i, unsigned int j, unsigned int q, unsigned int l, unsigned int param)
{
    if( param % model.pi.n_rows == q )
    {
        if( param / model.pi.n_rows == l)
            return net.adj(i,j)/(model.pi(q,l)) - (1-net.adj(i,j))/(1-model.pi(q,l));
    }
    
    return 0;
}

template<class membership_type>
inline
double maximum_step_in_direction(membership_type & membership, naive_bernoulli & model, naive_bernoulli::network & net, vec & direction)
{
    vec v = model.to_vector();

    double amax=1;

    for(unsigned int k=0; k<direction.n_elem; k++)
    {
        double b = direction(k)>0 ?
                                    (1-v(k))/direction(k)
                                  :
                                    -v(k)/direction(k);
        if(b<amax)
            amax=b;
    }

    return amax;            
    
}
