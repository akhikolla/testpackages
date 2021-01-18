
struct SBM
{
    mat Z;
    rowvec alpha;

    SBM(Rcpp::List & membership_from_R)
    {
        mat origZ = membership_from_R["Z"];
        Z=origZ;
        double tol = TOL_P/Z.n_rows;
        boundaries(Z,tol,1-tol);
        Z /= repmat( sum(Z,1), 1, Z.n_cols );
        alpha = sum(Z,0) / Z.n_rows;
    }

    SBM& operator=(const SBM& orig)
    {
        Z=orig.Z;
        alpha=orig.alpha;
        
        return *this;
    }
    
    inline
    double entropy()
    {
        return - accu( Z % log(Z) );
    }

    template<class model_type, class network_type>
    inline
    void e_step(model_type & model, network_type & net)
    {
        double tol = TOL_P/Z.n_rows;
        double step_size;
        unsigned int niter=0;
        do
        {
            // lZ the new log(Z) without renormalization
            mat lZ = repmat(log(alpha),Z.n_rows,1);
            
            // with a template, should be specialized by the model if possible
            e_fixed_step(*this, model, net, lZ);
            
            // This operation should change nothing theorically (but
            // pratically...)
            lZ -= repmat( max(lZ,1), 1, lZ.n_cols );

            // After this lZ is not the new log(Z) but the new Z without
            // normalization
            lZ = exp(lZ);

            // normalization
            lZ /= repmat( sum(lZ,1), 1, lZ.n_cols );

            boundaries(lZ,tol,1-tol);
            lZ /= repmat( sum(lZ,1), 1, lZ.n_cols );

            step_size = max(max(abs( Z-lZ )));

            Z=lZ;

            niter++;

            #ifdef DEBUG_E
                fprintf(stderr,"E iteration: %i %f\n",niter,step_size);
            #endif

        } while( step_size>TOL_F && niter<10);
    }

    inline
    double m_step()
    {
        alpha = sum(Z,0) / Z.n_rows;
        return accu(Z * log(alpha).t());
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["Z"] = Z;
        values["alpha"] = alpha;

        return values;
    }
        
};

struct SBM_sym : public SBM
{
    SBM_sym(Rcpp::List & membership_from_R) : SBM(membership_from_R) {}
};



