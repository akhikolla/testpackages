
struct LBM
{
    mat Z1;
    mat Z2;
    rowvec alpha1;
    rowvec alpha2;
    
    LBM(Rcpp::List & membership_from_R)
    {
        mat origZ1 = membership_from_R["Z1"];
        mat origZ2 = membership_from_R["Z2"];
        Z1=origZ1;
        Z2=origZ2;
        double tol1 = TOL_P/Z1.n_rows;
        double tol2 = TOL_P/Z2.n_rows;
        boundaries(Z1,tol1,1-tol1);
        boundaries(Z2,tol2,1-tol2);
        Z1 /= repmat( sum(Z1,1), 1, Z1.n_cols );
        Z2 /= repmat( sum(Z2,1), 1, Z2.n_cols );
        alpha1 = sum(Z1,0) / Z1.n_rows;
        alpha2 = sum(Z2,0) / Z2.n_rows;

    }

    LBM & operator=(const LBM & orig)
    {
        Z1=orig.Z1;
        Z2=orig.Z2;
        alpha1=orig.alpha1;
        alpha2=orig.alpha2;

        return *this;
    }

    inline
    double entropy() // TODO verif VB
    {
        return accu( Z1 % log(Z1) ) + accu(Z2 % log(Z2));
    }

    template<class model_type, class network_type>
    inline
    void e_step(model_type & model, network_type & net)
    {
        double tol1 = TOL_P/Z1.n_rows;
        double tol2 = TOL_P/Z2.n_rows;
        double step_size;
        unsigned int niter=0;

        do
        {
            mat lZ1 = repmat(log(alpha1),Z1.n_rows,1);
            mat lZ2 = repmat(log(alpha2),Z2.n_rows,1);

            e_fixed_step(*this, model, net, lZ1, lZ2);

            lZ1 -= repmat( mean(lZ1,1), 1, lZ1.n_cols );
            lZ2 -= repmat( mean(lZ2,1), 1, lZ2.n_cols );

            lZ1 -= repmat( max(lZ1,1), 1, lZ1.n_cols );
            lZ2 -= repmat( max(lZ2,1), 1, lZ2.n_cols );

            lZ1 = exp(lZ1);
            lZ2 = exp(lZ2);

            lZ1 /= repmat( sum(lZ1,1), 1, lZ1.n_cols );
            lZ2 /= repmat( sum(lZ2,1), 1, lZ2.n_cols );

            boundaries(lZ1,tol1,1-tol1);
            boundaries(lZ2,tol2,1-tol2);
            lZ1 /= repmat( sum(lZ1,1), 1, lZ1.n_cols );
            lZ2 /= repmat( sum(lZ2,1), 1, lZ2.n_cols );

            double step_size1 = max(max(abs( Z1-lZ1 )));
            double step_size2 = max(max(abs( Z2-lZ2 )));

            step_size = (step_size1>step_size2) ? step_size1 : step_size2;
            
            #ifdef DEBUG_E
                fprintf(stderr,"E iteration: %i %f [%f %f]\n",niter,step_size,step_size1,step_size2);
            #endif

            niter++;

            Z1 = lZ1;
            Z2 = lZ2;
        } while( step_size>TOL_F && niter<10);
    }

    inline
    double m_step()
    {
        alpha1 = sum(Z1,0) / Z1.n_rows;
        alpha2 = sum(Z2,0) / Z2.n_rows;

        return accu( Z1*log(alpha1).t() ) + accu( Z2*log(alpha2).t() );
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["Z1"] = Z1;
        values["alpha1"] = alpha1;
        values["Z2"] = Z2;
        values["alpha2"] = alpha2;

        return values;
    }



        

};


