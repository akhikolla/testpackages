

// struct for returning results of EM
template<class membership_type, class model_type>
struct result
{
    membership_type membership;
    model_type model;
    double PL;
    double H;
    
    template<class network_type>
    result(membership_type membership_init, network_type net) : 
        membership(membership_init),
        model(membership_init, net)
    {}

    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["membership"] = membership.export_to_R();
        values["model"] = model.export_to_R();
        values["PL"] = PL;
        values["H"] = H;

        return values;
    }
};

// em
template<class membership_type, class model_type, class network_type, bool real_EM>
inline
result<membership_type,model_type> em(
            membership_type & membership_init,
            network_type & net
        )
{
    result<membership_type,model_type> res(membership_init,net);

    res.H = res.membership.entropy();
    res.PL = res.membership.m_step();
    res.PL += m_step(res.membership, res.model, net);

    #ifdef DEBUG_EM
    std::cout << "M" << std::endl;
    fprintf(stderr,"J = %f\t H = %f\t PL = %f\n",res.H+res.PL,res.H,res.PL);
    #endif
    
    if(real_EM)
    {

        double J=res.H+res.PL;
        double old_J=0;
        
        do
        {
            res.membership.e_step (res.model, net);
            res.H = res.membership.entropy();
            
            #ifdef DEBUG_EM
            std::cout << "E" << std::endl;
            #endif
            
            res.PL = res.membership.m_step();
            res.PL += m_step<membership_type,model_type> (res.membership, res.model, net);
            
            old_J=J;
            
            #ifdef DEBUG_EM
            std::cout << "M" << std::endl;
            fprintf(stderr,"J = %f\t H = %f\t PL = %f\n",res.H+res.PL,res.H,res.PL);
            #endif

            J=res.H+res.PL;

        } while(J-old_J > TOL_EM);
        
        #ifdef DEBUG_EM
        fprintf(stderr,"EM finised\n");
        #endif
    }


    return res;
}

