
template<class model_type, class network_type>
inline
void e_fixed_step(SBM & membership, model_type & model, network_type & net, mat & lZ)
{
    
    for(unsigned int i=0; i<lZ.n_rows; i++)
        for(unsigned int j=0; j<lZ.n_rows; j++)
            if(i!=j)
                for(unsigned int q=0; q<lZ.n_cols; q++)
                    for(unsigned int l=0; l<lZ.n_cols; l++ )
                        lZ(i,q) += membership.Z(j,l) * (
                                                         logf(model,net,i,j,q,l)
                                                         +
                                                         logf(model,net,j,i,l,q)
                                                       );
}

template<class model_type, class network_type>
inline
void e_fixed_step(SBM_sym & membership, model_type & model, network_type & net, mat & lZ)
{
    
        for(unsigned int i=0; i<lZ.n_rows; i++)
            for(unsigned int j=0; j<lZ.n_rows; j++)
                if(i!=j)
                    for(unsigned int q=0; q<lZ.n_cols; q++)
                        for(unsigned int l=0; l<lZ.n_cols; l++ )
                            lZ(i,q) += membership.Z(j,l) * logf(model,net,i,j,q,l);
}

template<class model_type, class network_type>
inline
void e_fixed_step(LBM & membership, model_type & model, network_type &net, mat & lZ1, mat & lZ2)
{
    for(unsigned int i=0; i<lZ1.n_rows; i++)
        for(unsigned int j=0; j<lZ2.n_rows; j++)
            for(unsigned int q=0; q<lZ1.n_cols; q++)
                for(unsigned int l=0; l<lZ2.n_cols; l++ )
                {
                    double lfijql =logf(model,net,i,j,q,l);
                    lZ1(i,q) += membership.Z2(j,l) * lfijql;
                    lZ2(j,l) += membership.Z1(i,q) * lfijql;
                }
}

template<class membership_type, class model_type, class network_type>
inline
double m_step(membership_type & membership, model_type & model, network_type & net)
{

    mat pseudo_hessian(model.n_parameters,model.n_parameters);
    pseudo_hessian.eye();

    /* Be carreful
     * The gradient is the gradient of the opposite liklihood
     * The hessian is the hessian of the opposite liklihood
     */

    vec gradient = - grad(model,membership, net);
    double value = - PL(model,membership, net);
    #ifdef DEBUG_M
    fprintf(stderr,"M numerical optim: %f\n",value);
    gradient.print("gradient:");

    #endif

    double gain=0;
    unsigned int bfgs_iteration_counter=0;
    do
    {
        bfgs_iteration_counter++;
        vec direction = solve(pseudo_hessian, -gradient);
        
        #ifdef DEBUG_M
        print_model(model);
        pseudo_hessian.print("pseudo_hessian:");
        #endif
        
        double D_direction = mat(direction.t()*gradient)(0);
        
        if(D_direction>=0)
            break;

        #ifdef DEBUG_M
        direction.print("direction:");
        #endif

        double a=1.99 * maximum_step_in_direction(membership,model,net,direction);
        double new_value=0;

        unsigned int line_search_iteration_counter = 0;
        do
        {
            line_search_iteration_counter++;
            a *= .5;
            model_type model_added = copy_and_add(model,membership,a*direction);
            new_value = - PL(model_added,membership, net);
            #ifdef DEBUG_M
            fprintf(stderr,"line search: a=%f\tvalue=%f\told_value=%f\tobj=%f\n",a,new_value,value,value+.25*a*D_direction);
            #endif
        } while(new_value-value > .25 * a * D_direction && line_search_iteration_counter < LINE_SEARCH_ITER_MAX);

        if(line_search_iteration_counter == LINE_SEARCH_ITER_MAX)
            break;

        model = copy_and_add(model,membership,a*direction);

        vec new_gradient = - grad(model,membership, net);
        vec diff_grad = new_gradient-gradient;

        pseudo_hessian += diff_grad*diff_grad.t()/(a*mat(diff_grad.t()*direction)(0))
                          - (pseudo_hessian*direction*direction.t()*pseudo_hessian)/mat(direction.t()*pseudo_hessian*direction)(0);

        gain = value-new_value;
        gradient = new_gradient;
        value = new_value;
        #ifdef DEBUG_M
        fprintf(stderr,"M numerical optim: %f\n",value);
        gradient.print("gradient:");
        #endif

    } while( (gain > TOL_M || bfgs_iteration_counter < BFGS_ITER_MIN) && bfgs_iteration_counter < BFGS_ITER_MAX);

    return -value; /* value is the opposite of the pseudo likelihood */
}


