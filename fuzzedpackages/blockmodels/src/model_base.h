
template<class model_type, class network_type>
inline
vec grad_logf(model_type & model, network_type & net, unsigned int i, unsigned int j, unsigned int q, unsigned int l)
{
    vec g(model.n_parameters);
    for(unsigned int param=0; param<model.n_parameters; param++)
        g(param)=grad_logf(model, net, i,j,q,l,param);
    return g;
}
              
template<class model_type, class network_type>
inline
double PL(model_type & model, SBM & membership, network_type & net)
{
    double S=0;

    for(unsigned int i=0; i<membership.Z.n_rows; i++)
        for(unsigned int j=0; j<membership.Z.n_rows; j++)
            if(i!=j)
                for(unsigned int q=0; q<membership.Z.n_cols; q++)
                    for(unsigned int l=0; l<membership.Z.n_cols; l++)
                        S += logf(model,net,i,j,q,l) * membership.Z(i,q) * membership.Z(j,l);
    return S;
}

template<class model_type, class network_type>
inline
double PL(model_type & model, SBM_sym & membership, network_type & net)
{
    double S=0;

    for(unsigned int i=0; i<membership.Z.n_rows; i++)
        for(unsigned int j=i+1; j<membership.Z.n_rows; j++)
            for(unsigned int q=0; q<membership.Z.n_cols; q++)
                for(unsigned int l=0; l<membership.Z.n_cols; l++)
                    S += logf(model,net,i,j,q,l) * membership.Z(i,q) * membership.Z(j,l);
    return S;
}

template<class model_type, class network_type>
inline
double PL(model_type & model, LBM & membership, network_type & net)
{
    double S=0;

    for(unsigned int i=0; i<membership.Z1.n_rows; i++)
        for(unsigned int j=0; j<membership.Z2.n_rows; j++)
            if(i!=j)
                for(unsigned int q=0; q<membership.Z1.n_cols; q++)
                    for(unsigned int l=0; l<membership.Z2.n_cols; l++)
                        S += logf(model,net,i,j,q,l) * membership.Z1(i,q) * membership.Z2(j,l);
    return S;
}

template<class model_type, class network_type>
inline
vec grad(model_type & model, SBM & membership, network_type & net)
{
    vec G=zeros<vec>(model.n_parameters);

    for(unsigned int i=0; i<membership.Z.n_rows; i++)
        for(unsigned int j=0; j<membership.Z.n_rows; j++)
            if(i!=j)
                for(unsigned int q=0; q<membership.Z.n_cols; q++)
                    for(unsigned int l=0; l<membership.Z.n_cols; l++)
                        G += grad_logf(model,net,i,j,q,l) * membership.Z(i,q) * membership.Z(j,l);
    return G;
}

template<class model_type, class network_type>
inline
vec grad(model_type & model, SBM_sym & membership, network_type & net)
{
    vec G=zeros<vec>(model.n_parameters);

    for(unsigned int i=0; i<membership.Z.n_rows; i++)
        for(unsigned int j=i+1; j<membership.Z.n_rows; j++)
            for(unsigned int q=0; q<membership.Z.n_cols; q++)
                for(unsigned int l=0; l<membership.Z.n_cols; l++)
                        G += grad_logf(model,net,i,j,q,l) * membership.Z(i,q) * membership.Z(j,l);
    return G;
}

template<class model_type, class network_type>
inline
vec grad(model_type & model, LBM & membership, network_type & net)
{
    vec G=zeros<vec>(model.n_parameters);

    for(unsigned int i=0; i<membership.Z1.n_rows; i++)
        for(unsigned int j=0; j<membership.Z2.n_rows; j++)
            if(i!=j)
                for(unsigned int q=0; q<membership.Z1.n_cols; q++)
                    for(unsigned int l=0; l<membership.Z2.n_cols; l++)
                        G += grad_logf(model,net,i,j,q,l) * membership.Z1(i,q) * membership.Z2(j,l);
    return G;
}

template<class membership_type, class model_type>
inline
model_type copy_and_add(model_type & model, membership_type & membership, vec toadd)
{
    return model_type(membership,model.to_vector()+toadd);
}

template<class model_type>
inline
void print_model(model_type & model)
{
}
