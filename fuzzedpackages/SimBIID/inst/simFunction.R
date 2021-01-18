

    // initialise variables
    double tstar = 0.0, u_tmp = 0.0, totrate = 0.0, cumrate = 0.0;
    int i = 0;
    
    // initialise time and rates of the system
    double t = tstart;
    RATELINES1 
    
    // check states
    for(i = 0; i < u.size(); i++){
        if(u[i] < 0) {
            stop("Some states less than zero");
        }
    }
    
    // check rates
    for(i = 0; i < rates.size(); i++){
        if(R_FINITE(rates[i])) {
            if(rates[i] < 0.0) {
                stop("Some rates less than zero or non-finite");
            }
        } else {
            stop("Some rates are non-finite");
        }
    }

    // set up output vector
    MATCHCRIT0
    TSPAN0
    
    // sample next event time
    if(totrate > 0.0) {
        tstar = t + R::rexp(1.0 / totrate);
        AFTER_TSTAR0
        while(tstar < tstop){
            TSPAN1
            // update event type
            u_tmp = R::runif(0.0, totrate);
            RATELINES2  
            
            // check states
            for(i = 0; i < u.size(); i++){
                if(u[i] < 0) {
                    stop("Some states less than zero");
                }
            }
            
            // check rates
            for(i = 0; i < rates.size(); i++){
                if(R_FINITE(rates[i])) {
                    if(rates[i] < 0.0) {
                        stop("Some rates less than zero or non-finite");
                    }
                } else {
                    stop("Some rates are non-finite");
                }
            }
            
            // update time
            t = tstar;
            
            // sample next event time
            if(totrate > 0.0) {
                tstar = t + R::rexp(1.0 / totrate);
            } else {
                tstar = tstop;
            }
            AFTER_TSTAR1
            RATELINES3
        }
    }
    
    // record final event time
    if(totrate > 0.0) {
        t = tstop;
    }
    
    MATCHCRIT1
    
    // return output
    RETURNCRIT
}')
