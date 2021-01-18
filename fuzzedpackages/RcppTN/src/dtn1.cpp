// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

# include <Rcpp.h>

double dtn1(const double x,
            const double mean,
            const double sd,
            const double low,
            const double high
            ) {
    double zx = (x - mean) / sd ;
    double za = (low - mean) / sd ;
    double zb = (high - mean) / sd ;
    // Default to 0 for out of interval densitites
    double dens = 0.0 ;
    
    if ((zx >= za) & (zx <= zb)) {
        // If interval, calculate tru density        
        double q2 = ::Rf_pnorm5(zb, 0.0, 1.0, 1, 0) - ::Rf_pnorm5(za, 0.0, 1.0, 1, 0) ;
        double q1 = ::Rf_dnorm4(zx, 0.0, 1.0, 0) / sd ;
        dens = q1 / q2 ;
    }
    
    return(dens) ;
}
