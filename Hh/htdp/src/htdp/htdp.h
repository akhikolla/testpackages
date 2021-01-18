#ifndef HTDP_H
#define HTDP_H

#ifdef __cplusplus
extern "C" {
#endif

// SUBROUTINE MODEL
// Obtain parameters defining crustal motion model
void model_();

// SUBROUTINE SETTP
// Initialize transformation parameters from ITRF94 to other reference frames
void settp_();

// SUBROUTINE SETRF
// Initialize conversion table between reference frame identifiers
void setrf_();

// SUBROUTINE IYMDMJ
// Convert date to modified Julian date
void iymdmj_(const int *iyr, const int *imon, const int *iday, int *mjd);

// SUBROUTINE PREDV
// Predict velocity in iopt reference frame
void predv_(const double *ylat, const double *ylon, const double *eht,
            const double *date, const int *iopt, int *jregn,
            double *vn, double *ve, double *vu);

// SUBROUTINE NEWCOR
// Predict coordinates and displacements from time MIN1 to time MIN2
void newcor_(const double *ylat, const double *ylon, const double *htold,
             const int *min1, const int *min2, double *ylat3, double *ylon3,
             double *htnew, double *dn, double *de, double *du,
             const double *vn, const double *ve, const double *vu);

#ifdef __cplusplus
}
#endif

#endif // HTDP_H
