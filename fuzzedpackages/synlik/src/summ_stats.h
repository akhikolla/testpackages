/* The data -> statistics routines */
void slacf(double *acf,double *x,int *n,int *n_reps,int *max_lag,double *NAcode,int *correlation);
void slnlar(double *beta, double *x,int *n,int *n_reps,int *n_terms,
            int *lag,int *power,double *NAcode);
void order_reg(double *beta, double *x,double *z,int *n,int *n_reps,int *np,int *diff);
void blowC(double *n,double *theta,double *e,double *e1,int *burn_in,int *n_t, int *n_reps);

/* LAPACK based matrix routines, from mgcv's mat.c */

void mgcv_qr(double *x, int *r, int *c,int *pivot,double *tau);
void mgcv_qrqy(double *b,double *a,double *tau,int *r,int *c,int *k,int *left,int *tp);
void mgcv_backsolve(double *R,int *r,int *c,double *B,double *C, int *bc);
