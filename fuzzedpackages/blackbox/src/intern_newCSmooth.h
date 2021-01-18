#ifndef H_INTERN_NEWCSMOOTH
#define H_INTERN_NEWCSMOOTH

bool intern_newCSmooth( //header declaration
    double *xy,
    int *nrowxy, //with replicates
    int *ncol,
    int *nuniquerows, // required to tell the allocated size of c, d, D, u in *R*
    // double *maxSmoothness,
    //    double *c,
    //    double *d,
    //    double *D,
    //    double *u,
    //    double *lambda,
    double *GCV,
    //    double *covfnparam, //p+1 form
    //int *fneval,
    int *optimiseBool,
    int *verbosity
    //    double *initcovfnparam //p+1 form  // new arg 2015/10/24
);

#endif

