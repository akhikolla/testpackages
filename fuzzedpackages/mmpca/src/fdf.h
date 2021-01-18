void f_vxi(double * v, const double * x, const int p, const int k);
double f_obj(const double * theta, const double ** x, const double ** masks,
  const double * lambda, const int k, const int * inds, const int* p,
  const int m, const int n, const int len);
void d_obj(double * grad, const double * theta, const double ** x,
  const double ** masks, const double * lambda, const int k, const int * inds,
  const int* p, const int m, const int n, const int len, const double * indices,
  const int indices_len, const int num_threads);
void inv_v(double * xi, double * c_t, int n);
void init_parallel();
