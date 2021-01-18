#ifndef MATH_UTILS_H
#define MATH_UTILS_H
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define max_length 1000


inline int suminline(const std::vector<int> &data) {
	return accumulate(data.begin(), data.end(), 0);
}

inline double suminline(const std::vector<double> &data) {
	return accumulate(data.begin(), data.end(), 0.0);
}

inline double meaninline(const std::vector<double> &data) {
	return suminline(data) / (int)data.size();
}

inline double varinline(const std::vector<double> &data) {
	double sum1 = 0, sum2 = 0;
	for (int i = 0; i < (int)data.size(); i++) {
		sum1 += data[i];
		sum2 += data[i] * data[i];
	}
	return (sum2 - sum1 * sum1 / (int)data.size()) / (int)data.size();
}

inline double medianinline(const std::vector<double> &data){
    int N = (int) data.size();
    if (N%2 == 0){
        return((data[N / 2 - 1] + data[N / 2]) / 2.0);
    }
    else{
        return(data[(N - 1) / 2]);
    }
}

inline double quantileinline(const std::vector<double> &data, double quantile) {
	int index = (int)((data.size() - 1) * quantile);
	return data[index];
}


inline double ldnorminline(double x) {
  return -x*x/2.0 - 0.5 * log(2.0 * 3.1415926);
}

inline double ldnorminline(double x, double mu, double sigma = 1) {
  return ldnorminline((x-mu)/sigma) - log(sigma);
}

inline double ldtinline(double x, double df){
  return -0.5 * (1.0 + df) * log(x * x / df + 1.0) + lgamma(0.5 * (df + 1.0)) - lgamma(df / 2.0) - 0.5 * log(df * 3.1415926);
}

inline double ldscaleinvchisq(double sigma2, double nu, double s2){
    return log(s2 * nu / 2.0) * nu / 2.0 - std::lgamma(nu / 2.0) - s2 * nu / (2.0 * sigma2) - log(sigma2) * (nu / 2.0 + 1.0);
}

inline double ldnormmixinline(double y, std::vector<double> pi, std::vector<double> mu, std::vector<double> sigma2){
    int i, K = pi.size();
    std::vector<double> dens;
    dens.resize(K);
    double result = 0.0, maxdens = ldnorminline(y, mu[0], sqrt(sigma2[0]));
    dens[0] = maxdens;
    for (i = 1; i < K; i ++){
        dens[i] = ldnorminline(y, mu[i], sqrt(sigma2[i]));
        maxdens = MAX(maxdens, dens[i]);
    }
    for (i = 0; i < K; i ++){
        result += pi[i] * exp(dens[i] - maxdens);
    }
    result = log(result) + maxdens;
    return(result);
}

inline double runiforminline(){
  double x(1);
  x = R::runif(0,1);
  return(x);
}

inline double rnorminline() {
  double x(1);
  x = R::rnorm(0,1);
  return(x);
}

inline int randintinline(int low, int high){

    int temp = int((runiforminline() * (high - low + 1))) + low;
    if (temp > high) temp = high;
    if (temp < low) temp = low;
    return temp;
}

inline double randomchisq (double df){
  std::default_random_engine generator(runiforminline());
  std::chi_squared_distribution<double> distribution (df);
  return(distribution(generator));
}

inline double rtinline (const double df){
  double x = rnorminline(), y = randomchisq(df);
  return(x / sqrt(y / df));
}

inline double randomgamma (double df, double scale){
  std::default_random_engine generator(runiforminline());
  std::gamma_distribution<double> distribution (df, scale);
  return(distribution(generator));
}

inline std::vector<double> randomdirichlet (std::vector<double> alpha){
    int K = alpha.size();
    std::vector<double> rvtemp;
    if (K < 2){
        Rprintf("wrong class number for dirichlet rv.\n");
        stop(0);
    }
    else{
        int j;
        double cumsumtemp = 0.0;
        rvtemp.resize(K);
        for (j = 0; j < K; j ++){
            rvtemp[j] = randomgamma(alpha[j], 1.0);
            cumsumtemp += rvtemp[j];
        }
        for (j = 0; j < K; j ++){
            rvtemp[j] /= cumsumtemp;
        }
    }
    return(rvtemp);
}

std::string convertInt(int number)
{
    std::stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}

std::string convertDouble(double number)
{
    std::stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}

void convertNumericMat (Rcpp::NumericMatrix mat, std::vector<std::vector<double> > &vmat){

    int i, j, M = mat.nrow(), N = mat.ncol();
    vmat.resize(M);
    for(i = 0; i < M; i ++){
        vmat[i].resize(N);
        for(j = 0; j < N; j ++)
            vmat[i][j] = mat(i, j);
    }
    return;
}




#endif // MATH_UTILS_H







