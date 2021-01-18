// utils.cpp

#include "dmbc.h"

#define both_non_NA(a,b) (!ISNAN(a) && !ISNAN(b))

void logit(double* res, const double* p, int n){
  for(int i = 0; i < n; i++){
    res[i] = log((double) p[i]/(1 - p[i]));
  }
}

void expit(double* res, const double* x, int n){
  for(int i = 0; i < n; i++){
    res[i] = (double) 1/(1 + exp(-x[i]));
  }
}

void exp_vec(double* res, const double* x, int n){
  for(int i = 0; i < n; i++){
    res[i] = (double) exp(x[i]);
  }
}

// Euclidean distance
double euclidean(const double *x, int nr, int nc, int i1, int i2){
  double dev = 0, dist = 0;
  int count = 0;

  for(int j = 0; j < nc; j++){
    if(both_non_NA(x[i1], x[i2])){
      dev = (x[i1] - x[i2]);
      if(!ISNAN(dev)){
        dist += dev * dev;
        count++;
      }
    }
    i1 += nr;
    i2 += nr;
  }
  if(count == 0){
    return NA_REAL;
  }
  if(count != nc){
    dist /= ((double)count/nc);
  }
  return sqrt(dist);
}

// Euclidean distance
void dist(double* d, const double* x, int nr, int nc){
  size_t count = 0;
  for(int j = 0; j <= nr; j++){
    for(int i = (j + 1); i < nr; i++){
      d[count++] = euclidean(x, nr, nc, i, j);
    }
  }
}

bool any_na_nan(const arma::vec x, const int& n){
	bool res = false;

	for(int i = 0; i < n; i++){
		if(ISNAN(x(i))){
			res = true;
			break;
		}
	}

	return res;
}

void sample_no_rep(int n, double* p, int* perm, int nans, int* ans){
	double rT, mass, totalmass;
	int i, j, k, n1;

	for(i = 0; i < n; i++){
		perm[i] = i + 1;
	}
	revsort(p, perm, n);

	totalmass = 1;
	for(i = 0, n1 = (n - 1); i < nans; i++, n1--){
		rT = totalmass * unif_rand();
		mass = 0;
		for(j = 0; j < n1; j++){
			mass += p[j];
			if(rT <= mass){
				break;
			}
		}
		ans[i] = perm[j];
		totalmass -= p[j];
		for(k = j; k < n1; k++){
			p[k] = p[k + 1];
			perm[k] = perm[k + 1];
		}
	}
}

void tableC(int* counts, const int* x, int nelem, int ndistelem){
  int* xtmp = new int[nelem];
  for(int i = 0; i < nelem; i++){
    xtmp[i] = x[i];
  }
  // R_isort(xtmp, nelem);
  R_qsort_int(xtmp, 1, nelem);
  
  int xlead = 1;
  counts[xlead - 1] = 0;
  for(int i = 0; i < nelem; i++){
    if(xtmp[i] != xlead){
      xlead++;
      counts[xlead - 1] = 0;
    }
    counts[xlead - 1]++;
  }

  delete[] xtmp;
}

int factorial(const int& x){
  return R::gammafn(x + 1);
}

void permutations(int* perm, int n, int nperm, int byrow){
	int* v = new int[n];
	for(int i = 0; i < n; i++){
		v[i] = (i + 1);
	}

	int p = 0;
	do{
		if(byrow == 0){
			for(int i = 0; i < n; i++){
				perm[p + nperm*i] = v[i];
			}
		}
		else{
			for(int i = 0; i < n; i++){
				perm[i + n*p] = v[i];
			}
		}
		p++;
	} while (std::next_permutation(v, v + n));  // 'v' is the address of the array's first element, while 'v + n' that of its n-th element

	delete[] v;
}

void which_min(int* ans, const double* r, int n){
	double s = R_PosInf;
	int indx = NA_INTEGER;

	for(int i = 0; i < n; i++){
		// if(!ISNAN(r[i]) && (r[i] < s || indx == NA_INTEGER)){
		if(r[i] < s || indx == NA_INTEGER){
			s = r[i];
			indx = i;
		}
	}

	if(indx != NA_INTEGER){
		*ans = indx + 1;
	}
	else{
		*ans = 0;	
	}
}
