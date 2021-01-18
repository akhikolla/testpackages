// -*- mode: C++; c-indent-level: 2; c-basic-offset: 2; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

template <typename Real>
Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> tsolveAndMultiply(
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>& A,
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>& C) {
  Eigen::TriangularView<
      const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Upper>
      T = A.template triangularView<Eigen::Upper>();
  return T.template solve<Eigen::OnTheRight>(C);  // C %*% solve(A)
}

template <typename Real>
std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> QRdecomp(
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>&
        M) {  // for nrows >= ncols
  Eigen::HouseholderQR<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> qr =
      M.householderQr();
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> R_ =
      qr.matrixQR().template triangularView<Eigen::Upper>();
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Q_ = qr.householderQ();
  const size_t m = M.rows();
  const size_t n = M.cols();
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> R = R_.block(0, 0, n, n);
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Q = Q_.block(0, 0, m, n);
  return {Q, R};
}

template <typename Real>
Eigen::Matrix<Real, Eigen::Dynamic, 1> solve(
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>& A,
    const Eigen::Matrix<Real, Eigen::Dynamic, 1>& b) {
  return A.colPivHouseholderQr().solve(b);
}

template <typename Real>
int rankMatrix(const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>& M) {
  const Eigen::ColPivHouseholderQR<
      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>>
      qr(M);
  return qr.rank();
}

std::normal_distribution<double> gaussian(0.0, 1.0);

template <typename Real>
Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> gmatrix(
    const size_t nrows,
    const size_t ncols,
    std::default_random_engine& generator) {
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> G(nrows, ncols);
  for(size_t i = 0; i < nrows; i++) {
    for(size_t j = 0; j < ncols; j++) {
      G(i, j) = (Real)(gaussian(generator));
    }
  }
  return G;
}

template <typename Real>
const int sgn(const Real x) {
  return (0 < x) - (x < 0);
}

std::uniform_real_distribution<double> runif(0.0, 1.0);

const size_t spow(size_t base, size_t exp) {
  size_t result = 1;
  while(exp) {
    if(exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }
  return result;
}

template <typename Real>
const Real approx(const Real x, const size_t n) {
  return round(x * spow(10, n)) / spow(10, n);
}

template <typename Real>
const std::vector<Real> fidSample(
    const Eigen::Matrix<Real, Eigen::Dynamic, 1>& VT2,
    const Eigen::Matrix<Real, Eigen::Dynamic, 1>& VTsum,
    const Real L,
    const Real U,
    std::default_random_engine& generator) {
  Real ZZ, wt;                    // outputs
  const size_t p = VTsum.size();  // = VT2.size()
  std::vector<size_t> high, low, zero, zeronot;
  size_t lhigh = 0, llow = 0, lzero = 0;
  for(size_t i = 0; i < p; i++) {
    if(VT2.coeff(i) == 0) {
      lzero += 1;
      zero.push_back(i);
    } else {
      zeronot.push_back(i);
      if(VT2.coeff(i) > 0) {
        lhigh += 1;
        high.push_back(i);
      } else {
        llow += 1;
        low.push_back(i);
      }
    }
  }
  const size_t lzeronot = p - lzero;
  Real MAX, MIN;
  const Real infty = std::numeric_limits<Real>::infinity();
  if((lhigh > 0 && llow > 0) || lzero > 0) {
    Real temp;
    std::vector<int> UU(p);
    std::vector<int> LL(p);
    std::vector<int> SS(p);
    for(size_t i = 0; i < p; i++) {
      UU[i] = sgn<Real>(U - VTsum.coeff(i));
      LL[i] = sgn<Real>(L - VTsum.coeff(i));
      SS[i] = sgn<Real>(VT2.coeff(i));
    }
    if(lzero == p) {
      MAX = infty;
      MIN = -infty;
      bool anyUUpos = false, anyLLneg = false;
      size_t i = 0;
      while(i < p && !anyUUpos) {
        if(UU[i] > 0) {
          anyUUpos = true;
        }
        i += 1;
      }
      i = 0;
      while(i < p && !anyLLneg) {
        if(LL[i] < 0) {
          anyLLneg = true;
        }
        i += 1;
      }
      temp = (Real)(anyUUpos && anyLLneg);
    } else {
      size_t i = 0;
      bool c1 = true, c2 = true, d1 = true, d2 = true;
      while(i < lzero && c1) {
        if(UU[zero[i]] == -1 || LL[zero[i]] == -1) {
          c1 = false;
        }
        i += 1;
      }
      i = 0;
      while(i < lzero && c2) {
        if(UU[zero[i]] == 1 || LL[zero[i]] == 1) {
          c2 = false;
        }
        i += 1;
      }
      i = 0;
      while(i < p && d1) {
        if(SS[i] == -1) {
          d1 = false;
        }
        i += 1;
      }
      i = 0;
      while(i < p && d2) {
        if(SS[i] == 1) {
          d2 = false;
        }
        i += 1;
      }
      if((d1 && c1) || (d2 && c2)) {
        MAX = infty;
        MIN = infty;
        for(size_t i = 0; i < lzeronot; i++) {
          const size_t zni = zeronot[i];
          MIN =
              std::min(MIN, std::min((U - VTsum.coeff(zni)) / VT2.coeff(zni),
                                     (L - VTsum.coeff(zni)) / VT2.coeff(zni)));
        }
        temp = 1 - (atan(MIN) / PI + 0.5);
      } else if((d2 && c1) || (d1 && c2)) {
        MIN = -infty;
        MAX = -infty;
        for(size_t i = 0; i < lzeronot; i++) {
          const size_t zni = zeronot[i];
          MAX =
              std::max(MAX, std::max((U - VTsum.coeff(zni)) / VT2.coeff(zni),
                                     (L - VTsum.coeff(zni)) / VT2.coeff(zni)));
        }
        temp = atan(MAX) / PI + 0.5;
      } else {
        Real Hmax = -infty;
        Real Hmin = infty;
        for(size_t i = 0; i < lhigh; i++) {
          const size_t hi = high[i];
          Real xu = (U - VTsum.coeff(hi)) / VT2.coeff(hi);
          Real xl = (L - VTsum.coeff(hi)) / VT2.coeff(hi);
          Hmax = std::max(Hmax, std::max(xu, xl));
          Hmin = std::min(Hmin, std::min(xu, xl));
        }
        Real Lmax = -infty;
        Real Lmin = infty;
        for(size_t i = 0; i < llow; i++) {
          const size_t li = low[i];
          Real xu = (U - VTsum.coeff(li)) / VT2.coeff(li);
          Real xl = (L - VTsum.coeff(li)) / VT2.coeff(li);
          Lmax = std::max(Lmax, std::max(xu, xl));
          Lmin = std::min(Lmin, std::min(xu, xl));
        }
        Real bpos, tpos, bneg, tneg;
        if(approx<Real>(Lmin - Hmax, 12) >= 0) {
          bpos = -infty;
          tpos = Hmax;
          bneg = Lmin;
          tneg = infty;
        } else if(approx<Real>(Hmin - Lmax, 12) >= 0) {
          bpos = Hmin;
          tpos = infty;
          bneg = -infty;
          tneg = Lmax;
        } else {
          bpos = -infty;
          tpos = infty;
          bneg = -infty;
          tneg = infty;
        }
        Real Pprob, Nprob;
        if(tpos == infty) {
          Pprob = 1 - (atan(bpos) / PI + 0.5);
        } else {
          Pprob = atan(tpos) / PI + 0.5;
        }
        if(tneg == infty) {
          Nprob = 1 - (atan(bneg) / PI + 0.5);
        } else {
          Nprob = atan(tneg) / PI + 0.5;
        }
        temp = Pprob + Nprob;
        Pprob = Pprob / temp;
        Nprob = 1 - Pprob;
        if(runif(generator) <= Nprob) {
          MIN = bneg;
          MAX = tneg;
        } else {
          MIN = bpos;
          MAX = tpos;
        }
      }
    }
    const Real y = atan(MAX) / PI + 0.5;
    const Real x = atan(MIN) / PI + 0.5;
    const Real u = x + (y - x) * runif(generator);
    ZZ = tan(PI * (u - 0.5));
    const Real ZZ2 = ZZ * ZZ;
    wt = exp(-ZZ2 / 2) * (1 + ZZ2) * temp;
  } else {
    MAX = -infty;
    MIN = infty;
    for(size_t i = 0; i < p; i++) {
      const Real xu = (U - VTsum.coeff(i)) / VT2.coeff(i);
      const Real xl = (L - VTsum.coeff(i)) / VT2.coeff(i);
      MAX = std::max(MAX, std::max(xu, xl));
      MIN = std::min(MIN, std::min(xu, xl));
    }
    const Real y = atan(MAX) / PI + 0.5;
    const Real x = atan(MIN) / PI + 0.5;
    const Real u = x + (y - x) * runif(generator);
    ZZ = tan(PI * (u - 0.5));
    const Real ZZ2 = ZZ * ZZ;
    wt = exp(-ZZ2 / 2) * (1 + ZZ2) * (y - x);
  }

  const std::vector<Real> out = {ZZ, wt};
  return out;
}

const Eigen::VectorXi cppunique(Eigen::VectorXi& v) {
  int size = v.size();
  for(int i = 0; i < size; i++) {
    for(int j = i + 1; j < size;) {
      if(v.coeff(j) == v.coeff(i)) {
        for(int k = j; k + 1 < size; k++) {
          v(k) = v.coeff(k + 1);
        }
        size--;
      } else {
        j++;
      }
    }
  }
  return v.head(size);
}

const std::vector<std::vector<int>> cartesianProduct(
    const std::vector<std::vector<int>>& v) {
  std::vector<std::vector<int>> s = {{}};
  for(auto& u : v) {
    std::vector<std::vector<int>> r;
    for(auto& x : s) {
      for(auto y : u) {
        r.push_back(x);
        r.back().push_back(y);
      }
    }
    s.swap(r);
  }
  return s;
}

const std::vector<std::vector<int>> combinations(const std::vector<int>& C,
                                                 const int n) {
  std::vector<std::vector<int>> sets;
  for(size_t i = 0; i < C.size(); i++) {
    sets.push_back({C[i], C[i] + n});
  }
  return cartesianProduct(sets);
}

const Eigen::MatrixXi vv2matrix(std::vector<std::vector<int>> U,
                                const size_t nrow,
                                const size_t ncol) {
  Eigen::MatrixXi out(nrow, ncol);
  for(size_t i = 0; i < nrow; i++) {
    for(size_t j = 0; j < ncol; j++) {
      out(i, j) = U[j][i];
    }
  }
  return out;
}

const Eigen::VectorXi Vsort(Eigen::VectorXi& V) {
  std::sort(V.data(), V.data() + V.size());
  return V;
}

template <typename Real>
const std::vector<Real> Vcumsum(
    const Eigen::Matrix<Real, Eigen::Dynamic, 1>& vec) {
  std::vector<Real> out(vec.size());
  out[0] = vec.coeff(0);
  for(int i = 1; i < vec.size(); i++) {
    out[i] = out[i - 1] + vec.coeff(i);
  }
  return out;
}

template <class T>
std::vector<T> zero2n(T n) {
  std::vector<T> out(n);
  for(T i = 0; i < n; i++) {
    out[i] = i;
  }
  return out;
}

const std::vector<size_t> sample_int(const size_t n,
                                     std::default_random_engine& generator) {
  std::vector<size_t> elems = zero2n(n);
  std::shuffle(elems.begin(), elems.end(), generator);
  return elems;
}

template <typename Real>
const Real rchisq(int df, std::default_random_engine& generator) {
  Real out = 0;
  for(int i = 0; i < df; i++) {
    const Real toadd = (Real)(gaussian(generator));
    out += toadd * toadd;
  }
  return out;
}

const Eigen::MatrixXd umatrix(const size_t nrows,
                              const size_t ncols,
                              std::default_random_engine& generator) {
  Eigen::MatrixXd U(nrows, ncols);
  for(size_t i = 0; i < nrows; i++) {
    for(size_t j = 0; j < ncols; j++) {
      U(i, j) = runif(generator);
    }
  }
  return U;
}

template <typename Real>
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
pickCoordinates(
    const size_t Dim,
    const size_t N,
    const size_t fe,
    const std::vector<
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>&
        VT,
    const Eigen::MatrixXd& U) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> VTend(
      Dim, N);
  for(size_t i = 0; i < N; i++) {
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        VTi = VT[i];
    for(size_t j = 0; j < Dim; j++) {
      if(U.coeff(j, i) < 0.5) {
        if(j < fe) {
          VTend(j, i) = VTi.row(j).minCoeff();
        } else {
          Real x = VTi.row(j).minCoeff();
          if(x < 0) {
            x = 0;
          }
          VTend(j, i) = (double)x;
        }
      } else {
        if(j < fe) {
          VTend(j, i) = (double)VTi.row(j).maxCoeff();
        } else {
          Real x = VTi.row(j).maxCoeff();
          if(x < 0) {
            x = 0;
          }
          VTend(j, i) = (double)x;
        }
      }
    }
  }
  return VTend;
}

typedef struct GFI {
  Eigen::MatrixXd vertices;
  Eigen::VectorXd weights;
  Rcpp::NumericVector ess;
} GFI;

template <typename Real>
GFI gfilmm_(
    const Eigen::Matrix<Real, Eigen::Dynamic, 1>& L,
    const Eigen::Matrix<Real, Eigen::Dynamic, 1>& U,
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
        FE,
    const Eigen::SparseMatrix<Real, Eigen::RowMajor>& RE,
    const Eigen::MatrixXi& RE2,
    const Rcpp::IntegerVector E,
    const size_t N,
    const double thresh,
    const unsigned seed,
    const unsigned nthreads) {
  std::default_random_engine generator(seed);
  std::vector<std::default_random_engine> generators(nthreads);
  for(size_t t = 0; t < nthreads; t++) {
    std::default_random_engine gen(seed + (t + 1) * 2000000);
    generators[t] = gen;
  }
  Eigen::Matrix<Real, Eigen::Dynamic, 1> WTnorm(N);  // output:weights
  const size_t n = L.size();
  const size_t fe = FE.cols();  // si FE=NULL, passer une matrice n x 0
  const size_t re = RE2.cols();
  const size_t Dim = fe + re;
  const size_t Dimm1 = Dim - 1;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> VERTEX(
      Dim, N);  // output
  const Rcpp::IntegerVector Esum = Rcpp::cumsum(E);

  //-------- SET-UP ALGORITHM OBJECTS ------------------------------------------
  std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> Z(re);
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> weight =
      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>::Ones(n, N);
  Rcpp::NumericVector ESS(n, (double)N);
  std::vector<int> C;  // initial constraints
  std::vector<int> K;  // complement of C
  std::vector<
      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
      VT(N);  // vertices

  //-------- SAMPLE ALL Z's / SET-UP WEIGHTS -----------------------------------
  std::vector<
      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
      A(N);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
  for(size_t k = 0; k < N; k++) {
    A[k].resize(n, re);
  }

  for(size_t j = 0; j < re; j++) {
    Z[j] = gmatrix<Real>(E(j), N, generator);
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> M =
        RE.block(0, Esum(j) - E(j), n, E(j)) * Z[j];
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
    for(size_t k = 0; k < N; k++) {
      for(size_t i = 0; i < n; i++) {
        A[k](i, j) = M.coeff(i, k);
      }
    }
  }

  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> AA(n,
                                                                          Dim);
  AA << FE, A[0];
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> AT(0, Dim);
  int r = 0;
  for(size_t i = 0; i < n; i++) {
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Atemp(AT.rows() + 1,
                                                              Dim);
    Atemp << AT, AA.row(i);
    if(rankMatrix<Real>(Atemp) > r) {
      AT = Atemp;
      r += 1;
      C.push_back((int)i);
    } else {
      K.push_back((int)i);
    }
  }

  const Eigen::VectorXi K_start =
      Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(C.data(), Dim);
  Eigen::VectorXi JJ =
      Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(C.data(), Dim);

  for(size_t i = 0; i < n - Dim; i++) {
    const int K_i = K[i];
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
    for(size_t j = 0; j < N; j++) {
      Z[re - 1](K_i, j) = 0.0;
    }
  }

  //-------- FIND INITIAL VERTICES ---------------------------------------------
  const std::vector<std::vector<int>> USE = combinations(C, n);
  const size_t twoPowerDim = spow(2, Dim);
  const Eigen::MatrixXi tUSE = vv2matrix(USE, Dim, twoPowerDim);
  std::vector<Eigen::MatrixXi> CC(N, tUSE);  // constraints
  Eigen::Matrix<Real, Eigen::Dynamic, 1> b(2 * n);
  b << U, -L;
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> FEFE(
      2 * n, fe);
  FEFE << FE, -FE;

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
  for(size_t k = 0; k < N; k++) {
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> V(
        Dim, twoPowerDim);
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> AkAk(
        2 * n, re);
    AkAk << A[k], -A[k];
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> AA(
        2 * n, Dim);
    AA << FEFE, AkAk;
    for(size_t i = 0; i < twoPowerDim; i++) {
      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> AAuse(Dim, Dim);
      Eigen::Matrix<Real, Eigen::Dynamic, 1> buse(Dim);
      for(size_t j = 0; j < Dim; j++) {
        buse(j) = b.coeff(USE[i][j]);
        for(size_t l = 0; l < Dim; l++) {
          AAuse(j, l) = AA.coeff(USE[i][j], l);
        }
      }
      V.col(i) = solve<Real>(AAuse, buse);
    }
    VT[k] = V;
  }

  std::vector<size_t> VC(N, twoPowerDim);  // number of vertices

  //-------- MAIN ALGORITHM ----------------------------------------------------
  // double break_point = 10;
  const double lengthK = (double)(n - Dim);
  const size_t K_n = (size_t)(ceil(lengthK / 10.0));
  std::vector<Eigen::VectorXi> K_temp(K_n);
  const Eigen::VectorXi KV =
      Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(K.data(), n - Dim);
  for(size_t i = 0; i < K_n - 1; i++) {
    K_temp[i] = KV.segment(i * 10, 10);
  }
  K_temp[K_n - 1] = KV.tail(n - Dim - (K_n - 1) * 10);
  Eigen::VectorXi K1(0);
  for(size_t k_n = 0; k_n < K_n; k_n++) {
    K1.conservativeResize(K1.size() + K_temp[k_n].size());
    K1.tail(K_temp[k_n].size()) = K_temp[k_n];
    for(int ki = 0; ki < K_temp[k_n].size(); ki++) {
      int k = K_temp[k_n].coeff(ki);
      if(k_n > 0) {
        for(size_t i = 0; i < re; i++) {
          if(E(i) > Z[i].rows()) {
            const int nrowsZi = Z[i].rows();
            Z[i].conservativeResize(E(i), Eigen::NoChange);
            Z[i].bottomRows(E(i) - nrowsZi) =
                gmatrix<Real>(E(i) - nrowsZi, N, generator);
          }
        }
      }

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
      for(size_t i = 0; i < N; i++) {
#ifdef _OPENMP
        const unsigned thread = omp_get_thread_num();
#else
        const unsigned thread = 0;
#endif
        const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic,
                            Eigen::RowMajor>
            VTi = VT[i];
        const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic,
                            Eigen::RowMajor>
            VT1 = VTi.topRows(Dimm1);
        const Eigen::Matrix<Real, Eigen::Dynamic, 1> VT2 = VTi.row(Dimm1);
        Eigen::Matrix<Real, Eigen::Dynamic, 1> Z1t(re - 1);
        for(size_t j = 0; j < re - 1; j++) {
          Z1t(j) = Z[j](RE2.coeff(k, j), i);
        }
        Eigen::Matrix<Real, Eigen::Dynamic, 1> Z1(Dimm1);
        Z1 << FE.row(k).transpose(), Z1t;
        Eigen::Matrix<Real, Eigen::Dynamic, 1> VTsum = VT1.transpose() * Z1;
        std::vector<Real> sample = fidSample<Real>(
            VT2, VTsum, L.coeff(k), U.coeff(k), generators[thread]);
        const Real ZZ = sample[0];
        const Real wt = sample[1];
        Z[re - 1](k, i) = ZZ;
        weight(k, i) = wt;
        VTsum += ZZ * VT2;

        // fidVertex
        const Eigen::MatrixXi CCi = CC[i];
        const Real Lk = L.coeff(k);
        const Real Uk = U.coeff(k);
        const size_t p = VTsum.size();
        std::vector<size_t> whichl, checkl, whichu, checku, both;
        for(size_t h = 0; h < p; h++) {
          bool b = false;
          if(VTsum.coeff(h) >= Lk) {
            whichl.push_back(h);
            b = true;
          } else {
            checkl.push_back(h);
          }
          if(VTsum.coeff(h) <= Uk) {
            whichu.push_back(h);
            if(b) {
              both.push_back(h);
            }
          } else {
            checku.push_back(h);
          }
        }
        Eigen::MatrixXi CCtemp(Dim, 0);
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
            VTtemp(Dim, 0);
        int vert = 0;
        const size_t lcheckl = checkl.size();
        const size_t lwhichl = whichl.size();
        if(lcheckl != 0) {
          Eigen::MatrixXi CA(Dim, lcheckl);
          Eigen::MatrixXi CB(Dim, lwhichl);
          for(size_t h = 0; h < Dim; h++) {
            for(size_t j = 0; j < lcheckl; j++) {
              CA(h, j) = CCi.coeff(h, checkl[j]);
            }
            for(size_t j = 0; j < lwhichl; j++) {
              CB(h, j) = CCi.coeff(h, whichl[j]);
            }
          }
          Eigen::MatrixXi INT = Eigen::MatrixXi::Zero(2 * n, lcheckl);
          for(size_t ll = 0; ll < lcheckl; ll++) {
            for(size_t h = 0; h < Dim; h++) {
              INT(CA.coeff(h, ll), ll) = 1;
            }
          }
          Eigen::Matrix<Real, Eigen::Dynamic, 1> VTsum_cl(lcheckl);
          Eigen::Matrix<Real, Eigen::Dynamic, 1> VTsum_wl(lwhichl);
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VT1_cl(Dim,
                                                                     lcheckl);
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VT1_wl(Dim,
                                                                     lwhichl);
          for(size_t h = 0; h < lcheckl; h++) {
            VTsum_cl(h) = VTsum.coeff(checkl[h]);
            for(size_t j = 0; j < Dim; j++) {
              VT1_cl(j, h) = VTi.coeff(j, checkl[h]);
            }
          }
          for(size_t h = 0; h < lwhichl; h++) {
            VTsum_wl(h) = VTsum.coeff(whichl[h]);
            for(size_t j = 0; j < Dim; j++) {
              VT1_wl(j, h) = VTi.coeff(j, whichl[h]);
            }
          }
          for(size_t ii = 0; ii < p - lcheckl; ii++) {
            Eigen::MatrixXi INT2(Dim, lcheckl);
            for(size_t h = 0; h < Dim; h++) {
              for(size_t j = 0; j < lcheckl; j++) {
                INT2(h, j) = INT(CB.coeff(h, ii), j);
              }
            }
            for(size_t j = 0; j < lcheckl; j++) {
              int colSum = 0;
              for(size_t h = 0; h < Dim; h++) {
                colSum += INT2.coeff(h, j);
              }
              if(colSum == (int)Dimm1) {
                vert += 1;
                std::vector<int> inter(Dim);
                size_t m = 0;
                for(size_t h = 0; h < Dim; h++) {
                  if(INT2(h, j) == 1) {
                    inter[m] = CB.coeff(h, ii);
                    m += 1;
                  }
                }
                inter[Dimm1] = k + (int)n;
                CCtemp.conservativeResize(Eigen::NoChange, vert);
                for(size_t h = 0; h < Dim; h++) {
                  CCtemp(h, vert - 1) = inter[h];
                }
                Real lambda =
                    (Lk - VTsum_wl(ii)) / (VTsum_cl(j) - VTsum_wl(ii));
                VTtemp.conservativeResize(Eigen::NoChange, vert);
                for(size_t h = 0; h < Dim; h++) {
                  VTtemp(h, vert - 1) = lambda * VT1_cl.coeff(h, j) +
                                        (1 - lambda) * VT1_wl.coeff(h, ii);
                }
              }
            }
          }
        }
        size_t lchecku = checku.size();
        size_t lwhichu = whichu.size();
        if(lchecku != 0) {
          Eigen::MatrixXi CA(Dim, lchecku);
          Eigen::MatrixXi CB(Dim, lwhichu);
          for(size_t h = 0; h < Dim; h++) {
            for(size_t j = 0; j < lchecku; j++) {
              CA(h, j) = CCi.coeff(h, checku[j]);
            }
            for(size_t j = 0; j < lwhichu; j++) {
              CB(h, j) = CCi.coeff(h, whichu[j]);
            }
          }
          Eigen::MatrixXi INT = Eigen::MatrixXi::Zero(2 * n, lchecku);
          for(size_t ll = 0; ll < lchecku; ll++) {
            for(size_t h = 0; h < Dim; h++) {
              INT(CA.coeff(h, ll), ll) = 1;
            }
          }
          Eigen::Matrix<Real, Eigen::Dynamic, 1> VTsum_cu(lchecku);
          Eigen::Matrix<Real, Eigen::Dynamic, 1> VTsum_wu(lwhichu);
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VT1_cu(Dim,
                                                                     lchecku);
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VT1_wu(Dim,
                                                                     lwhichu);
          for(size_t h = 0; h < lchecku; h++) {
            VTsum_cu(h) = VTsum.coeff(checku[h]);
            for(size_t j = 0; j < Dim; j++) {
              VT1_cu(j, h) = VTi.coeff(j, checku[h]);
            }
          }
          for(size_t h = 0; h < lwhichu; h++) {
            VTsum_wu(h) = VTsum.coeff(whichu[h]);
            for(size_t j = 0; j < Dim; j++) {
              VT1_wu(j, h) = VTi.coeff(j, whichu[h]);
            }
          }
          for(size_t ii = 0; ii < p - lchecku; ii++) {
            Eigen::MatrixXi INT2(Dim, lchecku);
            for(size_t h = 0; h < Dim; h++) {
              for(size_t j = 0; j < lchecku; j++) {
                INT2(h, j) = INT(CB.coeff(h, ii), j);
              }
            }
            for(size_t j = 0; j < lchecku; j++) {
              int colSum = 0;
              for(size_t h = 0; h < Dim; h++) {
                colSum += INT2.coeff(h, j);
              }
              if(colSum == (int)Dimm1) {
                vert += 1;
                std::vector<int> inter(Dim);
                size_t m = 0;
                for(size_t h = 0; h < Dim; h++) {
                  if(INT2.coeff(h, j) == 1) {
                    inter[m] = CB.coeff(h, ii);
                    m += 1;
                  }
                }
                inter[Dimm1] = k;
                CCtemp.conservativeResize(Eigen::NoChange, vert);
                for(size_t h = 0; h < Dim; h++) {
                  CCtemp(h, vert - 1) = inter[h];
                }
                Real lambda = (Uk - VTsum_wu.coeff(ii)) /
                              (VTsum_cu.coeff(j) - VTsum_wu.coeff(ii));
                VTtemp.conservativeResize(Eigen::NoChange, vert);
                for(size_t h = 0; h < Dim; h++) {
                  VTtemp(h, vert - 1) = lambda * VT1_cu.coeff(h, j) +
                                        (1 - lambda) * VT1_wu.coeff(h, ii);
                }
              }
            }
          }
        }
        size_t lboth = both.size();
        if(lboth > 0) {
          for(size_t j = 0; j < lboth; j++) {
            vert += 1;
            CCtemp.conservativeResize(Eigen::NoChange, vert);
            VTtemp.conservativeResize(Eigen::NoChange, vert);
            for(size_t h = 0; h < Dim; h++) {
              CCtemp(h, vert - 1) = CCi.coeff(h, both[j]);
              VTtemp(h, vert - 1) = VTi.coeff(h, both[j]);
            }
          }
        }
        VC[i] = vert;
        CC[i] = CCtemp;
        VT[i] = VTtemp;
        if(vert == 0) {
          weight(k, i) = 0;
        }
      }

      Eigen::Matrix<Real, Eigen::Dynamic, 1> WT(N);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
      for(size_t j = 0; j < N; j++) {
        WT(j) = weight.col(j).prod();
      }
      const Real WTsum = WT.sum();
      if(WTsum == 0) {
        Rcpp::Rcout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        Rcpp::Rcout << "!!!   Error: possible underflow. !!!\n";
        Rcpp::Rcout << "!!!   Potential resolutions:     !!!\n";
        Rcpp::Rcout << "!!!   1.  Run more simulations   !!!\n";
        Rcpp::Rcout << "!!!   2.  Widen data intervals   !!!\n";
        Rcpp::Rcout << "!!!   3.  Center data            !!!\n";
        Rcpp::Rcout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        throw Rcpp::exception("Something bad occured.");
      }
      WTnorm = WT / WTsum;
      ESS(k) = (double)(1 / (WTnorm.dot(WTnorm)));

      JJ.conservativeResize(JJ.size() + 1);
      const int lenJJ = JJ.size();
      JJ(lenJJ - 1) = k;

      if(ESS(k) < thresh && k < K[n - Dim - 1]) {
        std::vector<size_t> N_sons(N, 0);
        const std::vector<Real> dist = Vcumsum<Real>(N * WTnorm);
        const double aux = runif(generator);
        size_t j = 0;
        for(size_t i = 0; i < N; i++) {
          const double u = aux + i;
          while(u > dist[j]) {
            j += 1;
          }
          N_sons[j] += 1;
        }
        Eigen::SparseMatrix<Real, Eigen::RowMajor> REJJ(lenJJ, Esum(re - 1));
        for(int jj = 0; jj < lenJJ; jj++) {
          REJJ.row(jj) = RE.row(JJ.coeff(jj));
        }
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
            FEJJ(lenJJ, fe);
        if(fe > 0) {
          for(int jj = 0; jj < lenJJ; jj++) {
            FEJJ.row(jj) = FE.row(JJ.coeff(jj));
          }
        }
        std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> ZZ(re);
        std::vector<size_t> Nsons_sum(N);
        Nsons_sum[0] = 0;
        for(size_t i = 1; i < N; i++) {
          Nsons_sum[i] = Nsons_sum[i - 1] + N_sons[i - 1];
        }
        for(size_t ii = 0; ii < re; ii++) {
          ZZ[ii].resize(E(ii), Nsons_sum[N - 1] + N_sons[N - 1]);
        }
        std::vector<size_t> VCVC(N, 0);
        std::vector<Eigen::MatrixXi> CCCC(N);
        std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic,
                                  Eigen::RowMajor>>
            VTVT(N);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
        for(size_t i = 0; i < N; i++) {
#ifdef _OPENMP
          const unsigned thread = omp_get_thread_num();
#else
          const unsigned thread = 0;
#endif
          const size_t Nsons_i = N_sons[i];
          if(Nsons_i) {
            std::vector<size_t> VCtemp(Nsons_i, VC[i]);
            std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>>
                Ztemp(re);
            const size_t copy = Nsons_i - 1;
            const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic,
                                Eigen::RowMajor>
                VTi = VT[i];
            std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic,
                                      Eigen::RowMajor>>
                VTtemp(Nsons_i, VTi);
            for(size_t ii = 0; ii < re; ii++) {
              Ztemp[ii].resize(E(ii), Nsons_i);
              for(size_t j = 0; j < Nsons_i; j++) {
                Ztemp[ii].col(j) = Z[ii].col(i);
              }
            }
            if(copy) {
              const std::vector<size_t> ord =
                  sample_int(re, generators[thread]);
              for(size_t j = 0; j < re; j++) {
                const size_t kk = ord[j];
                for(size_t ii = 0; ii < copy; ii++) {
                  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> XX(lenJJ,
                                                                         0);
                  for(size_t jj = 0; jj < re; jj++) {
                    if(jj != kk) {
                      XX.conservativeResize(Eigen::NoChange, XX.cols() + 1);
                      const Eigen::Matrix<Real, Eigen::Dynamic, 1> newcol =
                          REJJ.block(0, Esum(jj) - E(jj), lenJJ, E(jj)) *
                          Ztemp[jj].col(ii);
                      XX.rightCols(1) = newcol;
                    }
                  }
                  Eigen::VectorXi v(lenJJ);
                  for(int jj = 0; jj < lenJJ; jj++) {
                    v(jj) = RE2(JJ.coeff(jj), kk);
                  }
                  const Eigen::VectorXi vv = cppunique(v);
                  Eigen::Matrix<Real, Eigen::Dynamic, 1> Z1_(vv.size());
                  for(int jj = 0; jj < vv.size(); jj++) {
                    Z1_(jj) = Ztemp[kk].coeff(vv(jj), ii);
                  }
                  const Eigen::SparseMatrix<Real> CO2_ =
                      REJJ.block(0, Esum(kk) - E(kk), lenJJ, E(kk));
                  std::vector<int> pcolsums;
                  for(int jj = 0; jj < E(kk); jj++) {
                    const Real colsum = CO2_.col(jj).sum();
                    if(colsum > 0) {
                      pcolsums.push_back(jj);
                    }
                  }
                  const size_t ncolsCO2 = pcolsums.size();
                  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> CO2(
                      lenJJ, ncolsCO2);
                  for(size_t jj = 0; jj < ncolsCO2; jj++) {
                    CO2.col(jj) = CO2_.col(pcolsums[jj]);
                  }
                  std::vector<int> Z00;
                  for(int jj = 0; jj < Z1_.size(); jj++) {
                    if(Z1_.coeff(jj) != 0) {
                      Z00.push_back(jj);
                    }
                  }
                  const size_t lenZ1 = Z00.size();
                  Eigen::Matrix<Real, Eigen::Dynamic, 1> Z1(lenZ1);
                  for(size_t jj = 0; jj < lenZ1; jj++) {
                    Z1(jj) = Z1_.coeff(Z00[jj]);
                  }
                  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> XXX(
                      lenJJ, Dimm1);
                  if(fe > 0) {
                    XXX << FEJJ, XX;
                  } else {
                    XXX = XX;
                  }
                  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> MAT(
                      lenJJ, ((int)Dim) - 1 + ncolsCO2);
                  MAT << -XXX, CO2;
                  const Eigen::FullPivLU<
                      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>>
                      lu(MAT);
                  const int rk = lu.rank();
                  if(rk < MAT.cols()) {
                    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>
                        NUL = lu.kernel();
                    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>
                        n1 = NUL.topRows(NUL.rows() - ncolsCO2);
                    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>
                        n2 = NUL.bottomRows(ncolsCO2);
                    const std::vector<
                        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>>
                        QR = QRdecomp<Real>(n2);
                    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>
                        O2 = QR[0];
                    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>
                        R = QR[1];
                    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>
                        O1 = tsolveAndMultiply<Real>(R, n1);
                    const int rankO2 = O2.cols();
                    const Eigen::Matrix<Real, Eigen::Dynamic, 1> a =
                        O2.transpose() * Z1;
                    const Eigen::Matrix<Real, Eigen::Dynamic, 1> O2a = O2 * a;
                    const Eigen::Matrix<Real, Eigen::Dynamic, 1> tau_ =
                        Z1 - O2a;
                    const Real b = sqrt(tau_.dot(tau_));
                    const Eigen::Matrix<Real, Eigen::Dynamic, 1> tau = tau_ / b;
                    const Real bb =
                        sqrt(rchisq<Real>(lenZ1 - rankO2, generators[thread]));
                    const Real bbb = b / bb;
                    Eigen::Matrix<Real, Eigen::Dynamic, 1> aa(rankO2);
                    for(int jj = 0; jj < rankO2; jj++) {
                      aa(jj) = (Real)(gaussian(generators[thread]));
                    }
                    const Eigen::Matrix<Real, Eigen::Dynamic, 1> O2aa = O2 * aa;
                    const Eigen::Matrix<Real, Eigen::Dynamic, 1> bbtau =
                        bb * tau;
                    const Eigen::Matrix<Real, Eigen::Dynamic, 1> MM3 =
                        O2aa + bbtau;
                    for(size_t jj = 0; jj < lenZ1; jj++) {
                      Ztemp[kk](Z00[jj], ii) = MM3(jj);
                    }
                    const Eigen::Matrix<Real, Eigen::Dynamic, 1> M2 =
                        O1 * (bbb * aa - a);
                    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>
                        VTtemp_ii = VTtemp[ii];
                    for(size_t jj = 0; jj < VC[i]; jj++) {
                      for(size_t vert = 0; vert < Dim; vert++) {
                        if(vert < fe + kk) {
                          VTtemp[ii](vert, jj) =
                              VTtemp_ii.coeff(vert, jj) -
                              VTtemp_ii.coeff(fe + kk, jj) * M2.coeff(vert);
                        } else if(vert > fe + kk) {
                          VTtemp[ii](vert, jj) =
                              VTtemp_ii.coeff(vert, jj) -
                              VTtemp_ii.coeff(fe + kk, jj) * M2.coeff(vert - 1);
                        }
                      }
                      VTtemp[ii](fe + kk, jj) =
                          bbb * VTtemp_ii.coeff(fe + kk, jj);
                    }
                  } else {
                    const Real b = sqrt(Z1.dot(Z1));
                    const Eigen::Matrix<Real, Eigen::Dynamic, 1> tau = Z1 / b;
                    const Real bb =
                        sqrt(rchisq<Real>(lenZ1, generators[thread]));
                    for(size_t jj = 0; jj < lenZ1; jj++) {
                      Ztemp[kk](Z00[jj], ii) = bb * tau.coeff(jj);
                    }
                    VTtemp[ii].row(fe + kk) *= b / bb;
                  }
                }
              }
            }
            for(size_t ii = 0; ii < re; ii++) {
              ZZ[ii].block(0, Nsons_sum[i], E(ii), Nsons_i) = Ztemp[ii];
            }
            size_t d = 0;
            for(size_t ii = 0; ii < i; ii++) {
              d += N_sons[ii];
            }
            const size_t dd = d + Nsons_i;
            for(size_t ii = d; ii < dd; ii++) {
              VCVC[ii] = VCtemp[ii - d];
            }
            for(size_t kk = 0; kk < Nsons_i; kk++) {
              VTVT[kk + d] = VTtemp[kk];
              CCCC[kk + d] = CC[i];
            }
          }
        }
        Z = ZZ;
        VT = VTVT;
        VC = VCVC;
        CC = CCCC;
        weight =
            Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>::Ones(n, N);
      }
    }

    //------------determine signs ----------------------------------------------
    Eigen::MatrixXi signs = Eigen::MatrixXi::Zero(re, N);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
    for(size_t i = 0; i < N; i++) {
      const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
          VTi = VT[i];
      for(size_t j = 0; j < re; j++) {
        const Eigen::Matrix<Real, Eigen::Dynamic, 1> row = VTi.row(fe + j);
        bool allneg = true;
        size_t jj = 0;
        while(jj < VC[i] && allneg) {
          allneg = row.coeff(jj) < 0;
          jj += 1;
        }
        if(allneg) {
          signs(j, i) = -1;
        }
      }
    }

    //--------FINAL RESAMPLE ---------------------------------------------------
    std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> ZZ(
        re);  // resized below
    std::vector<Eigen::VectorXi> nn(re);
    std::vector<int> lengths_nn(re);
    std::vector<
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        VTVT(N);
    Eigen::VectorXi n1_(K1.size() + (int)Dim);
    n1_ << K1, K_start;
    const Eigen::VectorXi n1 = Vsort(n1_);
    for(size_t ii = 0; ii < re; ii++) {
      Eigen::VectorXi vec(n1.size());
      for(int jj = 0; jj < n1.size(); jj++) {
        vec(jj) = RE2.coeff(n1(jj), ii);
      }
      nn[ii] = cppunique(vec);
      lengths_nn[ii] = nn[ii].size();
      ZZ[ii].resize(lengths_nn[ii], N);
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
    for(size_t i = 0; i < N; i++) {
#ifdef _OPENMP
      const unsigned thread = omp_get_thread_num();
#else
      const unsigned thread = 0;
#endif
      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
          VTtemp = VT[i];
      std::vector<Eigen::Matrix<Real, Eigen::Dynamic, 1>> Ztemp(re);
      for(size_t ii = 0; ii < re; ii++) {
        Ztemp[ii].resize(lengths_nn[ii]);
        for(int iii = 0; iii < lengths_nn[ii]; iii++) {
          Ztemp[ii](iii) = Z[ii].coeff(nn[ii].coeff(iii), i);
        }
      }
      const std::vector<size_t> ord = sample_int(re, generators[thread]);
      for(size_t j = 0; j < re; j++) {
        const size_t kk = ord[j];
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> XX(n, 0);
        for(size_t jj = 0; jj < re; jj++) {
          if(jj != kk) {
            XX.conservativeResize(Eigen::NoChange, XX.cols() + 1);
            const Eigen::Matrix<Real, Eigen::Dynamic, 1> newcol =
                RE.block(0, Esum(jj) - E(jj), n, lengths_nn[jj]) * Ztemp[jj];
            XX.rightCols(1) = newcol;
          }
        }
        const Eigen::Matrix<Real, Eigen::Dynamic, 1> Z1 = Ztemp[kk];
        const int lenZ1 = Z1.size();
        const int ncolsCO2 = lengths_nn[kk];
        const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> CO2 =
            RE.block(0, Esum(kk) - E(kk), n, ncolsCO2);
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> XXX(n, Dimm1);
        if(fe > 0) {
          XXX << FE, XX;
        } else {
          XXX = XX;
        }
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> MAT(
            n, Dimm1 + ncolsCO2);
        MAT << -XXX, CO2;
        const Eigen::FullPivLU<
            Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>>
            lu(MAT);
        const int rk = lu.rank();
        if(rk < MAT.cols()) {
          const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> NUL =
              lu.kernel();
          const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> n1 =
              NUL.topRows(NUL.rows() - ncolsCO2);
          const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> n2 =
              NUL.bottomRows(ncolsCO2);
          const std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>>
              QR = QRdecomp<Real>(n2);
          const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> O2 = QR[0];
          const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> R = QR[1];
          const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> O1 =
              tsolveAndMultiply<Real>(R, n1);
          const int rankO2 = O2.cols();
          const Eigen::Matrix<Real, Eigen::Dynamic, 1> a = O2.transpose() * Z1;
          const Eigen::Matrix<Real, Eigen::Dynamic, 1> O2a = O2 * a;
          const Eigen::Matrix<Real, Eigen::Dynamic, 1> tau_ = Z1 - O2a;
          const Real b = sqrt(tau_.dot(tau_));
          const Eigen::Matrix<Real, Eigen::Dynamic, 1> tau = tau_ / b;
          const Real bb =
              sqrt(rchisq<Real>(lenZ1 - rankO2, generators[thread]));
          const Real bbb = b / bb;
          Eigen::Matrix<Real, Eigen::Dynamic, 1> aa(rankO2);
          for(int jj = 0; jj < rankO2; jj++) {
            aa(jj) = (Real)(gaussian(generators[thread]));
          }
          const Eigen::Matrix<Real, Eigen::Dynamic, 1> O2aa = O2 * aa;
          const Eigen::Matrix<Real, Eigen::Dynamic, 1> bbtau = bb * tau;
          Ztemp[kk] = O2aa + bbtau;
          const Eigen::Matrix<Real, Eigen::Dynamic, 1> M2 = O1 * (bbb * aa - a);
          for(size_t jj = 0; jj < VC[i]; jj++) {
            for(size_t vert = 0; vert < Dim; vert++) {
              if(vert < fe + kk) {
                VTtemp(vert, jj) -= VTtemp.coeff(fe + kk, jj) * M2.coeff(vert);
              } else if(vert > fe + kk) {
                VTtemp(vert, jj) -=
                    VTtemp.coeff(fe + kk, jj) * M2.coeff(vert - 1);
              }
            }
            VTtemp(fe + kk, jj) *= bbb;
          }
        } else {
          const Real b = sqrt(Z1.dot(Z1));
          const Eigen::Matrix<Real, Eigen::Dynamic, 1> tau = Z1 / b;
          const Real bb = sqrt(rchisq<Real>(lenZ1, generators[thread]));
          Ztemp[kk] = bb * tau;
          VTtemp.row(fe + kk) *= b / bb;
        }
      }
      for(size_t ii = 0; ii < re; ii++) {
        ZZ[ii].col(i) = Ztemp[ii];
      }
      VTVT[i] = VTtemp;
    }

    Z = ZZ;
    VT = VTVT;

    //---------flip negatives --------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
    for(size_t i = 0; i < N; i++) {
      for(size_t j = 0; j < re; j++) {
        if(signs(j, i) == -1) {
          VT[i].row(fe + j) *= -1;
          Z[j].col(i) *= -1;
        }
      }
    }

    if(k_n == K_n - 1) {  // if finished pick coordinates
      const Eigen::MatrixXd unif = umatrix(Dim, N, generator);
      VERTEX = pickCoordinates<Real>(Dim, N, fe, VT, unif);
    }
  }

  GFI out;
  out.weights = WTnorm.template cast<double>();
  out.vertices = VERTEX;
  out.ess = ESS;

  return out;
}

// [[Rcpp::export]]
Rcpp::List gfilmm_double(
    const Eigen::VectorXd& L,
    const Eigen::VectorXd& U,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& FE,
    const Eigen::SparseMatrix<double>& RE,
    const Eigen::MatrixXi& RE2,
    const Rcpp::IntegerVector E,
    const size_t N,
    const double thresh,
    const unsigned seed,
    const unsigned nthreads) {
  const GFI gfi =
      gfilmm_<double>(L, U, FE, RE, RE2, E, N, thresh, seed, nthreads);

  Rcpp::List out = Rcpp::List::create(Rcpp::Named("VERTEX") = gfi.vertices,
                                      Rcpp::Named("WEIGHT") = gfi.weights);
  out.attr("ESS") = gfi.ess;
  return out;
}

// [[Rcpp::export]]
Rcpp::List gfilmm_long(
    const Eigen::VectorXd& L,
    const Eigen::VectorXd& U,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& FE,
    const Eigen::SparseMatrix<double>& RE,
    const Eigen::MatrixXi& RE2,
    const Rcpp::IntegerVector E,
    const size_t N,
    const double thresh,
    const unsigned seed,
    const unsigned nthreads) {
  const GFI gfi = gfilmm_<long double>(
      L.cast<long double>(), U.cast<long double>(), FE.cast<long double>(),
      RE.cast<long double>(), RE2, E, N, thresh, seed, nthreads);

  Rcpp::List out = Rcpp::List::create(Rcpp::Named("VERTEX") = gfi.vertices,
                                      Rcpp::Named("WEIGHT") = gfi.weights);
  out.attr("ESS") = gfi.ess;
  return out;
}
