#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// #include <conio.h>

//   -*-   -*-   -*-

unsigned seed = 0;

void SetRNGSeed(unsigned newSeed) {
    seed = newSeed;
}

int RandomInteger() {
    unsigned prevSeed = seed;
    seed = (1664525 * seed) + 1013904223;
    return prevSeed;
}

float RandomFloat() {
    float tmp = 0.0f;
    *((unsigned *)&tmp) = (seed & 8388607) | 1065353216;
    seed = (1664525 * seed) + 1013904223;
    return tmp - 1.0f;
}

//   -*-   -*-   -*-

bool CholeskyRankOneUpdate(mat &L, vec v) {
    double *vPtr = v.memptr();
    for (int i = 0; i < L.n_rows; ++i) {
        double *LPtr = L.colptr(i);
        double r = sqrt((LPtr[i] * LPtr[i]) + (vPtr[i] * vPtr[i]));
        double c = LPtr[i] / r;
        double s = -vPtr[i] / r;
        LPtr[i] = r;
        for (int j = i + 1; j < L.n_cols; ++j) {
            double tmp1 = (LPtr[j] * c) - (vPtr[j] * s);
            double tmp2 = (LPtr[j] * s) + (vPtr[j] * c);
            LPtr[j] = tmp1;
            vPtr[j] = tmp2;
        }
    }
    return true;
}

//   -*-   -*-   -*-

bool CholeskyRankOneDowndate(mat &L, vec v) {
    double *vPtr = v.memptr();
    for (int i = 0; i < L.n_rows; ++i) {
        double *LPtr = L.colptr(i);
        double r = sqrt((LPtr[i] * LPtr[i]) - (vPtr[i] * vPtr[i]));
        double c = LPtr[i] / r;
        double s = -vPtr[i] / r;
        LPtr[i] = r;
        for (int j = i + 1; j < L.n_cols; ++j) {
            double tmp1 = (LPtr[j] * c) + (vPtr[j] * s);
            double tmp2 = (LPtr[j] * s) + (vPtr[j] * c);
            LPtr[j] = tmp1;
            vPtr[j] = tmp2;
        }
    }
    return true;
}

//   -*-   -*-   -*-

// [[Rcpp::export]]
void UpdateMeansForQuadraticFunction(Rcpp::List &res) {
    Rcpp::List mL = res["means"];
    Rcpp::List dirsL = res["directions"];
    Rcpp::List coeffsL = res["coefficients"];
    int maxClusters = mL.length();
    for (int i = 0; i < maxClusters; ++i) {
        if (mL[i] != R_NilValue) {
            vec m = mL[i];
            int dir = dirsL[i];
            vec coeffs = coeffsL[i];
            vec tmp(m.n_elem - 1);
            if (dir > 0) tmp.subvec(0, dir - 1) = m.subvec(0, dir - 1);
            if (dir < m.n_elem - 1) tmp.subvec(dir, m.n_elem - 2) = m.subvec(dir + 1, m.n_elem - 1);
            double value = 0.0;
            for (int k = 0; k < (coeffs.n_elem - 1) / 2; ++k) {
                value += coeffs(k) * tmp(k) * tmp(k);
                value += coeffs(((coeffs.n_elem - 1) / 2) + k) * tmp(k);
            }
            value += coeffs(coeffs.n_elem - 1);
            m[dir] = value;
            mL[i] = m;
        }
    }
    res["means"] = mL;
}

//   -*-   -*-   -*-

// [[Rcpp::export]]
Rcpp::List CalculateEllipsesOfConfidenceForQuadraticFunction(Rcpp::List res, double confidence, int segments) {
    Rcpp::List covsL = res["covariances"];
    Rcpp::List mL = res["means"];
    Rcpp::List dirsL = res["directions"];
    Rcpp::List coeffsL = res["coefficients"];
    int maxClusters = covsL.length();
    Rcpp::List ellipsesL(maxClusters);
    Rcpp::List axesL(maxClusters);
    for (int i = 0; i < maxClusters; ++i) {
        if (covsL[i] != R_NilValue) {
            mat cov = covsL[i];
            vec m = mL[i];
            int dir = dirsL[i];
            vec coeffs = coeffsL[i];
            mat points(segments + 1, 2);
            double lambda1 = cov(0, 0);
            double lambda2 = cov(1, 1);
            vec eAxes(2);
            eAxes(0) = sqrt(lambda1 * R::qchisq(confidence, 2, 1, 0));
            eAxes(1) = sqrt(lambda2 * R::qchisq(confidence, 2, 1, 0));
            for (int j = 0; j < segments; ++j) {
                double x = eAxes(0) * cos(((2.0 * M_PI) / segments) * j);
                double y = eAxes(1) * sin(((2.0 * M_PI) / segments) * j);
                points(j, 0) = x;
                points(j, 1) = y;
                points(j, (dir + 1) & 1) += m((dir + 1) & 1);
                double tmp = points(j, (dir + 1) & 1);
                double value = (coeffs(0) * tmp * tmp) + (coeffs(1) * tmp) + coeffs(2);
                points(j, dir) += value;
            }
            points(segments, 0) = points(0, 0);
            points(segments, 1) = points(0, 1);
            ellipsesL[i] = points;
            for (int j = 0; j <= segments; ++j) {
                points(j, (dir + 1) & 1) = -eAxes((dir + 1) & 1) + (((2.0 * eAxes((dir + 1) & 1)) / segments) * j) + m((dir + 1) & 1);
                double tmp = points(j, (dir + 1) & 1);
                double value = (coeffs(0) * tmp * tmp) + (coeffs(1) * tmp) + coeffs(2);
                points(j, dir) = value;
            }
            axesL[i] = points;
        } else {
            ellipsesL[i] = R_NilValue;
            axesL[i] = R_NilValue;
        }
    };
    return Rcpp::List::create(ellipsesL, axesL);
}

//   -*-   -*-   -*-

// [[Rcpp::export]]
Rcpp::List CalculateEllipsoidsOfConfidenceForQuadraticFunction(Rcpp::List res, double confidence, int gridRes) {
    Rcpp::List covsL = res["covariances"];
    Rcpp::List mL = res["means"];
    Rcpp::List dirsL = res["directions"];
    Rcpp::List coeffsL = res["coefficients"];
    int maxClusters = covsL.length();
    Rcpp::List ellipsoidsL(maxClusters);
    Rcpp::List ellipsoidsNL(maxClusters);
    Rcpp::List ellipsoidsML(maxClusters);
    Rcpp::List ellipsoidsMNL(maxClusters);
    for (int i = 0; i < maxClusters; ++i) {
        if (covsL[i] != R_NilValue) {
            vec eigval;
            mat eigvec;
            mat cov = covsL[i];
            eig_sym(eigval, eigvec, cov);
            double percentile = R::qchisq(confidence, 3, 1, 0);
            vec eAxes(3);
            for (int j = 0; j < 3; ++j) eAxes(j) = sqrt(eigval(j) * percentile);
            vec m = mL[i];
            int dir = dirsL[i];
            vec coeffs = coeffsL[i];

            mat verts(2 + ((gridRes - 1) * gridRes), 3);
            verts(0, 0) = 0.0; verts(0, 1) = 0.0; verts(0, 2) = eAxes(2);
            for (int j = 1; j < gridRes; ++j) {
                for (int k = 0; k < gridRes; ++k) {
                    verts(1 + ((j - 1) * gridRes) + k, 0) = eAxes(0) * sin((M_PI / gridRes) * j) * cos(((2 * M_PI) / gridRes) * k);
                    verts(1 + ((j - 1) * gridRes) + k, 1) = eAxes(1) * sin((M_PI / gridRes) * j) * sin(((2 * M_PI) / gridRes) * k);
                    verts(1 + ((j - 1) * gridRes) + k, 2) = eAxes(2) * cos((M_PI / gridRes) * j);
                }
            }
            verts(1 + ((gridRes - 1) * gridRes), 0) = 0.0;
            verts(1 + ((gridRes - 1) * gridRes), 1) = 0.0;
            verts(1 + ((gridRes - 1) * gridRes), 2) = -eAxes(2);

            verts = (eigvec * verts.t()).t();
            for (int j = 0; j < verts.n_rows; ++j) {
                vec tmp = verts.row(j).t();
                for (int k = 0; k < 3; ++k) {
                    if (k != dir) tmp(k) += m(k);
                }
                verts.row(j) = tmp.t();
                tmp.shed_row(dir);
                double value = (coeffs(0) * tmp(0) * tmp(0)) + (coeffs(1) * tmp(1) * tmp(1)) +
                    (coeffs(2) * tmp(0)) + (coeffs(3) * tmp(1)) +
                    coeffs(4);
                verts(j, dir) += value;
            }

            umat faces((gridRes * 2) + (((gridRes - 2) * gridRes) * 2), 3);
            for (int j = 0; j < gridRes; ++j) {
                faces(j, 0) = 0;
                faces(j, 1) = 1 + ((j + 1) % gridRes);
                faces(j, 2) = 1 + j;
            }
            for (int j = 0; j < gridRes - 2; ++j) {
                for (int k = 0; k < gridRes; ++k) {
                    faces(gridRes + (((j * gridRes) + k) * 2), 0) = 1 + (j * gridRes) + k;
                    faces(gridRes + (((j * gridRes) + k) * 2), 1) = 1 + (j * gridRes) + ((k + 1) % gridRes);
                    faces(gridRes + (((j * gridRes) + k) * 2), 2) = 1 + ((j + 1) * gridRes) + ((k + 1) % gridRes);
                    faces(gridRes + (((j * gridRes) + k) * 2) + 1, 0) = 1 + (j * gridRes) + k;
                    faces(gridRes + (((j * gridRes) + k) * 2) + 1, 1) = 1 + ((j + 1) * gridRes) + ((k + 1) % gridRes);
                    faces(gridRes + (((j * gridRes) + k) * 2) + 1, 2) = 1 + ((j + 1) * gridRes) + k;
                }
            }
            for (int j = 0; j < gridRes; ++j) {
                faces(gridRes + (((gridRes - 2) * gridRes) * 2) + j, 0) = 1 + ((gridRes - 2) * gridRes) + j;
                faces(gridRes + (((gridRes - 2) * gridRes) * 2) + j, 1) = 1 + ((gridRes - 2) * gridRes) + ((j + 1) % gridRes);
                faces(gridRes + (((gridRes - 2) * gridRes) * 2) + j, 2) = 1 + ((gridRes - 1) * gridRes);
            }

            mat fNormals(faces.n_rows, 3);
            for (int j = 0; j < faces.n_rows; ++j) {
                unsigned i1 = faces(j, 0);
                unsigned i2 = faces(j, 1);
                unsigned i3 = faces(j, 2);
                vec v1 = verts.row(i1).t();
                vec v2 = verts.row(i2).t();
                vec v3 = verts.row(i3).t();
                vec N = normalise(cross(v2 - v1, v3 - v1));
                fNormals.row(j) = N.t();
            }

            mat vNormals(verts.n_rows, 3, fill::zeros);
            for (int j = 0; j < faces.n_rows; ++j) {
                unsigned i1 = faces(j, 0);
                unsigned i2 = faces(j, 1);
                unsigned i3 = faces(j, 2);
                vNormals.row(i1) += fNormals.row(j);
                vNormals.row(i2) += fNormals.row(j);
                vNormals.row(i3) += fNormals.row(j);
            }

            for (int j = 0; j < vNormals.n_rows; ++j) vNormals.row(j) = normalise(vNormals.row(j));

            mat eFaces(faces.n_rows * 3, 3);
            mat eNormals(faces.n_rows * 3, 3);
            for (int j = 0; j < faces.n_rows; ++j) {
                unsigned i1 = faces(j, 0);
                unsigned i2 = faces(j, 1);
                unsigned i3 = faces(j, 2);
                eFaces.row(j * 3) = verts.row(i1);
                eFaces.row((j * 3) + 1) = verts.row(i2);
                eFaces.row((j * 3) + 2) = verts.row(i3);
                eNormals.row(j * 3) = vNormals.row(i1);
                eNormals.row((j * 3) + 1) = vNormals.row(i2);
                eNormals.row((j * 3) + 2) = vNormals.row(i3);
            }

            ellipsoidsL[i] = eFaces;
            ellipsoidsNL[i] = eNormals;

            cov.shed_row(dir);
            cov.shed_col(dir);
            eig_sym(eigval, eigvec, cov);
            for (int j = 0; j < 2; ++j) eAxes(j) = sqrt(eigval(j) * percentile);

            mat vertsM(1 + (gridRes * gridRes), 2);
            verts(0, 0) = 0.0; verts(0, 1) = 0.0;
            for (int j = 1; j <= gridRes; ++j) {
                for (int k = 0; k < gridRes; ++k) {
                    vertsM(1 + ((j - 1) * gridRes) + k, 0) = ((eAxes(0) / gridRes) * j) * cos(((2 * M_PI) / gridRes) * k);
                    vertsM(1 + ((j - 1) * gridRes) + k, 1) = ((eAxes(1) / gridRes) * j) * sin(((2 * M_PI) / gridRes) * k);
                }
            }

            vertsM = (eigvec * vertsM.t()).t();
            vertsM.insert_cols(dir, 1, true);
            for (int j = 0; j < vertsM.n_rows; ++j) {
                vec tmp = vertsM.row(j).t();
                for (int k = 0; k < 3; ++k) {
                    if (k != dir) tmp(k) += m(k);
                }
                vertsM.row(j) = tmp.t();
                tmp.shed_row(dir);
                double value = (coeffs(0) * tmp(0) * tmp(0)) + (coeffs(1) * tmp(1) * tmp(1)) +
                    (coeffs(2) * tmp(0)) + (coeffs(3) * tmp(1)) +
                    coeffs(4);
                vertsM(j, dir) += value;
            }

            umat facesM(gridRes + (((gridRes - 1) * gridRes) * 2), 3);
            for (int j = 0; j < gridRes; ++j) {
                facesM(j, 0) = 0;
                facesM(j, 1) = 1 + ((j + 1) % gridRes);
                facesM(j, 2) = 1 + j;
            }
            for (int j = 0; j < gridRes - 1; ++j) {
                for (int k = 0; k < gridRes; ++k) {
                    facesM(gridRes + (((j * gridRes) + k) * 2), 0) = 1 + (j * gridRes) + k;
                    facesM(gridRes + (((j * gridRes) + k) * 2), 1) = 1 + (j * gridRes) + ((k + 1) % gridRes);
                    facesM(gridRes + (((j * gridRes) + k) * 2), 2) = 1 + ((j + 1) * gridRes) + ((k + 1) % gridRes);
                    facesM(gridRes + (((j * gridRes) + k) * 2) + 1, 0) = 1 + (j * gridRes) + k;
                    facesM(gridRes + (((j * gridRes) + k) * 2) + 1, 1) = 1 + ((j + 1) * gridRes) + ((k + 1) % gridRes);
                    facesM(gridRes + (((j * gridRes) + k) * 2) + 1, 2) = 1 + ((j + 1) * gridRes) + k;
                }
            }

            mat fNormalsM(facesM.n_rows, 3);
            for (int j = 0; j < facesM.n_rows; ++j) {
                unsigned i1 = facesM(j, 0);
                unsigned i2 = facesM(j, 1);
                unsigned i3 = facesM(j, 2);
                vec v1 = vertsM.row(i1).t();
                vec v2 = vertsM.row(i2).t();
                vec v3 = vertsM.row(i3).t();
                vec N = normalise(cross(v2 - v1, v3 - v1));
                fNormalsM.row(j) = N.t();
            };

            mat vNormalsM(vertsM.n_rows, 3, fill::zeros);
            for (int j = 0; j < facesM.n_rows; ++j) {
                unsigned i1 = facesM(j, 0);
                unsigned i2 = facesM(j, 1);
                unsigned i3 = facesM(j, 2);
                vNormalsM.row(i1) += fNormalsM.row(j);
                vNormalsM.row(i2) += fNormalsM.row(j);
                vNormalsM.row(i3) += fNormalsM.row(j);
            }

            for (int j = 0; j < vNormalsM.n_rows; ++j) vNormalsM.row(j) = normalise(vNormalsM.row(j));

            mat eMFaces(facesM.n_rows * 3, 3);
            mat eMNormals(facesM.n_rows * 3, 3);
            for (int j = 0; j < facesM.n_rows; ++j) {
                unsigned i1 = facesM(j, 0);
                unsigned i2 = facesM(j, 1);
                unsigned i3 = facesM(j, 2);
                eMFaces.row(j * 3) = vertsM.row(i1);
                eMFaces.row((j * 3) + 1) = vertsM.row(i2);
                eMFaces.row((j * 3) + 2) = vertsM.row(i3);
                eMNormals.row(j * 3) = vNormalsM.row(i1);
                eMNormals.row((j * 3) + 1) = vNormalsM.row(i2);
                eMNormals.row((j * 3) + 2) = vNormalsM.row(i3);
            }

            ellipsoidsML[i] = eMFaces;
            ellipsoidsMNL[i] = eMNormals;
        } else {
            ellipsoidsL[i] = R_NilValue;
            ellipsoidsNL[i] = R_NilValue;
            ellipsoidsML[i] = R_NilValue;
            ellipsoidsMNL[i] = R_NilValue;
        }
    }
    return Rcpp::List::create(ellipsoidsL, ellipsoidsNL, ellipsoidsML, ellipsoidsMNL);
}

//   -*-   -*-   -*-

SEXP afCECLloyd(
    const mat &points,
    int maxClusters,
    ivec labels,
    double cardMin,
    double costThreshold,
    int minIterations,
    int maxIterations,
    const mat &values,
    bool interactive
) {
    int *card = NULL;
    vec **m = NULL;
    mat **sigma = NULL;
    mat ***A = NULL;
    vec ***b = NULL;
    double **sumOfSquares = NULL;
    double *E = NULL;
    int *bestAxes = NULL;
    vec **coeffs = NULL;
    double *probabilities = NULL;

    int pointsNum = points.n_cols;
    int dim = points.n_rows;
    int dimA = values.n_rows / dim;

    Rcpp::List res;

    card = new int[maxClusters];
    std::set<int> activeClusters;
    std::set<int> inactiveClusters;
    for (int i = 0; i < maxClusters; ++i) card[i] = 0;
    for (int i = 0; i < pointsNum; ++i) ++card[labels[i]];
    for (int i = 0; i < maxClusters; ++i) {
        if ((card[i] >= dim + 1) && (((double)card[i]) / pointsNum >= cardMin)) activeClusters.insert(activeClusters.end(), i);
        else {
            if (card[i] > 0) inactiveClusters.insert(inactiveClusters.end(), i);
        }
    }

    m = new vec*[maxClusters];
    sigma = new mat*[maxClusters];
    A = new mat**[maxClusters];
    b = new vec**[maxClusters];
    sumOfSquares = new double*[maxClusters];
    E = new double[maxClusters];
    bestAxes = new int[maxClusters];
    coeffs = new vec*[maxClusters];
    probabilities = new double[maxClusters];
    for (int i = 0; i < maxClusters; ++i) {
        m[i] = NULL;
        sigma[i] = NULL;
        A[i] = NULL;
        b[i] = NULL;
        sumOfSquares[i] = NULL;
        coeffs[i] = NULL;
    }
    for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) {
        m[*it] = new vec(dim, fill::zeros);
        sigma[*it] = new mat(dim, dim, fill::zeros);
        A[*it] = new mat*[dim];
        b[*it] = new vec*[dim];
        sumOfSquares[*it] = new double[dim];
        E[*it] = 0.0;
        coeffs[*it] = new vec(dimA, fill::zeros);
        for (int i = 0; i < dim; ++i) {
            A[*it][i] = NULL;
            b[*it][i] = NULL;
        }
        for (int i = 0; i < dim; ++i) {
            A[*it][i] = new mat(dimA, dimA, fill::zeros);
            b[*it][i] = new vec(dimA, fill::zeros);
            sumOfSquares[*it][i] = 0.0;
        }
    }

    for (int i = 0; i < pointsNum; ++i) {
        int cl = labels[i];
        if (activeClusters.find(cl) != activeClusters.end()) {
            vec p = points.col(i);
            *m[cl] += p;
            for (int j = 0; j < dim; ++j) {
                vec v = values.submat(j * dimA, i, size(dimA, 1)); // Zamienic na liste!
                *A[cl][j] += v * v.t();
                *b[cl][j] += p[j] * v;
                sumOfSquares[cl][j] += p[j] * p[j];
            }
        }
    }
    for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) {
        *m[*it] /= card[*it];
        probabilities[*it] = ((double)card[*it]) / pointsNum;
    }

    for (int i = 0; i < pointsNum; ++i) {
        int cl = labels[i];
        if (activeClusters.find(cl) != activeClusters.end())
            *sigma[cl] += (points.col(i) - *m[cl]) * (points.col(i) - *m[cl]).t();
    }
    for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end();) {
        *sigma[*it] /= card[*it];
        double EMin = INFINITY;
        int bestAxis;
        double bestVar;
        vec bestCoeff;
        try {
            for (int i = 0; i < dim; ++i) {
                vec x = solve(*A[*it][i], *b[*it][i]);
                double var = sumOfSquares[*it][i];
                var -= ((vec)(2.0 * x.t() * *b[*it][i]))[0];
                var += ((mat)(x.t() * *A[*it][i] * x)).at(0, 0);
                var /= card[*it];

                mat sigmaTmp = *sigma[*it];
                sigmaTmp(i, span::all).fill(0.0);
                sigmaTmp(span::all, i).fill(0.0);
                sigmaTmp(i, i) = var;

                double prob = probabilities[*it];
                double ETmp = prob * (-log(prob) + (0.5 * ((dim * log(2.0 * M_PI * M_E)) + log(det(sigmaTmp)))));
                if (!is_finite(ETmp)) throw 0;
                if (ETmp < EMin) {
                    EMin = ETmp;
                    bestAxis = i;
                    bestVar = var;
                    bestCoeff = x;
                }
            }
            E[*it] = EMin;
            (*sigma[*it])(bestAxis, span::all).fill(0.0);
            (*sigma[*it])(span::all, bestAxis).fill(0.0);
            (*sigma[*it])(bestAxis, bestAxis) = bestVar;
            bestAxes[*it] = bestAxis;
            *coeffs[*it] = bestCoeff;
            ++it;
        } catch (...) {
            std::set<int>::iterator itTmp = it;
            ++itTmp;
            inactiveClusters.insert(*it);
            activeClusters.erase(it);
            it = itTmp;
        }
    }

    //   -*-   -*-   -*-

    for (std::set<int>::iterator it1 = inactiveClusters.begin(); it1 != inactiveClusters.end(); ++it1) {
        for (int i = 0; i < pointsNum; ++i) {
            int cl = labels[i];
            if (*it1 == cl) {
                int bestCl;
                double EMin = INFINITY;
                for (std::set<int>::iterator it2 = activeClusters.begin(); it2 != activeClusters.end(); ++it2) {
                    for (int j = 0; j < dim; ++j) {
                        vec p = points.col(i) - *m[*it2];
                        p[bestAxes[*it2]] = points.col(i)(bestAxes[*it2]) - ((mat)(coeffs[*it2]->t() * values.submat(bestAxes[*it2] * dimA, i, size(dimA, 1))))(0, 0);
                        double prob = probabilities[*it2];
                        double ETmp = -log(prob) + (0.5 * ((dim * log(2.0 * M_PI)) + log(det(*sigma[*it2])) + ((mat)(p.t() * inv(*sigma[*it2]) * p))(0, 0)));
                        if (!is_finite(ETmp)) throw 0;
                        if (ETmp < EMin) {
                            EMin = ETmp;
                            bestCl = *it2;
                        }
                    }
                }
                labels[i] = bestCl;
                ++card[labels[i]];
            }
        }
    }
    inactiveClusters.clear();

    for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) {
        m[*it]->fill(0.0);
        sigma[*it]->fill(0.0);
        for (int i = 0; i < dim; ++i) {
            A[*it][i]->fill(0.0);
            b[*it][i]->fill(0.0);
            sumOfSquares[*it][i] = 0.0;
        }
    }

    for (int i = 0; i < maxClusters; ++i) card[i] = 0;
    for (int i = 0; i < pointsNum; ++i) ++card[labels[i]];

    for (int i = 0; i < pointsNum; ++i) {
        int cl = labels[i];
        if (activeClusters.find(cl) != activeClusters.end()) {
            vec p = points.col(i);
            *m[cl] += p;
            for (int j = 0; j < dim; ++j) {
                vec v = values.submat(j * dimA, i, size(dimA, 1)); // Zamienic na liste!
                *A[cl][j] += v * v.t();
                *b[cl][j] += p[j] * v;
                sumOfSquares[cl][j] += p[j] * p[j];
            }
        }
    }
    for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) {
        *m[*it] /= card[*it];
        probabilities[*it] = ((double)card[*it]) / pointsNum;
    }

    for (int i = 0; i < pointsNum; ++i) {
        int cl = labels[i];
        if (activeClusters.find(cl) != activeClusters.end())
            *sigma[cl] += (points.col(i) - *m[cl]) * (points.col(i) - *m[cl]).t();
    }
    for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end();) {
        *sigma[*it] /= card[*it];
        double EMin = INFINITY;
        int bestAxis;
        double bestVar;
        vec bestCoeff;

        for (int i = 0; i < dim; ++i) {
            vec x = solve(*A[*it][i], *b[*it][i]);
            double var = sumOfSquares[*it][i];
            var -= ((vec)(2.0 * x.t() * *b[*it][i]))[0];
            var += ((mat)(x.t() * *A[*it][i] * x)).at(0, 0);
            var /= card[*it];

            mat sigmaTmp = *sigma[*it];
            sigmaTmp(i, span::all).fill(0.0);
            sigmaTmp(span::all, i).fill(0.0);
            sigmaTmp(i, i) = var;

            double prob = probabilities[*it];
            double ETmp = prob * (-log(prob) + (0.5 * ((dim * log(2.0 * M_PI * M_E)) + log(det(sigmaTmp)))));
            if (!is_finite(ETmp)) throw 0;
            if (ETmp < EMin) {
                EMin = ETmp;
                bestAxis = i;
                bestVar = var;
                bestCoeff = x;
            }
        }
        E[*it] = EMin;
        (*sigma[*it])(bestAxis, span::all).fill(0.0);
        (*sigma[*it])(span::all, bestAxis).fill(0.0);
        (*sigma[*it])(bestAxis, bestAxis) = bestVar;
        bestAxes[*it] = bestAxis;
        *coeffs[*it] = bestCoeff;
        ++it;
    }

    //   -*-   -*-   -*-

    double ETotalOld;
    double ETotalNew = 0.0;
    for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) ETotalNew += E[*it];
    int numberOfIterations = 0;
    do {
        ETotalOld = ETotalNew;

        for (int i = 0; i < pointsNum; ++i) {
            int bestCl;
            double EMin = INFINITY;
            for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) {
                vec p = points.col(i) - *m[*it];
                p[bestAxes[*it]] = points.col(i)(bestAxes[*it]) - ((mat)(coeffs[*it]->t() * values.submat(bestAxes[*it] * dimA, i, size(dimA, 1))))(0, 0);
                double prob = probabilities[*it];
                double ETmp = -log(prob) + (0.5 * ((dim * log(2.0 * M_PI)) + log(det(*sigma[*it])) + ((mat)(p.t() * inv(*sigma[*it]) * p))(0, 0)));
                if (!is_finite(ETmp)) throw 0;
                if (ETmp < EMin) {
                    EMin = ETmp;
                    bestCl = *it;
                }
            }
            labels[i] = bestCl;
        }

        activeClusters.clear();
        inactiveClusters.clear();
        for (int i = 0; i < maxClusters; ++i) card[i] = 0;
        for (int i = 0; i < pointsNum; ++i) ++card[labels[i]];
        for (int i = 0; i < maxClusters; ++i) {
            if ((card[i] >= dim + 1) && (((double)card[i]) / pointsNum >= cardMin)) activeClusters.insert(activeClusters.end(), i);
            else {
                if (card[i] > 0) inactiveClusters.insert(inactiveClusters.end(), i);
            }
        }

        // !!! REMOVAL !!!
        for (std::set<int>::iterator it1 = inactiveClusters.begin(); it1 != inactiveClusters.end(); ++it1) {
            for (int i = 0; i < pointsNum; ++i) {
                int cl = labels[i];
                if (*it1 == cl) {
                    int bestCl;
                    double EMin = INFINITY;
                    for (std::set<int>::iterator it2 = activeClusters.begin(); it2 != activeClusters.end(); ++it2) {
                        for (int j = 0; j < dim; ++j) {
                            vec p = points.col(i) - *m[*it2];
                            p[bestAxes[*it2]] = points.col(i)(bestAxes[*it2]) - ((mat)(coeffs[*it2]->t() * values.submat(bestAxes[*it2] * dimA, i, size(dimA, 1))))(0, 0);
                            double prob = probabilities[*it2];
                            double ETmp = -log(prob) + (0.5 * ((dim * log(2.0 * M_PI)) + log(det(*sigma[*it2])) + ((mat)(p.t() * inv(*sigma[*it2]) * p))(0, 0)));
                            if (!is_finite(ETmp)) throw 0;
                            if (ETmp < EMin) {
                                EMin = ETmp;
                                bestCl = *it2;
                            }
                        }
                    }
                    labels[i] = bestCl;
                    ++card[labels[i]];
                }
            }
            card[*it1] = 0;
        }
        inactiveClusters.clear();
        // !!! REMOVAL !!!

        // !!!
        for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) {
            m[*it]->fill(0.0);
            sigma[*it]->fill(0.0);
            for (int i = 0; i < dim; ++i) {
                A[*it][i]->fill(0.0);
                b[*it][i]->fill(0.0);
                sumOfSquares[*it][i] = 0.0;
            }
        }
        // !!!

        for (int i = 0; i < pointsNum; ++i) {
            int cl = labels[i];
            if (activeClusters.find(cl) != activeClusters.end()) {
                vec p = points.col(i);
                *m[cl] += p;
                for (int j = 0; j < dim; ++j) {
                    vec v = values.submat(j * dimA, i, size(dimA, 1)); // Zamienic na liste!
                    *A[cl][j] += v * v.t();
                    *b[cl][j] += p[j] * v;
                    sumOfSquares[cl][j] += p[j] * p[j];
                }
            }
        }
        for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) {
            *m[*it] /= card[*it];
            probabilities[*it] = ((double)card[*it]) / pointsNum;
        }

        for (int i = 0; i < pointsNum; ++i) {
            int cl = labels[i];
            if (activeClusters.find(cl) != activeClusters.end())
                *sigma[cl] += (points.col(i) - *m[cl]) * (points.col(i) - *m[cl]).t();
        }
        for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end();) {
            *sigma[*it] /= card[*it];
            double EMin = INFINITY;
            int bestAxis;
            double bestVar;
            vec bestCoeff;
            try {
                for (int i = 0; i < dim; ++i) {
                    vec x = solve(*A[*it][i], *b[*it][i]);
                    double var = sumOfSquares[*it][i];
                    var -= ((vec)(2.0 * x.t() * *b[*it][i]))[0];
                    var += ((mat)(x.t() * *A[*it][i] * x)).at(0, 0);
                    var /= card[*it];

                    mat sigmaTmp = *sigma[*it];
                    sigmaTmp(i, span::all).fill(0.0);
                    sigmaTmp(span::all, i).fill(0.0);
                    sigmaTmp(i, i) = var;

                    double prob = ((double)card[*it]) / pointsNum;
                    double ETmp = prob * (-log(prob) + (0.5 * ((dim * log(2.0 * M_PI * M_E)) + log(det(sigmaTmp)))));
                    if (!is_finite(ETmp)) throw 0;
                    if (ETmp < EMin) {
                        EMin = ETmp;
                        bestAxis = i;
                        bestVar = var;
                        bestCoeff = x;
                    }
                }
                E[*it] = EMin;
                (*sigma[*it])(bestAxis, span::all).fill(0.0);
                (*sigma[*it])(span::all, bestAxis).fill(0.0);
                (*sigma[*it])(bestAxis, bestAxis) = bestVar;
                bestAxes[*it] = bestAxis;
                *coeffs[*it] = bestCoeff;
                ++it;
            } catch (...) {
                std::set<int>::iterator itTmp = it;
                ++itTmp;
                inactiveClusters.insert(*it);
                activeClusters.erase(it);
                it = itTmp;
            }
        }

        ETotalNew = 0.0;
        for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) ETotalNew += E[*it];
        if (!is_finite(ETotalNew)) throw "A numerical overflow occurred.";

        if (interactive) {
            Rcpp::List cardL(maxClusters);
            Rcpp::List mL(maxClusters);
            Rcpp::List covL(maxClusters);
            Rcpp::List coeffsL(maxClusters);
            Rcpp::List dirsL(maxClusters);
            for (int i = 0; i < maxClusters; ++i) {
                if (activeClusters.find(i) != activeClusters.end()) {
                    cardL[i] = card[i];
                    mL[i] = *m[i];
                    covL[i] = *sigma[i];
                    coeffsL[i] = *coeffs[i];
                    dirsL[i] = bestAxes[i];
                } else {
                    cardL[i] = 0;
                    mL[i] = R_NilValue;
                    covL[i] = R_NilValue;
                    coeffsL[i] = R_NilValue;
                    dirsL[i] = R_NilValue;
                }
            }
            res.push_back(Rcpp::List::create(
                Rcpp::Named("cardinalities", cardL),
                Rcpp::Named("number_of_clusters", (int)activeClusters.size()),
                Rcpp::Named("labels", labels),
                Rcpp::Named("means", mL),
                Rcpp::Named("covariances", covL),
                Rcpp::Named("coefficients", coeffsL),
                Rcpp::Named("directions", dirsL),
                Rcpp::Named("cost_total", ETotalNew)
            ));
        }

        ++numberOfIterations;
    } while ((numberOfIterations < minIterations) || ((numberOfIterations < maxIterations) && (ETotalNew - ETotalOld < costThreshold)));

    //   -*-   -*-   -*-

    if (!interactive) {
        Rcpp::List cardL(maxClusters);
        Rcpp::List mL(maxClusters);
        Rcpp::List covL(maxClusters);
        Rcpp::List coeffsL(maxClusters);
        Rcpp::List dirsL(maxClusters);
        for (int i = 0; i < maxClusters; ++i) {
            if (activeClusters.find(i) != activeClusters.end()) {
                cardL[i] = card[i];
                mL[i] = *m[i];
                covL[i] = *sigma[i];
                coeffsL[i] = *coeffs[i];
                dirsL[i] = bestAxes[i];
            } else {
                cardL[i] = 0;
                mL[i] = R_NilValue;
                covL[i] = R_NilValue;
                coeffsL[i] = R_NilValue;
                dirsL[i] = R_NilValue;
            }
        }
        res = Rcpp::List::create(
            Rcpp::Named("cardinalities", cardL),
            Rcpp::Named("number_of_clusters", (int)activeClusters.size()),
            Rcpp::Named("labels", labels),
            Rcpp::Named("means", mL),
            Rcpp::Named("covariances", covL),
            Rcpp::Named("coefficients", coeffsL),
            Rcpp::Named("directions", dirsL),
            Rcpp::Named("cost_total", ETotalNew)
        );
    }
    return res;
}

//   -*-   -*-   -*-

void afCECHartiganRemoveCluster(
    int clToRemove,
    int pointsNum,
    ivec &labels,
    std::set<int> &activeClusters,
    const mat &points,
    int dim,
    vec **m,
    int *card,
    mat ***sigmaL,
    mat ***AL,
    const mat &values,
    int dimA,
    vec ***b,
    mat ***A,
    double **sumOfSquares,
    double *E,
    int *bestAxes,
    vec **coeffs,
    double *variances
) {
    for (int j = 0; j < pointsNum; ++j) {
        int cl = labels[j];
        if (clToRemove == cl) {
            double dEMinIncl = INFINITY;
            int bestClIncl;
            int bestAxisIncl;
            double EMinIncl;
            vec bestCoeffsIncl;
            double bestVarIncl;

            for (std::set<int>::iterator it2 = activeClusters.begin(); it2 != activeClusters.end(); ++it2) {
                vec v = (points.submat(1, j, size(dim - 1, 1)) - m[*it2]->subvec(1, dim - 1)) / sqrt(((double)(card[*it2] + 1.0)));
                for (int k = 0; k < dim; ++k) {
                    mat sigmaLTmp = *sigmaL[*it2][k];
                    vec vTmp = v;
                    CholeskyRankOneUpdate(sigmaLTmp, vTmp);
                    sigmaLTmp *= sqrt(card[*it2] / (card[*it2] + 1.0));
                    double det = 1.0;
                    for (int l = 0; l < dim - 1; ++l) det *= sigmaLTmp.at(l, l);

                    mat ALTmp = *AL[*it2][k];
                    vec vA = values.submat(k * dimA, j, size(dimA, 1));
                    vec vATmp = vA;
                    CholeskyRankOneUpdate(ALTmp, vATmp);
                    vec bTmp = *b[*it2][k] + (points.at(k, j) * vA);
                    vec y = solve(trimatl(ALTmp), bTmp);
                    vec x = solve(trimatu(ALTmp.t()), y);
                    mat ATmp = *A[*it2][k] + (vA * vA.t());
                    double var = sumOfSquares[*it2][k] + (points.at(k, j) * points.at(k, j));
                    var -= ((vec)(2.0 * x.t() * bTmp))[0];
                    var += ((mat)(x.t() * ATmp * x)).at(0, 0);
                    var /= (card[*it2] + 1);

                    double prob = (card[*it2] + 1.0) / pointsNum;
                    double EIncl = prob * (-log(prob) + (0.5 * ((dim * log(2.0 * M_PI * M_E)) + log(det * det * var))));
                    if (!is_finite(EIncl)) throw "A numerical overflow occurred.";
                    if (EIncl - E[*it2] < dEMinIncl) {
                        dEMinIncl = EIncl - E[*it2];
                        bestClIncl = *it2;
                        bestAxisIncl = k;
                        EMinIncl = EIncl;
                        bestCoeffsIncl = x;
                        bestVarIncl = var;
                    }

                    if (k < dim - 1) v[k] = (points.at(k, j) - (*m[*it2])[k]) / sqrt(((double)(card[*it2] + 1.0)));
                }
            }
            if (dEMinIncl < INFINITY) {
                vec v = (points.submat(1, j, size(dim - 1, 1)) - m[bestClIncl]->subvec(1, dim - 1)) / sqrt(((double)(card[bestClIncl] + 1.0)));
                for (int k = 0; k < dim; ++k) {
                    vec vTmp = v;
                    CholeskyRankOneUpdate(*sigmaL[bestClIncl][k], vTmp);
                    *sigmaL[bestClIncl][k] *= sqrt(card[bestClIncl] / (card[bestClIncl] + 1.0));

                    vec vA = values.submat(k * dimA, j, size(dimA, 1));
                    vec vATmp = vA;
                    CholeskyRankOneUpdate(*AL[bestClIncl][k], vATmp);
                    *b[bestClIncl][k] += points.at(k, j) * vA;
                    *A[bestClIncl][k] += (vA * vA.t());

                    sumOfSquares[bestClIncl][k] += points.at(k, j) * points.at(k, j);

                    if (k < dim - 1) v[k] = (points.at(k, j) - (*m[bestClIncl])[k]) / sqrt(((double)(card[bestClIncl] + 1.0)));
                }
                *m[bestClIncl] = ((*m[bestClIncl] * card[bestClIncl]) + points.col(j)) / (card[bestClIncl] + 1.0);
                bestAxes[bestClIncl] = bestAxisIncl;
                E[bestClIncl] = EMinIncl;
                ++card[bestClIncl];
                labels[j] = bestClIncl;
                *coeffs[bestClIncl] = bestCoeffsIncl;
                variances[bestClIncl] = bestVarIncl;
            }
        }
    }
}

//   -*-   -*-   -*-

Rcpp::List afCECHartigan (
    const mat &points,
    int maxClusters,
    ivec labels,
    double cardMin,
    double costThreshold,
    int minIterations,
    int maxIterations,
    const mat &values,
    bool interactive
) {
    int *card = NULL;
    vec **m = NULL;
    mat **sigma = NULL;
    mat ***sigmaL = NULL;
    mat ***A = NULL;
    mat ***AL = NULL;
    vec ***b = NULL;
    double **sumOfSquares = NULL;
    double *E = NULL;
    int *bestAxes = NULL;
    vec **coeffs = NULL;
    double *variances = NULL;

    mat *badCl; // !!! !!! !!!
    int varErr = 0; // !!! !!! !!!

    int pointsNum = points.n_cols;
    int dim = points.n_rows;
    int dimA = values.n_rows / dim;

    Rcpp::List res;

    card = new int[maxClusters];
    std::set<int> activeClusters;
    std::set<int> inactiveClusters;
    for (int i = 0; i < maxClusters; ++i) card[i] = 0;
    for (int i = 0; i < pointsNum; ++i) ++card[labels[i]];
    for (int i = 0; i < maxClusters; ++i) {
        if ((card[i] >= dim + 1) && (((double)card[i]) / pointsNum >= cardMin)) activeClusters.insert(activeClusters.end(), i);
    }
    if (activeClusters.size() == 0) {
        for (int i = 0; i < pointsNum; ++i) labels[i] = 0;
        activeClusters.insert(activeClusters.end(), 0);
        card[0] = pointsNum;
    }

    m = new vec*[maxClusters];
    sigma = new mat*[maxClusters];
    sigmaL = new mat**[maxClusters];
    A = new mat**[maxClusters];
    AL = new mat**[maxClusters];
    b = new vec**[maxClusters];
    sumOfSquares = new double*[maxClusters];
    E = new double[maxClusters];
    bestAxes = new int[maxClusters];
    coeffs = new vec*[maxClusters];
    variances = new double[maxClusters];
    for (int i = 0; i < maxClusters; ++i) {
        m[i] = NULL;
        sigma[i] = NULL;
        sigmaL[i] = NULL;
        A[i] = NULL;
        AL[i] = NULL;
        b[i] = NULL;
        sumOfSquares[i] = NULL;
        coeffs[i] = NULL;
        variances[i] = 0.0;
    }
    for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) {
        m[*it] = new vec(dim, fill::zeros);
        sigma[*it] = new mat(dim, dim, fill::zeros);
        sigmaL[*it] = new mat*[dim];
        A[*it] = new mat*[dim];
        AL[*it] = new mat*[dim];
        b[*it] = new vec*[dim];
        sumOfSquares[*it] = new double[dim];
        E[*it] = 0.0;
        coeffs[*it] = new vec(dimA, fill::zeros);
        for (int i = 0; i < dim; ++i) {
            sigmaL[*it][i] = NULL;
            A[*it][i] = NULL;
            AL[*it][i] = NULL;
            b[*it][i] = NULL;
        }
        for (int i = 0; i < dim; ++i) {
            sigmaL[*it][i] = new mat(dim - 1, dim - 1);
            A[*it][i] = new mat(dimA, dimA, fill::zeros);
            AL[*it][i] = new mat(dimA, dimA);
            b[*it][i] = new vec(dimA, fill::zeros);
            sumOfSquares[*it][i] = 0.0;
        }
    }

    for (int i = 0; i < pointsNum; ++i) {
        int cl = labels[i];
        if (activeClusters.find(cl) != activeClusters.end()) {
            vec p = points.col(i);
            *m[cl] += p;
            for (int j = 0; j < dim; ++j) {
                vec v = values.submat(j * dimA, i, size(dimA, 1)); // Zamienic na liste!
                *A[cl][j] += v * v.t();
                *b[cl][j] += p[j] * v;
                sumOfSquares[cl][j] += p[j] * p[j];
            }
        }
    }
    for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) *m[*it] /= card[*it];

    for (int i = 0; i < pointsNum; ++i) {
        int cl = labels[i];
        if (activeClusters.find(cl) != activeClusters.end())
            *sigma[cl] += (points.col(i) - *m[cl]) * (points.col(i) - *m[cl]).t();
    }
    for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end();) {
        *sigma[*it] /= card[*it];
        uvec ind(dim - 1);
        for (int i = 0; i < dim - 1; ++i) ind[i] = i + 1;
        double EMin = INFINITY;
        int bestAxis;
        double bestVar;
        vec bestCoeff;
        try {
            for (int i = 0; i < dim; ++i) {
                *sigmaL[*it][i] = chol(sigma[*it]->submat(ind, ind), "lower");
                double det = 1.0;
                for (int j = 0; j < dim - 1; ++j) det *= sigmaL[*it][i]->at(j, j);

                *AL[*it][i] = chol(*A[*it][i], "lower");
                vec y = solve(trimatl(*AL[*it][i]), *b[*it][i]);
                vec x = solve(trimatu(AL[*it][i]->t()), y);
                double var = sumOfSquares[*it][i];
                var -= ((vec)(2.0 * x.t() * *b[*it][i]))[0];
                var += ((mat)(x.t() * *A[*it][i] * x)).at(0, 0);
                var /= card[*it];

                double prob = ((double)card[*it]) / pointsNum;
                double ETmp = prob * (-log(prob) + (0.5 * ((dim * log(2.0 * M_PI * M_E)) + log(det * det * var))));
                if (!is_finite(ETmp)) throw 0;
                if (ETmp < EMin) {
                    EMin = ETmp;
                    bestAxis = i; // było: *it;
                    bestCoeff = x;
                    bestVar = var;
                }

                if (i < dim - 1) ind[i] = i;
            }
            E[*it] = EMin;
            bestAxes[*it] = bestAxis;
            *coeffs[*it] = bestCoeff;
            variances[*it] = bestVar;
            ++it;
        } catch (...) {
            std::set<int>::iterator itTmp = it;
            ++itTmp;
            inactiveClusters.insert(*it);
            activeClusters.erase(it);
            it = itTmp;
        }
    }

    //   -*-   -*-   -*-

    //printf("---------------------------------\n");

    double ETotalOld;
    double ETotalNew = 0.0; // !!! !!! !!!
    for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) ETotalNew += E[*it];
    int numberOfIterations = 0;
    do {
        //printf("%d\n", numberOfIterations);

        ETotalOld = ETotalNew;
        for (int i = 0; i < pointsNum; ++i) {
            //printf("%d\n", i);

            int cl = labels[i];
            if (activeClusters.find(cl) != activeClusters.end()) {
                double dEMinExcl = INFINITY;
                int bestAxisExcl;
                double EMinExcl;
                vec bestCoeffsExcl;
                double bestVarExcl;
                try {
                    vec v = (points.submat(1, i, size(dim - 1, 1)) - m[cl]->subvec(1, dim - 1)) / sqrt(((double)(card[cl] - 1.0)));
                    for (int j = 0; j < dim; ++j) {
                        mat sigmaLTmp = *sigmaL[cl][j];
                        vec vTmp = v;
                        CholeskyRankOneDowndate(sigmaLTmp, vTmp);
                        sigmaLTmp *= sqrt(card[cl] / (card[cl] - 1.0));
                        double det = 1.0;
                        for (int k = 0; k < dim - 1; ++k) det *= sigmaLTmp.at(k, k);

                        mat ALTmp = *AL[cl][j];
                        vec vA = values.submat(j * dimA, i, size(dimA, 1));
                        vec vATmp = vA;
                        CholeskyRankOneDowndate(ALTmp, vATmp);
                        vec bTmp = *b[cl][j] - (points.at(j, i) * vA);
                        vec y = solve(trimatl(ALTmp), bTmp);
                        vec x = solve(trimatu(ALTmp.t()), y);
                        mat ATmp = *A[cl][j] - (vA * vA.t());
                        double var = sumOfSquares[cl][j] - (points.at(j, i) * points.at(j, i));
                        var -= ((vec)(2.0 * x.t() * bTmp))[0];
                        var += ((mat)(x.t() * ATmp * x)).at(0, 0);
                        var /= (card[cl] - 1);

                        double prob = (card[cl] - 1.0) / pointsNum;
                        double EExcl = prob * (-log(prob) + (0.5 * ((dim * log(2.0 * M_PI * M_E)) + log(det * det * var))));
                        if (!is_finite(EExcl)) throw 0;
                        if (EExcl - E[cl] < dEMinExcl) {
                            dEMinExcl = EExcl - E[cl];
                            bestAxisExcl = j;
                            EMinExcl = EExcl;
                            bestCoeffsExcl = x;
                            bestVarExcl = var;
                        }

                        if (j < dim - 1) v[j] = (points.at(j, i) - (*m[cl])[j]) / sqrt(((double)(card[cl] - 1.0)));
                    }
                } catch (...) {
                    //printf("WYWALAMY %d\n", cl);
                    activeClusters.erase(activeClusters.find(cl));
                    afCECHartiganRemoveCluster(
                        cl,
                        pointsNum,
                        labels,
                        activeClusters,
                        points,
                        dim,
                        m,
                        card,
                        sigmaL,
                        AL,
                        values,
                        dimA,
                        b,
                        A,
                        sumOfSquares,
                        E,
                        bestAxes,
                        coeffs,
                        variances
                    );
                }
                if (activeClusters.find(cl) != activeClusters.end()) {
                    double dEMinIncl = INFINITY;
                    int bestClIncl;
                    int bestAxisIncl;
                    double EMinIncl;
                    vec bestCoeffsIncl;
                    double bestVarIncl;
                    for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) {
                        if ((*it != cl) && (activeClusters.find(*it) != activeClusters.end())) {
                            vec v = (points.submat(1, i, size(dim - 1, 1)) - m[*it]->subvec(1, dim - 1)) / sqrt(((double)(card[*it] + 1.0)));
                            for (int j = 0; j < dim; ++j) {
                                mat sigmaLTmp = *sigmaL[*it][j];
                                vec vTmp = v;
                                CholeskyRankOneUpdate(sigmaLTmp, vTmp);
                                sigmaLTmp *= sqrt(card[*it] / (card[*it] + 1.0));
                                double det = 1.0;
                                for (int k = 0; k < dim - 1; ++k) det *= sigmaLTmp.at(k, k);

                                mat ALTmp = *AL[*it][j];
                                vec vA = values.submat(j * dimA, i, size(dimA, 1));
                                vec vATmp = vA;
                                CholeskyRankOneUpdate(ALTmp, vATmp);
                                vec bTmp = *b[*it][j] + (points.at(j, i) * vA);
                                vec y = solve(trimatl(ALTmp), bTmp);
                                vec x = solve(trimatu(ALTmp.t()), y);
                                mat ATmp = *A[*it][j] + (vA * vA.t());
                                double var = sumOfSquares[*it][j] + (points.at(j, i) * points.at(j, i));
                                var -= ((vec)(2.0 * x.t() * bTmp))[0];
                                var += ((mat)(x.t() * ATmp * x)).at(0, 0);
                                var /= (card[*it] + 1);

                                double prob = (card[*it] + 1.0) / pointsNum;
                                double EIncl = prob * (-log(prob) + (0.5 * ((dim * log(2.0 * M_PI * M_E)) + log(det * det * var))));
                                if (!is_finite(EIncl)) throw "A numerical overflow occurred.";
                                if (EIncl - E[*it] < dEMinIncl) {
                                    dEMinIncl = EIncl - E[*it];
                                    bestClIncl = *it;
                                    bestAxisIncl = j;
                                    EMinIncl = EIncl;
                                    bestCoeffsIncl = x;
                                    bestVarIncl = var;
                                }

                                if (j < dim - 1) v[j] = (points.at(j, i) - (*m[*it])[j]) / sqrt(((double)(card[*it] + 1.0)));
                            }

                        }
                    }
                    if (dEMinExcl + dEMinIncl < 0.0) {
                        vec v = (points.submat(1, i, size(dim - 1, 1)) - m[cl]->subvec(1, dim - 1)) / sqrt(((double)(card[cl] - 1.0)));
                        for (int j = 0; j < dim; ++j) {
                            vec vTmp = v;
                            CholeskyRankOneDowndate(*sigmaL[cl][j], vTmp);
                            *sigmaL[cl][j] *= sqrt(card[cl] / (card[cl] - 1.0));

                            vec vA = values.submat(j * dimA, i, size(dimA, 1));
                            vec vATmp = vA;
                            CholeskyRankOneDowndate(*AL[cl][j], vATmp);
                            *b[cl][j] -= points.at(j, i) * vA;
                            *A[cl][j] -= (vA * vA.t());

                            sumOfSquares[cl][j] -= points.at(j, i) * points.at(j, i);

                            if (j < dim - 1) v[j] = (points.at(j, i) - (*m[cl])[j]) / sqrt(((double)(card[cl] - 1.0)));
                        }
                        v = (points.submat(1, i, size(dim - 1, 1)) - m[bestClIncl]->subvec(1, dim - 1)) / sqrt(((double)(card[bestClIncl] + 1.0)));
                        for (int j = 0; j < dim; ++j) {
                            vec vTmp = v;
                            CholeskyRankOneUpdate(*sigmaL[bestClIncl][j], vTmp);
                            *sigmaL[bestClIncl][j] *= sqrt(card[bestClIncl] / (card[bestClIncl] + 1.0));

                            vec vA = values.submat(j * dimA, i, size(dimA, 1));
                            vec vATmp = vA;
                            CholeskyRankOneUpdate(*AL[bestClIncl][j], vATmp);
                            *b[bestClIncl][j] += points.at(j, i) * vA;
                            *A[bestClIncl][j] += (vA * vA.t());

                            sumOfSquares[bestClIncl][j] += points.at(j, i) * points.at(j, i); //Było: sumOfSquares[cl][j] += points.at(j, i) * points.at(j, i);

                            if (j < dim - 1) v[j] = (points.at(j, i) - (*m[bestClIncl])[j]) / sqrt(((double)(card[bestClIncl] + 1.0)));
                        }
                        *m[cl] = ((*m[cl] * card[cl]) - points.col(i)) / (card[cl] - 1.0);
                        *m[bestClIncl] = ((*m[bestClIncl] * card[bestClIncl]) + points.col(i)) / (card[bestClIncl] + 1.0);
                        bestAxes[cl] = bestAxisExcl;
                        bestAxes[bestClIncl] = bestAxisIncl;
                        E[cl] = EMinExcl;
                        E[bestClIncl] = EMinIncl;
                        --card[cl];
                        ++card[bestClIncl];
                        labels[i] = bestClIncl;
                        variances[cl] = bestVarExcl;
                        variances[bestClIncl] = bestVarIncl;
                        *coeffs[cl] = bestCoeffsExcl;
                        *coeffs[bestClIncl] = bestCoeffsIncl;
                        if (((double)card[cl]) / pointsNum < cardMin) {
                            //printf("WYWALAMY %d\n", cl);
                            activeClusters.erase(activeClusters.find(cl));
                            afCECHartiganRemoveCluster(
                                cl,
                                pointsNum,
                                labels,
                                activeClusters,
                                points,
                                dim,
                                m,
                                card,
                                sigmaL,
                                AL,
                                values,
                                dimA,
                                b,
                                A,
                                sumOfSquares,
                                E,
                                bestAxes,
                                coeffs,
                                variances
                            );
                        }
                    }
                }
            } else {
                //printf("WYWALAMY %d\n", cl);
                afCECHartiganRemoveCluster(
                    cl,
                    pointsNum,
                    labels,
                    activeClusters,
                    points,
                    dim,
                    m,
                    card,
                    sigmaL,
                    AL,
                    values,
                    dimA,
                    b,
                    A,
                    sumOfSquares,
                    E,
                    bestAxes,
                    coeffs,
                    variances
                );
            }
        }
        ETotalNew = 0.0;
        for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) ETotalNew += E[*it];
        if (!is_finite(ETotalNew)) throw "A numerical overflow occurred.";

        if (interactive) {
            Rcpp::List cardL(maxClusters);
            Rcpp::List mL(maxClusters);
            Rcpp::List covL(maxClusters);
            Rcpp::List GLMDesignMatrixL(maxClusters); // !!!
            Rcpp::List GLMResponseVectorL(maxClusters); // !!!
            Rcpp::List coeffsL(maxClusters);
            Rcpp::List dirsL(maxClusters);
            Rcpp::List EL(maxClusters); // !!!
            for (int i = 0; i < maxClusters; ++i) {
                if (activeClusters.find(i) != activeClusters.end()) {
                    uvec ind(dim - 1);
                    for (int j = 0; j < bestAxes[i]; ++j) ind(j) = j;
                    for (int j = bestAxes[i] + 1; j < dim; ++j) ind(j - 1) = j;

                    (*sigma[i])(bestAxes[i], span::all).fill(0.0);
                    (*sigma[i])(span::all, bestAxes[i]).fill(0.0);
                    (*sigma[i])(bestAxes[i], bestAxes[i]) = variances[i];
                    sigma[i]->submat(ind, ind) = *sigmaL[i][bestAxes[i]] * sigmaL[i][bestAxes[i]]->t();

                    cardL[i] = card[i];
                    mL[i] = *m[i];
                    covL[i] = *sigma[i];

                    // !!!
                    Rcpp::List tmp1(dim);
                    Rcpp::List tmp2(dim);
                    for (int j = 0; j < dim; ++j) {
                        tmp1[j] = *A[i][j];
                        tmp2[j] = *b[i][j];
                    }
                    GLMDesignMatrixL[i] = tmp1;
                    GLMResponseVectorL[i] = tmp2;
                    // !!!

                    coeffsL[i] = *coeffs[i];
                    dirsL[i] = bestAxes[i];
                    EL[i] = E[i]; // !!!
                } else {
                    cardL[i] = 0;
                    mL[i] = R_NilValue;
                    covL[i] = R_NilValue;
                    GLMDesignMatrixL[i] = R_NilValue; // !!!
                    GLMResponseVectorL[i] = R_NilValue; // !!!
                    coeffsL[i] = R_NilValue;
                    dirsL[i] = R_NilValue;
                    EL[i] = R_NilValue; // !!!
                }
            }
            res.push_back(Rcpp::List::create(
                Rcpp::Named("cardinalities", cardL),
                Rcpp::Named("number_of_clusters", (int)activeClusters.size()),
                Rcpp::Named("labels", labels),
                Rcpp::Named("means", mL),
                Rcpp::Named("covariances", covL),
                Rcpp::Named("GLM_design_matrix", GLMDesignMatrixL), // !!!
                Rcpp::Named("GLM_response_vector", GLMResponseVectorL), // !!!
                Rcpp::Named("coefficients", coeffsL),
                Rcpp::Named("directions", dirsL),
                Rcpp::Named("cost_total", ETotalNew),
                Rcpp::Named("cost", EL) // !!!
            ));
        }

        ++numberOfIterations;
    } while ((numberOfIterations < minIterations) || ((numberOfIterations < maxIterations) && (ETotalNew - ETotalOld < costThreshold)));

    //   -*-   -*-   -*-

    if (!interactive) {
        Rcpp::List cardL(maxClusters);
        Rcpp::List mL(maxClusters);
        Rcpp::List covL(maxClusters);
        Rcpp::List GLMDesignMatrixL(maxClusters); // !!!
        Rcpp::List GLMResponseVectorL(maxClusters); // !!!
        Rcpp::List coeffsL(maxClusters);
        Rcpp::List dirsL(maxClusters);
        Rcpp::List EL(maxClusters); // !!!
        for (int i = 0; i < maxClusters; ++i) {
            if (activeClusters.find(i) != activeClusters.end()) {
                uvec ind(dim - 1);
                for (int j = 0; j < bestAxes[i]; ++j) ind(j) = j;
                for (int j = bestAxes[i] + 1; j < dim; ++j) ind(j - 1) = j;

                (*sigma[i])(bestAxes[i], span::all).fill(0.0);
                (*sigma[i])(span::all, bestAxes[i]).fill(0.0);
                (*sigma[i])(bestAxes[i], bestAxes[i]) = variances[i];
                sigma[i]->submat(ind, ind) = *sigmaL[i][bestAxes[i]] * sigmaL[i][bestAxes[i]]->t();

                cardL[i] = card[i];
                mL[i] = *m[i];
                covL[i] = *sigma[i];

                // !!!
                Rcpp::List tmp1(dim);
                Rcpp::List tmp2(dim);
                for (int j = 0; j < dim; ++j) {
                    tmp1[j] = *A[i][j];
                    tmp2[j] = *b[i][j];
                }
                GLMDesignMatrixL[i] = tmp1;
                GLMResponseVectorL[i] = tmp2;
                // !!!

                coeffsL[i] = *coeffs[i];
                dirsL[i] = bestAxes[i];
                EL[i] = E[i]; // !!!
            } else {
                cardL[i] = 0;
                mL[i] = R_NilValue;
                covL[i] = R_NilValue;
                GLMDesignMatrixL[i] = R_NilValue; // !!!
                GLMResponseVectorL[i] = R_NilValue; // !!!
                coeffsL[i] = R_NilValue;
                dirsL[i] = R_NilValue;
                EL[i] = R_NilValue; // !!!
            }
        }
        res = Rcpp::List::create(
            Rcpp::Named("cardinalities", cardL),
            Rcpp::Named("number_of_clusters", (int)activeClusters.size()),
            Rcpp::Named("labels", labels),
            Rcpp::Named("means", mL),
            Rcpp::Named("covariances", covL),
            Rcpp::Named("GLM_design_matrix", GLMDesignMatrixL), // !!!
            Rcpp::Named("GLM_response_vector", GLMResponseVectorL), // !!!
            Rcpp::Named("coefficients", coeffsL),
            Rcpp::Named("directions", dirsL),
            Rcpp::Named("cost_total", ETotalNew),
            Rcpp::Named("cost", EL) // !!!
        );
    }
    return res;
}

//   -*-   -*-   -*-

// [[Rcpp::export]]
Rcpp::List afCECCppRoutine (
    const arma::mat &points,
    int maxClusters,
    const SEXP &initialLabels,
    double cardMin,
    double costThreshold,
    int minIterations,
    int maxIterations,
    int numberOfStarts,
    const std::string &method,
    const arma::mat &values,
    bool interactive
) try {
    int pointsNum = points.n_cols;
    int dim = points.n_rows;

    if (dim < 2) throw "The points have to be of the dimension at least 2.";

    if (maxClusters < 1) throw "The maximum number of clusters should be a positive value.";
    if (maxClusters > pointsNum) throw "The maximum number of clusters cannot exceed the number of points.";

    if (numberOfStarts < 1) throw "The number of starts should be a positive value.";

    std::string initialLabelsStr;
    imat initialLabelsMat;
    try {
        switch (TYPEOF(initialLabels)) {
            case STRSXP : {
                initialLabelsStr = Rcpp::as<std::string>(initialLabels);
                if ((initialLabelsStr != "random") && (initialLabelsStr != "k-means++"))
                    throw "The initialLabels parameter as a string directive should take a value of either \"random\" or \"k-means++\".";
                break;
            }
            case INTSXP : {
                if (numberOfStarts == 1) {
                    ivec initialLabelsVec = Rcpp::as<ivec>(initialLabels);
                    initialLabelsMat = imat(initialLabelsVec);
                } else
                    initialLabelsMat = Rcpp::as<imat>(initialLabels);
                if ((initialLabelsMat.n_rows != pointsNum) || (initialLabelsMat.n_cols != numberOfStarts))
                    throw "The initialLabels parameter as a matrix should be of the dimensions (number of points x number of starts).";
                for (int i = 0; i < initialLabelsMat.n_rows; ++i) {
                    for (int j = 0; j < initialLabelsMat.n_cols; ++j) {
                        if ((initialLabelsMat.at(i, j) < 0) || (initialLabelsMat.at(i, j) >= maxClusters))
                            throw "The initialLabels parameter as a matrix should have elements of the range [0, maxClusters - 1].";
                    }
                }
                break;
            }
            default : {
                throw "The initialLabels parameter should be either a string directive or a matrix of the type integer.";
            }
        }
    } catch (std::bad_alloc &e) {
        throw "Not enough memory.";
    } catch (std::exception &e) {
        throw "The initialLabels parameter should be either a string directive or a matrix of the type integer.";
    }

    if (maxIterations < 0)
        throw "The maximum number of iterations should be a non-negative value.";
    if ((method != "Lloyd") && (method != "Hartigan"))
        throw "The method parameter should take a value of either \"Lloyd\" or \"Hartigan\".";

    ivec labels(pointsNum);
    double EMin = INFINITY;
    Rcpp::List bestRes;
    for (int i = 0; i < numberOfStarts; ++i) {
        try {
            switch (TYPEOF(initialLabels)) {
                case STRSXP : {
                    if (initialLabelsStr == "random") {
                        //for (int j = 0; j < pointsNum; ++j) labels[j] = (RandomInteger() & 2147483647) % maxClusters;
                        for (int j = 0; j < pointsNum; ++j) labels[j] = (RandomInteger() & 2147483647) % maxClusters;
                    } else {
                        // We perform the k-means++ algorithm.
                        mat centers(dim, maxClusters);
                        vec cdf(pointsNum);
                        centers.col(0) = points.col((RandomInteger() & 2147483647) % pointsNum);
                        for (int j = 1; j < maxClusters; ++j) {
                            for (int k = 0; k < pointsNum; ++k) {
                                vec p = points.col(k);
                                double dMin = INFINITY;
                                for (int l = 0; l < j; ++l) {
                                    double d = norm(p - centers.col(l));
                                    if (d < dMin) dMin = d;
                                }
                                if (k == 0) cdf[k] = dMin * dMin;
                                else
                                    cdf[k] = cdf[k - 1] + (dMin * dMin);
                            }
                            cdf /= cdf[pointsNum - 1];
                            if (!cdf.is_finite()) throw "A numerical overflow occurred.";
                            // In order to apply the inverse transform sampling we have to generate an uniformly distributed
                            // pseudorandom value belonging to the interval (0, 1).
                            double value = RandomFloat();
                            // If the generated value is equal to 0, we generate it one more time. In the next step we are
                            // guaranteed to get a value different than 0.
                            if (value == 0.0) value = RandomFloat();
                            // We apply the inverse transform sampling.
                            int lB = 0;
                            int uB = pointsNum - 1;
                            while (uB > lB) {
                                int pivot = (lB + uB) >> 1;
                                if (value <= cdf[pivot]) uB = pivot;
                                else
                                    lB = pivot + 1;
                            }
                            centers.col(j) = points.col(lB);
                        }
                        for (int j = 0; j < pointsNum; ++j) {
                            vec p = points.col(j);
                            double dMin = INFINITY;
                            int dMinInd;
                            for (int k = 0; k < maxClusters; ++k) {
                                double d = norm(p - centers.col(k));
                                if (d < dMin) {
                                    dMin = d;
                                    dMinInd = k;
                                }
                            }
                            if (is_finite(dMin)) labels[j] = dMinInd;
                            else
                                throw "A numerical overflow occurred.";
                        }
                    }
                    break;
                }
                case INTSXP : {
                    labels = initialLabelsMat.col(i);
                    break;
                }
            }
            Rcpp::List res;
            if (method == "Lloyd") res = afCECLloyd(points, maxClusters, labels, cardMin, costThreshold, minIterations, maxIterations, values, interactive);
            else {
                res = afCECHartigan(points, maxClusters, labels, cardMin, costThreshold, minIterations, maxIterations, values, interactive);

                //CafCECHartigan afCECHartigan(points, maxClusters, labels, cardMin, minIterations, maxIterations, 0.0, values, interactive);
                //res = afCECHartigan.res;
            }
            if (res.length() > 0) {
                if (!interactive) {
                    double E = res["cost_total"];
                    if (E < EMin) {
                        EMin = E;
                        bestRes = res;
                    }
                } else {
                    Rcpp::List tmp = res[res.length() - 1];
                    double E = tmp["cost_total"];
                    if (E < EMin) {
                        EMin = E;
                        bestRes = res;
                    }
                }
            }
        } catch (std::bad_alloc &e) {
        } catch (const char *e) {
        } catch (...) {
        }
    }

    //if (EMin == INFINITY) throw "A clustering hasn't been found during any of the starts.";
    return bestRes;
} catch (const char *e) {
    return Rcpp::List();
} catch (...) {
    return Rcpp::List();
}
