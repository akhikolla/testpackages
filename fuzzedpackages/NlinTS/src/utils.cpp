#include<Rcpp.h>
#include "../inst/include/utils.h"

/******************************/
double get_random()
{
    static std::default_random_engine e;
    e. seed (5);
    static std::uniform_real_distribution<double> dis(0.0, 0.1); // rage 0 - 1
    return dis(e);
}

/******************************/
vector <double> random_vector (unsigned long size)
{
    vector<double> nums;
     for (unsigned long i = 0; i < size; ++i)
         nums. push_back (get_random() * 2.0 / sqrt (size));

     return nums;
}

/************ get column *************/
std::vector<double> get_col (const std::vector<std::vector<double>> & A, unsigned long j)
{
    std::vector<double> res;
    for (unsigned long i = 0; i < A. size (); ++i)
        res. push_back (A[i][j]);

    return res;
}

/*********** Transpose ***************/
std::vector<std::vector<double>> Transpose (const std::vector<std::vector<double>> & M)
{
    unsigned long n = M.size(), m = M[0].size ();
    std::vector<std::vector<double>> T (m);
    for (unsigned int i = 0 ; i < m ; ++i)
    {
        T[i] = std::vector<double> (n);
         for (unsigned long j = 0 ; j < n ; ++j)
             T[i][j] = M[j][i];
    }
        return T;
}

/********************************/
MatD reshape (const VectD & A, unsigned axis)
{
    MatD RES;

    if (axis == 0)
    {
        for (double a : A)
            RES.push_back ({a});
    }
    else if (axis == 1)
    {
        RES. resize (1);
        for (double a : A)
            RES[0].push_back (a);
    }
    else
    {
        Rcpp::Rcout << "Wring value for exis, choices in [0, 1] \n";
    }

    return RES;
}

/***************************************/
VectD matrix_dot (const VectD & A, const VectD & B)
{
    if (A. size () != B.size ())
    {
        Rcpp::Rcout << "Error when multiplying two vectors, they have not the same size. \n";
        Rcpp::stop ("\n.");
    }
    VectD res;

    res.resize (A.size());

    for (unsigned long i = 0 ; i < A.size() ; ++i)
        res[i]  = (A[i] * B[i]) ;

    return res;
}

/***************************************/
void matrix_dot (MatD & A, double a)
{
    for (unsigned long i = 0 ; i < A.size () ; ++i)
        for (unsigned long j = 0 ; j < A[0].size () ; ++j)
            A[i][j] *= a ;
}
/***************************************/
VectD matrix_dot (const MatD & A, const VectD & B)
{
    VectD res;
    res.resize(A.size ());

    for (unsigned long i = 0 ; i < A.size () ; ++i)
        for (unsigned long j = 0 ; j < B.size () ; ++j)
            res[i] += (A[i][j] * B[j]) ;

    return res;
}

/************************************************/
MatD matrix_sum (const MatD & A, const MatD & B, const MatD & C, const MatD & D)
{
    if (A. size () != B.size () or B. size () != C.size () or C. size () != D.size ())
    {
        Rcpp::Rcout << "Error when summing matrix, they have not the same length. \n";
        Rcpp::stop ("\n.");
    }
    if (A[0]. size () != B[0].size () or B[0]. size () != C[0].size () or C[0]. size () != D[0].size ())
    {
        Rcpp::Rcout << "Error when summing matrix, they have not the same number of columns. \n";
        Rcpp::stop ("\n.");
    }

    MatD RES (A.size ());

    for (unsigned long i = 0 ; i < A.size() ; ++i)
    {
        for (unsigned long j = 0 ; j < A[0].size() ; ++j)
            RES[i]. push_back (A[i][j] + B[i][j] + C[i][j] + D[i][j]);
     }

    return RES;
}

/************************************************/
VectD matrix_sum (const VectD & A, const VectD & B, const VectD & C, const VectD & D)
{
    if (A. size () != B.size () or B. size () != C.size () or C. size () != D.size ())
    {
        Rcpp::Rcout << "Error when summing 4 vectors, they have not the same length. \n";
        Rcpp::stop ("\n.");
    }

    VectD Res (A.size ());

    for (unsigned long i = 0 ; i < A.size() ; ++i)
    {
            Res[i] = A[i] + B[i] + C[i] + D[i];
     }

    return Res;
}

/********************************************/
double sum_vect (const vector<double> & Vect)
{
    if (0 == Vect.size ())
    {
       Rcpp::Rcout << ("Vector of size null");
       Rcpp::stop ("\n.");
    }
    double sum = 0;
     for (auto val: Vect)
          sum += val;

     return sum;
}

/***************************************/
MatD matrix_dot (const MatD & A, const MatD & B)
{
    MatD RES (A.size ());

        for (unsigned long i = 0 ; i < A.size () ; ++i)
        {
            RES[i]. resize (B[0].size ());
            for (unsigned long j = 0 ; j < B[0].size () ; ++j)
            {
                RES[i][j] = sum_vect (matrix_dot (A[i], get_col (B, j))) ;
            }
        }
    return RES;
}

/***************************************/
VectD matrix_sum (const VectD & A, const VectD & B)
{
    if (A. size () != B.size ())
    {
        Rcpp::Rcout << "Error when summing two vectors, they have not the same size. \n";
        Rcpp::stop ("\n.");
    }
    VectD res;

    res.resize (A.size());

    for (unsigned long i = 0 ; i < A.size() ; ++i)
        res[i]  = (A[i] + B[i]) ;

    return res;
}

/***************************************/
MatD matrix_sum (const MatD & A, const MatD & B)
{
    if (A. size () != B.size () or A[0]. size () != B[0].size ())
    {
        Rcpp::Rcout << "\nError when summing two matrix, they don't have the same size. \n";
        Rcpp::stop ("\n.");
    }

    MatD RES (A);

    for (unsigned long i = 0 ; i < B.size() ; ++i)
        for (unsigned long j = 0 ; j < B[0].size() ; ++j)
            RES[i][j] += B[i][j];

    return RES;
}


/************************************/
double min_vect (const vector<double> & Vect)
{
    if (0 == Vect.size ())
    {
       Rcpp::Rcout << ("Vector of size null");
       Rcpp::stop ("\n.");
    }
    double min = Vect[0];
     for (auto Val: Vect)
         if ( Val < min)
             min = Val;
     return min;
}

/***********************************************************************************/
double max_vect (const vector<double> & Vect)
{
    if (0 == Vect.size ())
    {
       Rcpp::Rcout << ("Vector of size null");
       Rcpp::stop ("\n.");
    }
    double max = Vect[0];
     for (auto Val: Vect)
         if ( Val > max)
             max = Val;
     return max;
}

double mean_vect (const vector<double> & Vect)
{
    if (0 == Vect.size ())
    {
       Rcpp::Rcout << ("Vector of size null");
      Rcpp::stop ("\n.");
    }
    double mean = 0;
     for (auto val: Vect)
          mean += val;

     return mean / Vect. size ();
}


/************************************/
vector<vector<double>> Normalise (vector<vector<double>> & mat)
{
    if (mat. size () == 0)
    {
        Rcpp::Rcout << "Matrix of size null. \n";
        Rcpp::stop ("\n.");
    }

    mat = Transpose (mat);
    double min, max;
    vector<vector<double>> minMax (2);

    for (auto Vect = mat. begin () ; Vect != mat. end () ; ++Vect)
    {
        min = min_vect (*Vect);
        max = max_vect (*Vect);
        minMax[0].push_back  (min);
        minMax[1].push_back  (max);

        for (auto it = Vect->begin () ; it != Vect->end () ; ++it)
            (*it) = ((*it) - min) / (max - min);
    }
    mat = Transpose (mat);

    return (minMax);
}

/************* average loss function *******************/
vector <double> r_score (const vector<vector<double> > &real, const vector<vector<double> > &pred)
{
    if (pred. size () != real. size ())
    {
        Rcpp::Rcout << "Error in calculating the average_loss function, preds and real have not the same size. \n";
        Rcpp::stop ("\n.");
    }

    vector <double> scores (pred[0]. size (), 0); // final r-square scores
    vector <double> rss (pred[0]. size (), 0); // residual sum of square
    vector <double> tss (pred[0]. size (), 0); // total sum of square
    vector <double> mean_real (pred[0]. size ()); // mean of real values

    //compute the mean of real values
     for (unsigned j = 0; j < real[0]. size (); ++j)
         mean_real [j] = mean_vect (get_col(real, j));

    // compute rss and tss
     for (unsigned i = 0; i < pred. size (); ++i)
         for (unsigned j = 0; j < pred[0]. size (); ++j)
         {
            rss[j] +=  (pred[i][j] - real [i][j]) * (pred[i][j] - real [i][j]);
            tss[j] +=  (real [i][j] - mean_real[j]) * (real [i][j] - mean_real[j]);
        }

    // compute the socres (R-square scorer)
     for (unsigned j = 0; j < real[0]. size (); ++j)
         scores [j] = 1 - (rss[j] / tss[j]);

     return scores;
}


/**************** SIGMOID ACTIVATIOPN FUNCTION *******************/
double sigmoid (double x)
{
    return 1 / (1 + exp (-x));
}

/**************** DERIVATIVE OF SIGMOID FUNCTION *******************/
double deriveSigmoid (double x)
{
    return sigmoid (x) * (1 - sigmoid (x));
}

/**************** RELU ACTIVATIOPN FUNCTIOn *******************/
double relu (double x)
{
    if (x <= 0)
        return 0;
    else
        return x;
}
/**************** DERIVATIVE OF SIGMOID FUNCTION *******************/
double deriveRELU (double x)
{
    if (x <= 0)
        return 0;
    else
        return 1;
}

/************** tanh activation function **************/
double tanh_ac (double x)
{
    return (exp (x) - exp (-x)) / (exp (x) + exp (-x));
}

/**********************************/
double deriveTanh (double x)
{
    return (1 - (tanh_ac(x)*tanh_ac(x)));
}

/***************************************/
VectD sigmoid_v (const VectD & A)
{
     VectD res (A);
     for (double & a : res)
         a = sigmoid (a);
     return res;
}

/***************************************/
VectD relu_v (const VectD & A)
{
    VectD res (A);
    for (double & a : res)
        a = relu (a);
    return res;
}

/***************************************/
VectD tanh_v (const VectD & A)
{
    VectD res (A);
    for (double & a : res)
        a = tanh_ac (a);
    return res;
}
/***************************************/
VectD sigmoid_diff (const VectD & A)
{
     VectD res (A);
     for (double & a : res)
         a = deriveSigmoid (a);
     return res;
}

/***************************************/
VectD relu_diff (const VectD & A)
{
    VectD res (A);
    for (double & a : res)
        a = deriveRELU (a);
    return res;
}

/***************************************/
VectD tanh_diff (const VectD & A)
{
    VectD res (A);
    for (double & a : res)
        a = deriveTanh (a);
    return res;
}

/***********************************/
VectD vect_activation (const VectD & A, const string & f)
{
    VectD res;

    if (f == "linear")
        res = VectD (A);

    else if (f == "sigmoid")
        res = sigmoid_v (A);

    else if (f == "relu")
        res = relu_v (A);

    else if (f == "tanh")
        res = tanh_v (A);

    return res;
}

/***********************************/
VectD diff_activation (const VectD & A, const string & f)
{
    VectD res;
    if (f == "sigmoid")
        res = sigmoid_diff (A);

    else if (f == "relu")
        res = relu_diff (A);

    else if (f == "tanh")
        res = tanh_diff (A);

    else if (f == "linear")
        res = VectD (A. size (), 1);

    return res;
}
