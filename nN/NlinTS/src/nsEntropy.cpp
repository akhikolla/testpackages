//#include <Rcpp.h>
#include <cstdlib>
#include <math.h>
#include <algorithm>
#include "../inst/include/nsEntropy.h"


using namespace std;

/*****************************************************/
// Extract distinct values

//using namespace nsEntropy;
//namespace nsEntropy {

double myLOG (double x, std::string log)
{
  if (x == 0)
	return 0;

  if (log == "loge")
    return (log2 (x) / log2(EXP));
  else if (log == "log10")
    return (log2 (x) / log2(10));
  else if (log == "log2")
    return log2 (x);
  else
    return log2 (x);
}

/*****************************************************/
double digamma (double x)
{
  double a = 0, b, c;

  while (x <= 5)
  {
    a -= 1 / x;
    x += 1;
  }

  b = 1 / (x * x);

  c = b * (-1/12.0 +
      b * (1/120.0 +
      b * (-1/252.0 +
      b * (1/240.0 +
      b * (-1/132.0 +
      b * (691/32760.0 +
      b * (-1/12.0 +
      b * 3617/8160.0)))))));

  return (a + log (x) - 0.5 / x + c);

}

/*****************************************************/
	// Min and max values of columns of a matrix
VectD nsEntropy::minMax (const VectD & vect)
{
	VectD result (2);
	result [0] = vect [0];
	result [1] = vect [0];

	for (unsigned i = 1; i < vect. size (); ++i)
	{
		if (vect [i] < result [0])
			result [0] = vect [i];

		if (vect [i] > result [1])
			result [1] = vect [i];
	}
	return result;
}
/*****************************************************/
	// Min and max values of columns of a matrix
VectD nsEntropy::minMax (const vector<int> & vect)
{
	VectD result (2);
	result [0] = vect [0];
	result [1] = vect [0];

	for (unsigned i = 1; i < vect. size (); ++i)
	{
		if (vect [i] < result [0])
			result [0] = vect [i];

		if (vect [i] > result [1])
			result [1] = vect [i];
	}
	return result;
}
/*****************************************************/
MatD nsEntropy::minMax (const MatD & mat)
{
	MatD result (mat[0]. size ());
	for (unsigned j = 0; j < mat[0]. size (); ++j)
	{
		result [j] = minMax (getColumn (mat, j));
	}
	return result;
}

/*****************************************************/
void nsEntropy::normalize (MatD & mat)
{
	MatD min_max = minMax (mat);
	unsigned m = mat[0]. size (), n = mat. size ();

		for (unsigned j = 0; j < m; ++j)
		{
			if (min_max[j][0] != min_max [j][1])
			{
				for (unsigned i = 0; i < n; ++i)
					mat[i][j] = (mat[i][j] - min_max[j][0]) / (min_max[j][1] - min_max[j][0]);
			}

		}
}

/*****************************************************/
//make combinations from left to n of size k
void nsEntropy::generateKCombinations(vector<vector<unsigned>> & combins, vector<unsigned>& tmp, unsigned n, unsigned left, unsigned k)
{
    // record the combination
    if (k == 0) {
        combins.push_back(tmp);
        return;
    }

    for (unsigned i = left; i <= n; ++i)
    {
        tmp.push_back(i);
        generateKCombinations(combins, tmp, n, i + 1, k - 1);
        tmp.pop_back();
    }
}

/*****************************************************/
// make combinations from left to n all sizes from 1 to  n - left + 1
vector<vector<unsigned> > nsEntropy::generateAllCombinations(unsigned n, unsigned left)
{
    vector<vector<unsigned>> combins;
    for (unsigned k (1); k <= n - left + 1; ++k)
    {
         vector<unsigned> tmp;
         generateKCombinations (combins, tmp, n, left, k);
    }

    return combins;
}

/*****************************************************/
vector<int> nsEntropy::count (const std::vector<int> & X)
{
	std::vector<int> Vect (X);
	std::vector<int>::iterator it;
	std::sort (Vect.begin(), Vect.end());
	it = std::unique (Vect. begin (), Vect. end ());

		// Distinct values
	Vect.resize (std::distance (Vect.begin(), it));
	return Vect;
}

/*****************************************************/
MatInt nsEntropy::count (const MatInt & X)
{
	vector<vector<int>> Vect (X);
	vector<vector<int>>::iterator it;
	std::sort (Vect.begin(), Vect.end());
	it = std::unique (Vect. begin (), Vect. end ());

		// Distinct values
	Vect.resize (std::distance (Vect.begin(), it));
	return Vect;
}

/*****************************************************/
double nsEntropy::joinProba (vector<int> X, vector<int> Y, int x, int y)
{
	double J = 0;
	for (unsigned i = 0; i < X.size (); ++i)
	{
		if (X[i] == x and Y[i] == y)
			J++;
	}
	return (J / X.size ());
}

/*****************************************************/
double nsEntropy::joinProba (MatInt Y, VectInt y)
{
	double J = 0;
	unsigned j;

	for (unsigned i = 0; i < Y.size (); ++i)
	{
		for (j = 0; j < Y [0] .size (); ++j)
			if (Y[i][j] != y[j])
        break;

		if (j == Y [0].size ())
			J = J + 1.0;
	}
	return (J / Y.size ());
}

/*****************************************************/
double nsEntropy::Proba (vector<int> X, int x)
{
	double P = 0;
	for (unsigned i = 0; i < X.size (); ++i)
	{
		if (X[i] == x)
			P++;
	}
	return (P / X.size ());
}


/*****************************************************/
double nsEntropy::entropy (const VectInt & X, string log)
{
	double E = 0, x;

	vector<int> Vect = count (X);
	unsigned n = Vect.size ();
	for (unsigned i = 0; i < n; i++)
	{
			x = Proba (X, Vect[i]);
			if (x > 0)
				E += x  * myLOG (x, log) ;
	}

	return -E;
}

/****************************/
double nsEntropy::joinEntropy (const VectInt & X1, const VectInt & X2, string log)
{
	double J = 0;
	double x;

	vector<int> X = count (X1);
	vector<int> Y = count (X2);
	unsigned n = X. size ();
	unsigned m = Y. size ();

	for (unsigned i = 0; i < n; i++)
	{
		for (unsigned j = 0; j < m; j++) {
			x = joinProba (X1, X2, X[i], Y[j]);

			if (x > 0)
				J = J + x  * myLOG (x, log);
		}
	}
	return -J;
}

/*****************************************************/
double nsEntropy::joinEntropy (const MatInt & Mat, string log)
{
	double J = 0;
	double x;

	MatInt tuples = count (Mat);

	for (auto tuple : tuples)
	{
		x = joinProba (Mat, tuple);
		if (x > 0)
			J += x  * myLOG (x, log);
	}
	return -J;
}

/*****************************************************/
double nsEntropy::condEntropy (const VectInt & X, const VectInt & Y, string log)
{
	return joinEntropy (X, Y, log) - entropy (Y, log);
}

/*****************************************************/
double nsEntropy::condEntropy (const VectInt & X, const MatInt & Y, string log)
{
	MatInt M = Y;
	M .push_back (X);

	return (joinEntropy (M, log) - joinEntropy (Y, log));
}

/*****************************************************/
double nsEntropy::mutualInformation (const VectInt & X, const VectInt & Y, std::string log, bool normalize)
{
  double je = joinEntropy (X, Y, log);
  double eX = entropy (X, log), eY = entropy (Y, log);
  double mi = eX + eY  - je;

  if (normalize && max (eX, eY) > 0)
    mi /= max (eX, eY);

	return  mi;
}

/*****************************************************/
double nsEntropy::mutualInformation (const MatInt & X, std::string log, bool normalize)
{
  double mi = 0, max = 0;
  unsigned n_cols = unsigned (X [0]. size ());

  // compute the join entropy of all combinations of columns of X
  vector<vector<unsigned>> combins = nsEntropy::generateAllCombinations (unsigned (n_cols - 1), 0);

  for (auto & combin : combins)
  {
      if (combin. size () == 1)
      {
        double e =  entropy (getColumn (X, combin[0]), log);
        if (e > max)
          max = e;
        mi -= e;
      }
      else
      {
        mi += joinEntropy (getCols (X, combin), log) *  pow (-1, combin. size ());
      }
  }

	if (normalize && (max > 0))
		mi /= max;

	return  -mi;
}
/*****************************************************/
double nsEntropy::transferEntropy (const VectInt & X, const VectInt & Y, int p, int q, std::string log, bool normalize)
{
	double te, denom;
	MatInt Xp, Xm, Ym, XmYm, XpYm;
	Xm = lagg (X, p, 0);
	Xp = lagg (X, p, 1);
	Ym = lagg (Y, q, 0);

	XmYm = Xm;
	XpYm = Xp;

	// Resize the join matrix to have the same lenght if p # q
	if ((p - q) < 0)
	{
		XmYm. erase (XmYm.begin (), XmYm.begin()  + q - p);
		XpYm. erase (XpYm.begin (), XpYm.begin()  + q - p);
	}
	if ((p - q) > 0)
		Ym. erase (Ym.begin (), Ym.begin()  + p - q);

	unsigned N = Ym. size ();

	for (unsigned i = 0 ; i < N; ++i)
		for (unsigned j = 0 ; j < Ym [0]. size (); ++j){
			XmYm [i]. push_back ( Ym [i][j] );
			XpYm [i]. push_back ( Ym [i][j] );
		}

	denom = joinEntropy (Xp, log) - joinEntropy (Xm, log);
	te = denom - joinEntropy (XpYm, log) + joinEntropy (XmYm, log);

	// normalisation: deviding TE by max entropy of X: log (max(X) - min (X))
	if (normalize and denom != 0)
	{
      VectD min_max = minMax (X);
      double H0 = myLOG (abs (min_max[1] - min_max[0]) , log);
      
      if (H0 != 0)
      	te = te / H0;
  	}
	//te = te / denom; // this is another method but not precise

	return te;
}

/*********************************************************/
/*--------------- continuous variables ------------------*/
/*********************************************************/

/*****************************************************/
double nsEntropy::dist (double x, double y) {
	return abs (x - y);
}

double nsEntropy::dist (VectD X, VectD Y)
{
	double distance = 0;

	for (unsigned i = 0; i < X. size (); ++i){
		  if (distance < abs (X[i] - Y[i]))
			distance = abs (X[i] - Y[i]);
		}
	return  distance;
}


 /*********************************************************/

double nsEntropy::entropy (const VectD & V, int k, std::string log)
{
	double E = 0;
	unsigned N = V. size ();
	double sum = 0;
	double cd = 1;

	VectD distances = kNearest (V, k);


	for (unsigned i = 0; i < N; i ++){
		sum += myLOG (2 * distances [i], log);
		//cout << distances[i] << "   " << log2 (2 * distances [i]) << endl;
	}

	sum = sum / N;

	E = digamma (N) - digamma (k) +  sum  + myLOG (cd, log);
	return E ;
}

 /*********************************************************/
double nsEntropy::joinEntropy (const MatD & M, int k, std::string log)
{
	double E = 0;
	unsigned N = M. size ();
	double sum = 0;

	unsigned d = M[0]. size ();

	VectD distances = kNearest (M, k);

	for (unsigned i = 0; i < N; i ++)
		sum += myLOG (2 * distances [i], log);

	sum = sum * d / N;

	E = digamma (N) - digamma (k) +  sum ;
	return E ;
}

/*******************************************************/
VectInt nsEntropy::nbOfNeighborsInRectangle (const MatD & X, const MatD X1, const MatD X2,
						           const VectD & distances)
{
	unsigned N = X. size ();
	VectInt Nx (N, 0);

	double dLocalx, dLocaly, dx, dy;

	for (unsigned i = 0 ; i < N; ++i)
	{
		dx = 0;
		dy = 0;
		// Compute dx and dy : edge lengths of the hyper-rectangle
		for (unsigned j = 0 ; j < N; ++j)
		{
			dLocalx = dist (X1[i], X1[j]);
			dLocaly = dist (X2[i], X2[j]);

			if (dist (X[i], X[j]) <= distances [i] and j != i)
			{
				if (dx < dLocalx)
					dx = dLocalx;
				if (dy < dLocaly)
					dy = dLocaly;
			}
		}

		// Count the number of points in the hyper-rectangle
		for (unsigned j = 0 ; j < N; ++j)
		{
			if (dist (X1[i], X1[j]) <= (2 * dx) and  dist (X2[i], X2[j]) <= (2 * dy) )
					Nx[i] += 1;
		}
	}
	return Nx;
}

/*****************************************************/
// M: array of dimension (x, 2), computing mi between the two columns

double nsEntropy::mutualInformation (const MatD & M, int k, string alg, bool normalize)
{
	double mi = 0;
	unsigned N = M. size ();
	double sum = 0;

	//std::cout << "Size: " << N << " " << M[0]. size () << "\n";

	VectInt NX, NY;
	VectD X, Y;
	X = getColumn (M, 0);
	Y = getColumn (M, 1);

	VectD distances = kNearest (M, k);

	if (alg == "ksg1")
	{
		NX = computeNbOfNeighbors (X, distances, false);
		NY = computeNbOfNeighbors (Y, distances, false);

		for (unsigned i = 0; i < N; i ++){

			sum += digamma (NX[i] + 1) + digamma (NY[i] + 1);
		}

		sum = sum / N;
		mi = digamma (k) + digamma (N) - sum;

	}

	else if (alg == "ksg2")
	{
		VectD distances_x = kNearest (X, k);
		VectD distances_y = kNearest (Y, k);
		NX = computeNbOfNeighbors (X, distances_x, true);
		NY = computeNbOfNeighbors (Y, distances_y, true);

		for (unsigned i = 0; i < N; i ++)
			sum += digamma (NX[i]) + digamma (NY[i]);

		sum = sum  / N;
		mi = digamma (k) - (1.0 / k) + digamma (N) - sum;
	}

  // Normalizing mutual information by divide it by the joint entropy
  if (normalize)
  {
    double jointEn = 0;
    for (double d: distances)
      jointEn += d;
    jointEn *= (2.0 / distances. size ());
    jointEn +=  digamma (N) - digamma (k);
    mi = mi / jointEn;
  }
	return mi;
}

/***************************************************************/
// Tranfer entropy from Y to X
double nsEntropy::transferEntropy (const VectD & X, const VectD & Y, int p, int q, int k, bool normalize)
{
	double te, sum = 0;


	MatD Xm, Xp, Ym, XpYm, XmYm;
	VectInt NXm, NXmYm, NXpYm, NXp;

	Xm = lagg (X, p, 0);
	Ym = lagg (Y, q);
	Xp = lagg (X, p, 1);

	// Resize the join matrix to have the same lenght if p # q
	if ((p - q) < 0)
	{
		Xm. erase (Xm.begin (), Xm.begin()  + q - p);
		Xp. erase (Xp.begin (), Xp.begin()  + q - p);
	}
	if ((p - q) > 0)
		Ym. erase (Ym.begin (), Ym.begin()  + p - q);

	//the number of observation of lagged variables : n - max (p, q)
	unsigned N = Ym.size ();
	XpYm = Xp;
	XmYm = Xm;

	for (unsigned i = 0 ; i < N; ++i)
		for (unsigned j = 0 ; j < Ym [0]. size (); ++j){
			XmYm [i]. push_back ( Ym [i][j] );
			XpYm [i]. push_back ( Ym [i][j] );
		}

	//  distances from k neighbors of the join matrix Xp  (Xcurrent + Xpassed + Ypassed)
	VectD distances = kNearest (XmYm, k);

	// we count the number of local points relative to  the marginal matrices
	NXmYm = computeNbOfNeighbors (XmYm, distances);
	NXm = computeNbOfNeighbors (Xm, distances);
	NXp = computeNbOfNeighbors (Xp, distances);

	for (unsigned i = 0; i < N; i ++){
		sum += digamma (NXm[i] + 1) - digamma (NXmYm[i] + 1) - digamma (NXp[i] + 1);
	}

	te = digamma (k) + (sum / N) ;

	
	// Compute  NTE <- TE / (H0 - H(Xp|Xm,Ym))
	if (normalize == true)
	{
		double denom = 0, H0;

		//VectD Xt = getColumn (Xp, 0);
		VectD min_max = minMax (X);

		H0 = myLOG (abs (min_max[1] - min_max[0]) , "loge");


		for (unsigned i = 0; i < N; i ++)
			denom +=  myLOG (2*distances[i], "loge") + digamma (NXmYm[i] + 1);

		denom = H0 - ( (denom / N) - digamma (k) );

		if (denom != 0)
			te = te / denom;
	}

	return te;
}


/**************************************************************************/
double nsEntropy::transferEntropy_ksg (const VectD & X, const VectD & Y, int p, int q, int k)
{
	double te, sum = 0;

	VectInt NXm, NXmYm, NXpYm;
	MatD Xm, Ym, XpYm, XmYm, Xp;

	Xm = lagg (X, p, 0);
	Ym = lagg (Y, q);
	Xp = lagg (X, p, 1);

	// Resize the join matrix to have the same lenght if p # q
	if ((p - q) < 0)
	{
		Xm. erase (Xm.begin (), Xm.begin()  + q - p);
		Xp. erase (Xp.begin (), Xp.begin()  + q - p);
	}
	if ((p - q) > 0)
		Ym. erase (Ym.begin (), Ym.begin()  + p - q);

	//the number of observation of lagged variables : n - max (p, q)
	unsigned N = Ym.size ();
	XpYm = Xp;
	XmYm = Xm;

	for (unsigned i = 0 ; i < N; ++i)
		for (unsigned j = 0 ; j < Ym [0]. size (); ++j){
			XmYm [i]. push_back ( Ym [i][j] );
			XpYm [i]. push_back ( Ym [i][j] );
		}

	//  distances from k neighbors of the join matrix Xp  (Xcurrent + Xpassed + Ypassed)
	VectD distances = kNearest (XpYm, k);

	// We countthe number of points within the hyper-rectangle equal to the Cartesian prod of  XmYm, Xp and Xm
	NXmYm = nbOfNeighborsInRectangle (XpYm, Xm, Ym, distances);
	NXpYm = nbOfNeighborsInRectangle (XpYm, Xp, Xm, distances);
	NXm = nbOfNeighborsInRectangle (XpYm, Xm, Xm, distances);
	//NXm = computeNbOfNeighbors (Xm, distances, true);

	for (unsigned i = 0; i < N; i ++)
		sum += digamma (NXm[i]) - digamma (NXmYm[i]) - digamma (NXpYm[i]) + (1 / NXpYm[i]) + (1 / NXmYm[i]);


	te = digamma (k) - (2 / k) + sum / N; // - (lag * epsilon / N);


	return te;
}

//} // end namespace
