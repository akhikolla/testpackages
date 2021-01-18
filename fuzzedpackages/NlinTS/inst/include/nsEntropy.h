#ifndef NSENTROPY_H
#define NSENTROPY_H

#include <vector>
#include <string>

#define EXP exp(1)

using namespace std;

double myLOG (double x, std::string log = "log2");
double digamma (double x);

typedef vector<int> VectInt;
typedef vector<vector<int>> MatInt;
typedef vector<double> VectD;
typedef vector<VectD> MatD;

namespace nsEntropy{

		VectInt count (const std::vector<int> & X);
		MatInt count (const MatInt & X);

		VectD minMax (const VectD & vect);
		VectD minMax (const std::vector<int> & vect);
		MatD minMax (const MatD & mat);
		void normalize (MatD & mat);

		double entropy (const std::vector<int> & Vect, std::string log = "log");

		double joinProba (std::vector<int> X, std::vector<int> Y, int x, int y);
		double joinProba (MatInt Y, VectInt y);

		double Proba (std::vector<int> X,  int x);

		// make combinations from left to n all sizes from 1 to  n - left + 1
		void generateKCombinations(vector<vector<unsigned>> & combins, vector<unsigned>& tmp, unsigned n, unsigned left, unsigned k);

		// make combinations from left to n all sizes from 1 to  n - left
		vector<vector<unsigned>> generateAllCombinations(unsigned n, unsigned left);


		double joinEntropy (const std::vector<int> & X, const std::vector<int> & Y, std::string log = "log");
		double joinEntropy (const MatInt & Mat, std::string base = "log");


		double condEntropy (const VectInt & X, const VectInt & Y, std::string log = "log");
		double condEntropy (const VectInt & X, const MatInt & Y, std::string log = "log");

		double mutualInformation (const MatInt & X, std::string log = "log",  bool normalize = false);
		double mutualInformation (const std::vector<int> & X, const std::vector<int> & Y, std::string log = "log", bool normalize = false);

		double transferEntropy (const VectInt & X, const VectInt & Y, int p, int q, std::string log = "log", bool normalize = false);


		/*--------------- contnious variables ----------------*/
		double dist (double x, double y);
		double dist (VectD X, VectD Y);

		double entropy (const VectD & V, int k, std::string log = "loge");
		double joinEntropy (const MatD & M, int k, std::string log = "loge");

		double mutualInformation (const MatD & M, int k, std::string alg, bool normalize);
		double transferEntropy (const VectD & X, const VectD & Y, int p=1, int q=1, int k=3, bool normalize = true);
		double transferEntropy_ksg (const VectD & X, const VectD & Y, int p=1, int q=1, int k=3);


		/*template <class type>
		void show (const std::vector<type> & Vect);*/
		/*template <class type>
		void show (const std::vector<std::vector<type>> & matrix)
		{
			unsigned i = 1;
			for (auto & row : matrix)
			{
					cout << i++ << ". ";

				 cout << "[ ";
				 for (auto & val : row)
						 cout << val << "  ";
					cout << "]\n";
			}
		}*/

		/*********************************************************/
		// compute the number of neighbors of a given axis or set of axis
		template <class type>
		VectInt computeNbOfNeighbors (const vector<type> & X, VectD radius, bool equal = false)
		{
			unsigned N = X. size ();
			VectInt NAx (N, 0);

			// we count the number of points Xj whose distance from Xi is strictly less than radius[i],
			for (unsigned i = 0; i < N; i ++)
				for (unsigned j = 0; j < N; j ++)
				{
					if (j != i)
					{
						if (dist (X[i], X[j]) < radius [i] and equal == 0)
							NAx[i] += 1;
						else if (dist (X[i], X[j]) <= radius [i] and equal)
							NAx[i] += 1;
					}
				}
			return NAx;
		}

		/*********************************************************/
		VectInt nbOfNeighborsInRectangle (const MatD & X, const MatD X1, const MatD X2,
								           const VectD & distances);


		/*********************************************************/
		template <class type>
		MatD distanceMatrix (const vector <type> & V)
		{
			unsigned n = V.size ();
			MatD M (n);

			for (unsigned i = 0; i < n; i++)
				M[i]. resize (n, 0);

			for (unsigned i = 0; i < n - 1; i++)
				for (unsigned j = i + 1; j < n; j++){
					M[i][j] = dist (V[i], V[j]);
					M[j][i] = M[i][j];
				}

			return M;
		}

		 /*********************************************************/
		template <class type>
		vector<double> kNearest (const vector<type> & V, int k)
		{
			MatD distMat = distanceMatrix (V);
			vector<double> result (V.size ());

			for (unsigned i = 0; i < V. size (); i++)
			{
				std::sort (distMat[i].begin(), distMat[i].end());
				result [i] = distMat[i][k];
			}

			return result;
		}

		/***************************************************************/
		template <class type>
		vector<vector<type>> lagg (const vector<type> & V, unsigned  p, bool c = 0)
		{
		    unsigned int  N = V.size ();
		    vector<type> P ;

		    // current = 1 if we want to add the  current values of the variable with the lagged ones
		    vector<vector<type>> M (N - p);
		    for (unsigned j = 0 ; j < N - p ; j++)
		    	M[j]. resize (p + c, 0);

		    for (unsigned i = 0 ; i < N - p ; i++)
		    	for (unsigned j = 0; j < p + c ; ++j)
		    		M[i][j] = V[i + (p - 1 + c) - j];

		    return M;
		}

		/*********************************************************/
		template<typename type>
		vector<vector<type>> getCols (const vector<vector<type>> & M, const vector<unsigned> & cols)
		{
			vector<vector<type>> SmallM (M. size ()) ;

			for (unsigned i = 0; i < M. size (); ++i)
				for (auto & idx : cols)
						SmallM[i]. push_back (M[i][idx]);

			return SmallM;
		}

		/*********************************************************/
		template <typename type>
		vector<type> getColumn (const vector<vector<type>> & M, unsigned col)
		{
			vector<type> Vec (M. size ()) ;

			for (unsigned i = 0; i < M. size (); ++i)
						Vec[i] =  M[i][col];

			return Vec;
		}

}// end namespace


#endif
