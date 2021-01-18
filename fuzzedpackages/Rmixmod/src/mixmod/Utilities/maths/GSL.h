//****************************************************************************************|
//	file      : mixmod / src/mixmod/Utilities/Maths/Eigen.h
//	copyright : (C) MIXMOD Team - 2001-2013
//	email     : contact@mixmod.org
//========================================================================================
//	This file is part of MIXMOD (see <http://www.mixmod.org>)
//
//	MIXMOD is free software: you can redistribute it and/or modify it under the terms of
//	the GNU General Public License as published by the Free Software Foundation,
//	either version 3 of the License, or (at your option) any later version.
//
//	MIXMOD is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
//	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//	See the GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License along with MIXMOD.
//	If not, see <http://www.gnu.org/licenses/>.
//****************************************************************************************|

#ifndef mixmod_GSL_H
#define	mixmod_GSL_H

#include <gsl/gsl_linalg.h>
#include <cmath>

namespace XEM {
namespace MATH {

class DiagonalMatrix {

public:

	DiagonalMatrix(int dim)
	{
		_value = gsl_vector_alloc(dim);
	}

	~DiagonalMatrix()
	{
		gsl_vector_free(_value);
	}

	double* Store()
	{
		return _value->data;
	}

	double& operator[] (int i)
	{
		return _value->data[i];
		//return gsl_vector_get(_value, i); //Does not return a reference
	}

	int Nrow()
	{
		return _value->size;
	}

gsl_vector* getValue(){
	return _value;
}

// Fonction pour les tests unitaires, compare les termes de 2 matrices
bool isAlmostEqual(DiagonalMatrix & dm,double epsilon) const {
	bool res=true;
	for (int i=0;i<dm.Nrow();i++){
		if ( fabs( _value->data[i]-dm.getValue()->data[i] )>epsilon) res = false;
	}
	return res;
}

// Fonction pour les tests unitaires, remplissage d'une matrice à partir d'un fichier .txt
void fillMatrix(std::string file){
	std::ifstream fichiertmp(file,std::ios::in); // on ouvre en lecture
	int i=0;
 
	// On verifie que le fichier est bien ouvert ...
	if(!fichiertmp)
	{
		std::cerr << "[ERROR] Impossible d'ouvrir le fichier " << file << std::endl;
	}
	
	std::string ligne;
	int n=0;
	while(std::getline(fichiertmp, ligne)) n++;  //On lit chaque ligne du fichier que l'on stoke dans "ligne"
		

	std::ifstream fichier(file,std::ios::in); 
	while (std::getline(fichier, ligne)){
		std::stringstream ss(ligne);
		double tmp;
		ss >> tmp;
		_value->data[i]=tmp;i++;
	}
}

	gsl_vector* _value;
};

class Matrix {
	
public:
	
	Matrix(int nrow, int ncol)
	{
		_value = gsl_matrix_alloc(nrow, ncol);
	}
	
	~Matrix()
	{
		gsl_matrix_free(_value);
	}

	double* Store()
	{
		return _value->data;
	}
	
	double& operator() (int i, int j)
	{
		return _value->data[i*(_value->size2)+j];
		//return gsl_matrix_get(_value, i, j); //Does not return a reference
	}

	double* GetRow(int index)
	{
		return gsl_matrix_row(_value, index).vector.data;
	}
	
	int Nrow()
	{
		return _value->size1;
	}
	
	int Ncol()
	{
		return _value->size2;
	}
  
  gsl_matrix* getValue() {
    return _value;
  }

// Fonction pour les tests unitaires, compare les valeurs de 2 matrices (on ne tient pas compte du signe
// car on n'utilise cette fonction que pour le calcul du SVD et les vecteurs propres peuvent avoir un signe différent
// suivant la bibliothèque.
bool isAlmostEqual(Matrix & m,double epsilon) const {
	bool res=true;
	int ncol = _value->size2;
	for(int j=0; j<ncol ;j++){
		for(int i=0; i<ncol ; i++){
			if ( fabs( fabs((_value)->data[i*(_value->size2)+j])-fabs(m.getValue()->data[i*(_value->size2)+j]) )>epsilon) res = false;
		}
	}
	return res;
}
	
// Fonction pour les tests unitaires, remplissage d'une matrice à partir d'un fichier .txt
void fillMatrix(std::string file){
	std::ifstream fichiertmp(file,std::ios::in); // on ouvre en lecture
	int i=0;
 
	// On verifie que le fichier est bien ouvert ...
	if(!fichiertmp)
	{
		std::cerr << "[ERROR] Impossible d'ouvrir le fichier " << file << std::endl;
	}
	std::string ligne;
	int n=0;
	while(std::getline(fichiertmp, ligne)) n++;  //On lit chaque ligne du fichier que l'on stoke dans "ligne"
  _value = gsl_matrix_alloc(n, n);

	std::ifstream fichier(file,std::ios::in); 
	while (std::getline(fichier, ligne)){
		int j=0;
		// Construction pour une ligne
		std::stringstream ss(ligne);
		// Extraction des données pour une ligne
		while(!ss.eof()){
			double tmp;
			ss >> tmp;
			(_value)->data[i*(_value->size2)+j]=tmp;
			j++;
			if (j==n) break;
		}
	i++;
	}
}
	
	gsl_matrix* _value;
};

class SymmetricMatrix {
	
public:
	
// nrow == ncol
SymmetricMatrix(int nrow)
{
	_value = gsl_matrix_alloc(nrow, nrow);
}

SymmetricMatrix(int nrow,double* store){
  _value = gsl_matrix_alloc(nrow, nrow);
  updateData(store);
}
	

~SymmetricMatrix()
{
  gsl_matrix_free(_value);
  if (_store)
    delete[] _store;
}
	
gsl_matrix* getValue() {
  return _value;
} 

int Nrow(){
	return _value->size1;
}

double* Store()
{
  double* data = _value->data;
  int nrow = _value->size1;
  int i,j;
  int z=0;
  _store = new double [nrow*(nrow + 1)/2];
  for(j=0; j<nrow ;j++){
    for(i=0; i<j+1 ; i++){
      _store[z]=data[i*nrow+j];z++;
    }
  }
  return _store;
}
	
void updateData(double* store){
  int i,j;
  int z=0;
  int ncol = _value->size1;
  for(j=0; j<ncol ;j++){
    for(i=0; i<j+1 ; i++){
      _value->data[i*ncol + j]=store[z];
      _value->data[j*ncol + i]=store[z];
      z++;
    }
  }
}

// get determinant
double determinant(double* store)
{
  updateData(store);
  int nrow = _value->size1;
  gsl_permutation* permutation = gsl_permutation_alloc(nrow);
  int signum;
  gsl_matrix* copyValue = gsl_matrix_alloc(nrow, nrow); 
  gsl_matrix_memcpy(copyValue, _value);
  gsl_linalg_LU_decomp(copyValue, permutation, &signum);
  double logDet = exp(gsl_linalg_LU_lndet(copyValue));
  gsl_permutation_free(permutation);
  gsl_matrix_free(copyValue);
  return logDet; 	
}

// get inverse
SymmetricMatrix* Inverse(double* store)
{
  updateData(store);
  int nrow = _value->size1;
  gsl_permutation* permutation = gsl_permutation_alloc(nrow);
  int signum;
  gsl_matrix* LUdecomp = gsl_matrix_alloc_from_matrix(_value, 0, 0, nrow, nrow);
  gsl_linalg_LU_decomp(LUdecomp, permutation, &signum);
  SymmetricMatrix* inverse = new SymmetricMatrix(nrow);
  gsl_linalg_LU_invert(LUdecomp, permutation, inverse->_value);
  gsl_permutation_free(permutation);
  gsl_matrix_free(LUdecomp);
  return inverse;
}

// compute SVD (only matrices U and D, not V) 
void computeSVD(DiagonalMatrix* D, Matrix* U, double* store)
{
  updateData(store);
  int nrow = _value->size1;
  // copy current matrix into U
  for (int i=0; i<nrow; i++) {
    for (int j=0; j<nrow; j++)
      U->_value-> data[i*nrow + j] = _value-> data[i*nrow + j];
  }
  gsl_matrix* V = gsl_matrix_alloc(nrow, nrow);
  gsl_vector* work = gsl_vector_alloc(nrow);
  gsl_linalg_SV_decomp (U->_value, V, D->_value, work);
  gsl_matrix_free(V);
  gsl_vector_free(work);
}

// Fonction pour les tests unitaires, compare les termes de 2 matrices
bool isAlmostEqual(SymmetricMatrix & sm,double epsilon) const {
	bool res=true;
	int ncol = _value->size2;
	for(int j=0; j<ncol ;j++){
		for(int i=0; i<j+1 ; i++){
			if ( fabs( (_value)->data[i*(_value->size2)+j]-sm.getValue()->data[i*(_value->size2)+j] )>epsilon) res = false;
		}
	}
	return res;
}

// Fonction pour les tests unitaires, remplissage d'une matrice à partir d'un fichier .txt
void fillMatrix(std::string file){
	std::ifstream fichiertmp(file,std::ios::in); // on ouvre en lecture
	int i=0;
 
	// On verifie que le fichier est bien ouvert ...
	if(!fichiertmp)
	{
		std::cerr << "[ERROR] Impossible d'ouvrir le fichier " << file << std::endl;
	}
	
	std::string ligne;
	int n=0;
	while(std::getline(fichiertmp, ligne)) n++;  //On lit chaque ligne du fichier que l'on stoke dans "ligne"
		
  _value = gsl_matrix_alloc(n, n);

	std::ifstream fichier(file,std::ios::in); 
	while (std::getline(fichier, ligne)){
		int j=0;
		// Construction pour une ligne
		std::stringstream ss(ligne);
		// Extraction des données pour une ligne
		while(!ss.eof()){
			double tmp;
			ss >> tmp;
			(_value)->data[i*(_value->size2)+j]=tmp;
			j++;
			if (j==n) break;
		}
	i++;
	}
}
	

private:
	
	gsl_matrix* _value;
  double* _store;
};

}
}

#endif
