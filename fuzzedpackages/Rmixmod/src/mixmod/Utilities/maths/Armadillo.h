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

#ifndef mixmod_Armadillo_H
#define	mixmod_Armadillo_H

#include <armadillo>
#include <cmath>

namespace XEM {
namespace MATH {

class DiagonalMatrix {

public:

	DiagonalMatrix(int dim)
	{
		_value = new arma::Col<double>(dim);
	}

	~DiagonalMatrix()
	{
		delete _value;
	}

	double* Store()
	{
		return _value->memptr();
	}

	double& operator[] (int i)
	{
		return (*_value)[i];
	}

	int Nrow()
	{
		return _value->n_rows;
	}
  
  arma::Col<double>* getValue(){
    return _value;
  }

// Fonction pour les tests unitaires, compare les termes de 2 matrices
bool isAlmostEqual(DiagonalMatrix & dm,double epsilon) const {
	bool res=true;
	for (int i=0;i<dm.Nrow();i++){
		if ( fabs( _value->memptr()[i]-dm.getValue()->memptr()[i] )>epsilon) res = false;
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
		_value->memptr()[i]=tmp;i++;
	}
}

	arma::Col<double>* _value;
};

class Matrix {

public:

	Matrix(int nrow, int ncol)
	{
		_value = new arma::Mat<double>(nrow, ncol);
	}

	~Matrix()
	{
		delete _value;
	}

	double* Store()
	{
		return _value->memptr();
	}

	double& operator() (int i, int j)
	{
		return (*_value)(i,j);
	}

	double* GetRow(int index)
	{
		auto row = _value->row(index);
		// TODO: check correctness
		return row.colptr(0);
	}

	int Nrow()
	{
		return _value->n_rows;
	}

	int Ncol()
	{
		return _value->n_cols;
	}

  arma::Mat<double>* getValue(){
    return _value;
  }

// Fonction pour les tests unitaires, compare les valeurs de 2 matrices (on ne tient pas compte du signe
// car on n'utilise cette fonction que pour le calcul du SVD et les vecteurs propres peuvent avoir un signe différent
// suivant la bibliothèque.
bool isAlmostEqual(Matrix & m,double epsilon) const {
	bool res=true;
	int ncol = _value->n_cols;
	for(int j=0; j<ncol ;j++){
		for(int i=0; i<ncol ; i++){
			if ( fabs( fabs((_value)->memptr()[i*(_value->n_cols)+j])-fabs(m.getValue()->memptr()[i*(_value->n_cols)+j]) )>epsilon) res = false;
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
  _value = new arma::Mat<double>(n, n);

	std::ifstream fichier(file,std::ios::in); 
	while (std::getline(fichier, ligne)){
		int j=0;
		// Construction pour une ligne
		std::stringstream ss(ligne);
		// Extraction des données pour une ligne
		while(!ss.eof()){
			double tmp;
			ss >> tmp;
			(_value)->memptr()[i*(_value->n_cols)+j]=tmp;
			j++;
			if (j==n) break;
		}
	i++;
	}
}

	arma::Mat<double>* _value;
};

class SymmetricMatrix {

public:

// nrow == ncol
SymmetricMatrix(int nrow)
{
  _value = new arma::Mat<double>(nrow, nrow);
}

SymmetricMatrix(int nrow, double* store){
  _value = new arma::Mat<double>(nrow, nrow);
  updateData(store);
}


~SymmetricMatrix()
{
  if (_value)
    delete _value;
  if (_store)
    delete[] _store;

}

arma::Mat<double>* getValue() {
  return _value;
} 

int Nrow(){
	return _value->n_rows;
}

double* Store()
{
  double* data = _value->memptr();
  int nrow = _value->n_rows;
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
  int ncol = _value->n_rows;
  for(j=0; j<ncol ;j++){
    for(i=0; i<j+1 ; i++){
      _value->memptr()[i*ncol + j]=store[z];
      _value->memptr()[j*ncol + i]=store[z];
      z++;
    }
  }
}

// get determinant
double determinant(double* store)
{
  updateData(store);
	arma::Mat<double> L;
	arma::Mat<double> U;
	arma::Mat<double> P;
	arma::lu(L, U, P, *_value);
	double logDet = 0.0;
	for (int i=0; i<_value->n_rows; i++)
		logDet += log(fabs(U(i,i)));
  logDet = exp(logDet);
	return logDet;
}

// get inverse
SymmetricMatrix* Inverse(double* store)
{
  updateData(store);
  SymmetricMatrix* inverse = new SymmetricMatrix(_value->n_rows);
	arma::inv(*inverse->_value, *_value);
	return inverse;
}

// compute SVD (only matrices U and D, not V)
void computeSVD(DiagonalMatrix* D, Matrix* U, double* store)
{
  updateData(store);
	arma::Mat<double> V;
	arma::Mat<double> armaU;
	//TODO: tune 500
	if (_value->n_rows < 500) {
		// "small" matrix, standard method
		svd(armaU, *(D->_value), V,*_value, "standard");
	}
	else {
		// (very) "large" matrix, divide & conquer method
		svd(armaU, *(D->_value), V, *_value, "dc");
	}
  int64_t compt=0;
  for (int64_t i=0; i<_value->n_rows; i++) {
    for (int64_t j=0; j<_value->n_rows; j++) {
      U->Store()[compt] = armaU.memptr()[i+j*_value->n_rows];
      compt++;
    }
  }
}

// Fonction pour les tests unitaires, compare les termes de 2 matrices
bool isAlmostEqual(SymmetricMatrix & sm,double epsilon) const {
	bool res=true;
	int ncol = _value->n_cols;
	for(int j=0; j<ncol ;j++){
		for(int i=0; i<j+1 ; i++){
			if ( fabs( (_value)->memptr()[i*(_value->n_cols)+j]-sm.getValue()->memptr()[i*(_value->n_cols)+j] )>epsilon) res = false;
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
		
  _value = new arma::Mat<double>(n, n);

	std::ifstream fichier(file,std::ios::in); 
	while (std::getline(fichier, ligne)){
		int j=0;
		// Construction pour une ligne
		std::stringstream ss(ligne);
		// Extraction des données pour une ligne
		while(!ss.eof()){
			double tmp;
			ss >> tmp;
			(_value)->memptr()[i*(_value->n_cols)+j]=tmp;
			j++;
			if (j==n) break;
		}
	i++;
	}
}


private:

	arma::Mat<double>* _value;
  double* _store;
};

}
}

#endif
