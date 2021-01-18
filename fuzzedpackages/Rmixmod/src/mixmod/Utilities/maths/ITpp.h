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

#ifndef mixmod_ITpp_H
#define	mixmod_ITpp_H

#include <itpp/itbase.h>
#include <cmath>

namespace XEM {
namespace MATH {

class DiagonalMatrix {

public:

	DiagonalMatrix(int dim)
	{
		_value = new itpp::Vec<double>(dim);
	}

	~DiagonalMatrix()
	{
		delete _value;
	}

	double* Store()
	{
		return _value->_data();
	}

	double& operator[] (int i)
	{
		return (*_value)[i];
	}

	int Nrow()
	{
		return _value->size();
	}

  itpp::Vec<double>* getValue(){
    return _value;
  }

// Fonction pour les tests unitaires, compare les termes de 2 matrices
bool isAlmostEqual(DiagonalMatrix & dm,double epsilon) const {
	bool res=true;
	for (int i=0;i<dm.Nrow();i++){
		if ( fabs( _value->_data()[i]-dm.getValue()->_data()[i] )>epsilon) res = false;
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
		_value->_data()[i]=tmp;i++;
	}
}
  
	itpp::Vec<double>* _value;
};

class Matrix {
	
public:
	
	Matrix(int nrow, int ncol)
	{
		_value = new itpp::Mat<double>(nrow, ncol);
	}
	
	~Matrix()
	{
		delete _value;
	}

	double* Store()
	{
		return _value->_data();
	}
	
	double& operator() (int i, int j)
	{
		return (*_value)(i,j);
	}

	double* GetRow(int index)
	{
		return _value->get_row(index)._data();
	}
	
	int Nrow()
	{
		return _value->rows();
	}
	
	int Ncol()
	{
		return _value->cols();
	}

  itpp::Mat<double>* getValue(){
    return _value;
  }

// Fonction pour les tests unitaires, compare les valeurs de 2 matrices (on ne tient pas compte du signe
// car on n'utilise cette fonction que pour le calcul du SVD et les vecteurs propres peuvent avoir un signe différent
// suivant la bibliothèque.
bool isAlmostEqual(Matrix & m,double epsilon) const {
	bool res=true;
	int ncol = _value->cols();
	for(int j=0; j<ncol ;j++){
		for(int i=0; i<ncol ; i++){
			if ( fabs( fabs((_value)->_data()[i*(_value->cols())+j])-fabs(m.getValue()->_data()[i*(_value->cols())+j]) )>epsilon) res = false;
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

  _value = new itpp::Mat<double>(n, n);

	std::ifstream fichier(file,std::ios::in); 
	while (std::getline(fichier, ligne)){
		int j=0;
		// Construction pour une ligne
		std::stringstream ss(ligne);
		// Extraction des données pour une ligne
		while(!ss.eof()){
			double tmp;
			ss >> tmp;
			(_value)->_data()[i*(_value->cols())+j]=tmp;
			j++;
			if (j==n) break;
		}
	i++;
	}
}
	
	itpp::Mat<double>* _value;
};

class SymmetricMatrix {
	
public:
	
// nrow == ncol
SymmetricMatrix(int nrow)
{
	_value = new itpp::Mat<double>(nrow, nrow);
}

SymmetricMatrix(int nrow, double* store)
{
	_value = new itpp::Mat<double>(nrow, nrow);
  updateData(store);
}

~SymmetricMatrix()
{
	delete _value;
  if (_store)
    delete[] _store;
}

itpp::Mat<double>* getValue() {
  return _value;
} 

int Nrow(){
	return _value->rows();
}

double* Store()
{
  double* data = _value->_data();
  int nrow = _value->rows();
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
  int ncol = _value->rows();
  for(j=0; j<ncol ;j++){
    for(i=0; i<j+1 ; i++){
      _value->_data()[i*ncol + j]=store[z];
      _value->_data()[j*ncol + i]=store[z];
      z++;
    }
  }
}

// get determinant 
double determinant(double* store)
{
  updateData(store);
	itpp::Mat<double> L(_value->rows(), _value->cols());
	itpp::Mat<double> U(_value->rows(), _value->cols());
	itpp::ivec p(_value->rows());
	itpp::lu(*_value, L, U, p);
	double logDet = 0.0;
	for (int i=0; i<_value->rows(); i++)
		logDet += log(fabs(U(i,i)));
  logDet = exp(logDet);
	return logDet;
}
	
// get inverse
SymmetricMatrix* Inverse(double* store)
{
  updateData(store);
  SymmetricMatrix* inverse = new SymmetricMatrix(_value->rows());
	itpp::inv(*_value, *inverse->_value);
	return inverse;
}
	
// compute SVD (only matrices U and D, not V)
void computeSVD(DiagonalMatrix* D, Matrix* U, double* store)
{
  updateData(store);
	itpp::Mat<double> V_itpp(U->Ncol(), U->Ncol());
	itpp::Mat<double> U_itpp(U->Ncol(), U->Ncol());
	itpp::svd(*_value, U_itpp, *(D->_value), V_itpp);
  int64_t compt=0;
  for (int64_t i=0; i<_value->rows(); i++) {
    for (int64_t j=0; j<_value->rows(); j++) {
      U->Store()[compt] = U_itpp._data()[i+j*_value->rows()];
      compt++;
    }
  }

}

// Fonction pour les tests unitaires, compare les termes de 2 matrices
bool isAlmostEqual(SymmetricMatrix & sm,double epsilon) const {
	bool res=true;
	int ncol = _value->cols();
	for(int j=0; j<ncol ;j++){
		for(int i=0; i<j+1 ; i++){
			if ( fabs( (_value)->_data()[i*(_value->cols())+j]-sm.getValue()->_data()[i*(_value->cols())+j] )>epsilon) res = false;
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
		
	_value = new itpp::Mat<double>(n, n);

	std::ifstream fichier(file,std::ios::in); 
	while (std::getline(fichier, ligne)){
		int j=0;
		// Construction pour une ligne
		std::stringstream ss(ligne);
		// Extraction des données pour une ligne
		while(!ss.eof()){
			double tmp;
			ss >> tmp;
			(_value)->_data()[i*(_value->cols())+j]=tmp;
			j++;
			if (j==n) break;
		}
	i++;
	}
}
	
private:
	
	itpp::Mat<double>* _value;
  double* _store;
};

}
}

#endif
