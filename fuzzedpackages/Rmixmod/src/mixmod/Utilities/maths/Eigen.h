/***************************************************************************
                             SRC/mixmod/Utilities/maths/Eigen.h  description
    copyright            : (C) MIXMOD Team - 2001-2016
    email                : contact@mixmod.org
 ***************************************************************************/

/***************************************************************************
    This file is part of MIXMOD
    
    MIXMOD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MIXMOD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MIXMOD.  If not, see <http://www.gnu.org/licenses/>.

    All informations available on : http://www.mixmod.org
***************************************************************************/
#ifndef XEM_MATH_EIGEN_H
#define XEM_MATH_EIGEN_H

#ifdef RPACKAGE
#include <RcppEigen.h>
#else
#include <Eigen/Dense>
#include <Eigen/Core>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>

namespace XEM {
namespace MATH {

// TODO: copy constructors ?

class DiagonalMatrix {
	
public:
	
DiagonalMatrix(int dim) {
	_dim = dim;
	_value = new double[dim];
}
	
~DiagonalMatrix() {
	if (_value)
		delete [] _value;
}
	
double* Store() {
	return _value;
}
	
int Nrow() {
	return _dim;
}

double* getValue(){
	return _value;
}

// Fonction pour les tests unitaires, compare les termes de 2 matrices
bool isAlmostEqual(DiagonalMatrix & dm,double epsilon) const {
	bool res=true;
	for (int i=0;i<dm.Nrow();i++){
		if ( fabs( _value[i]-dm.getValue()[i] )>epsilon) res = false;
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
		
	_dim = n;

	std::ifstream fichier(file,std::ios::in); 
	while (std::getline(fichier, ligne)){
		std::stringstream ss(ligne);
		double tmp;
		ss >> tmp;
		_value[i]=tmp;i++;
	}
}
	
private:

	int64_t _dim;
	double* _value;
	
};

class Matrix {
	
public:
	

Matrix(int nrow, int ncol) {
	_value = new Eigen::MatrixXd(nrow, ncol);
}

~Matrix() {
	if (_value)
		delete _value;
}

double* Store() {
	return _value->data();
}
	
double* GetRow(int index) {
	return _value->data() + index * _value->cols();
}
	
int Nrow() {
	return _value->rows();
}
	
int Ncol() {
	return _value->cols();
}

Eigen::MatrixXd* getValue(){
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
			if ( fabs( fabs((_value)->data()[i*ncol+j])-fabs(m.getValue()->data()[i*ncol+j]) )>epsilon) res = false;
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
	_value->resize(n,n);

	std::ifstream fichier(file,std::ios::in); 
	while (std::getline(fichier, ligne)){
		int j=0;
		// Construction pour une ligne
		std::stringstream ss(ligne);
		// Extraction des données pour une ligne
		while(!ss.eof()){
			double tmp;
			ss >> tmp;
//Not (*_value)(i,j)=tmp; anymore to be consistent with change line 294 (Eigen does not stock value like NEWMAT, U is actually transpose(eigenU))
			(*_value)(j,i)=tmp;
			j++;
			if (j==n) break;
		}
	i++;
	}
}
	
private:
	
	Eigen::MatrixXd* _value;
	
};

class SymmetricMatrix {
	
public:
	
		
// nrow == ncol
SymmetricMatrix(int nrow) {
	_value = new Eigen::MatrixXd(nrow, nrow);
}
	
SymmetricMatrix(int nrow,double* store){
	_value = new Eigen::MatrixXd(nrow, nrow);
	updateData(store);
}
	
~SymmetricMatrix() {
	if (_value)
		delete _value;
  if (_store)
    delete[] _store;
}
	
Eigen::MatrixXd* getValue(){
	return _value;
}

int Nrow(){
	return _value->rows();
}


	
void setValue(Eigen::MatrixXd const& matEigen){
	*_value=matEigen;
}
	
double* Store() {
	double* data = _value->data();
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
	
void updateData(double* store) {
	int i,j;
	int z=0;
	int ncol = _value->cols();
	//Remplissage de la partie supérieure de la matrice
	for(j=0; j<ncol ;j++){
		for(i=0; i<j+1 ; i++) {
			(*_value)(i,j)=store[z];
			(*_value)(j,i)=store[z];
			z++;
		}
	}
	//(*_value) = (*_value).triangularView<Eigen::Upper>();
}
	
// get determinant
double determinant(double* store) {
	updateData(store);
  double det = _value->determinant();
  return det;
}
	
// get inverse
SymmetricMatrix* Inverse(double* store) {
	//TODO Chercher à optimiser le nombre d'objet
	updateData(store);
	int _s_pbDimension = _value->rows();
  SymmetricMatrix* invMat = new SymmetricMatrix(_s_pbDimension);
  Eigen::MatrixXd  eigenInverse = _value->inverse();
  invMat->setValue(eigenInverse);
  return invMat;
}
	
// compute SVD (only matrices U and D, not V)
void computeSVD(DiagonalMatrix* D, Matrix* U,double* store) {
	// TODO: Eigen::ComputeThinU ?
  updateData(store);

//Switch when n > 16
  if (_value->rows() < 16) {
    auto svd = new Eigen::JacobiSVD<Eigen::MatrixXd>(*_value, Eigen::ComputeFullU);
    auto eigenD = svd->singularValues();
    auto eigenU = svd->matrixU();
    for (int64_t i=0; i<eigenD.rows(); i++) 
      D->Store()[i] = eigenD.data()[i];
    //for (int64_t i=0; i<eigenU.rows()*eigenU.cols(); i++)
    //  U->Store()[i] = eigenU.data()[i];
    //Eigen does not stock value like NEWMAT, U is actually transpose(eigenU)
    int64_t compt=0;
    for (int64_t i=0; i<eigenU.rows(); i++) {
      for (int64_t j=0; j<eigenU.cols(); j++) {
        U->Store()[compt] = eigenU.data()[i+j*eigenU.cols()];
        compt++;
      }
    }
    delete svd;

  }
  else {
    auto svd = new Eigen::BDCSVD<Eigen::MatrixXd>(*_value, Eigen::ComputeFullU);
    auto eigenD = svd->singularValues();
    auto eigenU = svd->matrixU();
    for (int64_t i=0; i<eigenD.rows(); i++) 
      D->Store()[i] = eigenD.data()[i];
    //for (int64_t i=0; i<eigenU.rows()*eigenU.cols(); i++)
    //  U->Store()[i] = eigenU.data()[i];
    //Eigen does not stock value like NEWMAT, U is actually transpose(eigenU)
    int64_t compt=0;
    for (int64_t i=0; i<eigenU.rows(); i++) {
      for (int64_t j=0; j<eigenU.cols(); j++) {
        U->Store()[compt] = eigenU.data()[i+j*eigenU.cols()];
        compt++;
      }
    }
    delete svd;

  } 
}


// Fonction pour les tests unitaires, compare les termes de 2 matrices
bool isAlmostEqual(SymmetricMatrix & sm,double epsilon) const {
	bool res=true;
	int ncol = _value->cols();
	for(int j=0; j<ncol ;j++){
		for(int i=0; i<j+1 ; i++){
			if ( fabs( (_value)->data()[i*ncol+j] -sm.getValue()->data()[i*ncol+j] )>epsilon) res = false;
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
		
	_value->resize(n,n);

	std::ifstream fichier(file,std::ios::in); 
	while (std::getline(fichier, ligne)){
		int j=0;
		// Construction pour une ligne
		std::stringstream ss(ligne);
		// Extraction des données pour une ligne
		while(!ss.eof()){
			double tmp;
			ss >> tmp;
			(*_value)(i,j)=tmp;
			j++;
			if (j==n) break;
		}
	i++;
	}
}
	
private:
	
	Eigen::MatrixXd* _value;
	double* _store;
};


}
}


#endif
