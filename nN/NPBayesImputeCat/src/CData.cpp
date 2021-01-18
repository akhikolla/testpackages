/*
 * Copyright (C) 2007-2014 Daniel Manrique-Vallier
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 * Modified by Quanli Wang, 2014
 */ 
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <string>
#include "CData.h"
#include "margin_conditions.h"
#include "CArrayND.h"
using namespace std;

//Constructor
CData::CData(){
	//Constructor
}

//Destructor
CData::~CData(){
	free(this->levelsJ);
	free(this->cumLevelsJ);
	//free(this->x[0]);
	free(this->x);
	if (this->nZeroMC > 0) {
		//free(this->ZeroMC_IJ[0]);
		free(this->ZeroMC_IJ);
	}
}

void CData::Process_MC(){
	if (nZeroMC == 0) { return; }
	// takes ZeroMC_IJ and convet it to an exclusive list
	int size;
	int** newzeros = MC_to_MCpartition(ZeroMC_IJ, J, nZeroMC, levelsJ, size);
	//replace the old list
	//free(this->ZeroMC_IJ[0]);
	free(this->ZeroMC_IJ);
	this->ZeroMC_IJ = newzeros;
	this->nZeroMC = size;
}

void CData::SetData(int *x_flat, int J, int n, int* ZeroMC_flat, int nZeroMC, int *levels){
	//have to set: J, levels, L, x, n
	this->J = J;
	this->n = n;
	//Allocate arrays
	this->levelsJ = new int[J]; 
	this->cumLevelsJ = new int[J+1]; 
	this->x = (int**)CArrayND<int>::flat2arrayND(new int[J * n], sizeof(int),2, n, J);
	//Copy arrays
	memcpy(this->x[0], x_flat, sizeof(int) * J * n);
	memcpy(levelsJ, levels, J * sizeof(int));
	
	this->cumLevelsJ[0] = 0;
	this->L = 0;
	for (int i = 0; i < J; i++) {
		this->cumLevelsJ[i+1] = this->cumLevelsJ[i] + levelsJ[i];
		if (levelsJ[i] > this->L) {this->L = levelsJ[i];}
	}
	
	this->nZeroMC = nZeroMC;
	if (nZeroMC >0 ) {
		this->ZeroMC_IJ = (int**)CArrayND<int>::flat2arrayND( new int[nZeroMC * J], sizeof(int),2, nZeroMC, J); 
		memcpy(ZeroMC_IJ[0], ZeroMC_flat, sizeof(int) * nZeroMC * J);
		if(!check_disjoint_MC(ZeroMC_IJ, nZeroMC, J)){
			Process_MC();
		}
	}
}
void CData::UpdateX(std::vector<int>& x_flat) {
  std::copy(x_flat.begin(), x_flat.end(), this->x[0]); 
}

void CData::SetData(std::vector<int>& x_flat, int J, int n, std::vector<int>& ZeroMC_flat, int nZeroMC, std::vector<int>& levels){
	//have to set: J, levels, L, x, n
	this->J = J;
	this->n = n;
	//Allocate arrays
	this->levelsJ = new int[J]; 
	this->cumLevelsJ = new int[J+1]; 
	this->x = (int**)CArrayND<int>::flat2arrayND(new int[J * n], sizeof(int),2, n, J);

	//Copy arrays
	std::copy(x_flat.begin(), x_flat.end(), this->x[0]); 
	std::copy(levels.begin(), levels.end(), levelsJ); 

	//memcpy(this->x[0], &x_flat[0], sizeof(int) * J * n);
	//memcpy(levelsJ, &levels[0], J * sizeof(int));
	
	this->cumLevelsJ[0] = 0;
	this->L = 0;
	for (int i = 0; i < J; i++) {
		this->cumLevelsJ[i+1] = this->cumLevelsJ[i] + levelsJ[i];
		if (levelsJ[i] > this->L) {this->L = levelsJ[i];}
	}
	
	this->nZeroMC = nZeroMC;	
	if (nZeroMC > 0) {
		this->ZeroMC_IJ = (int**)CArrayND<int>::flat2arrayND( new int[nZeroMC * J], sizeof(int),2, nZeroMC, J); 
		std::copy(ZeroMC_flat.begin(), ZeroMC_flat.end(), ZeroMC_IJ[0]); 
		//memcpy(ZeroMC_IJ[0], ZeroMC_flat, sizeof(int) * nZeroMC * J);
		if(!check_disjoint_MC(ZeroMC_IJ, nZeroMC, J)){
			Process_MC();
		}
	}
}

