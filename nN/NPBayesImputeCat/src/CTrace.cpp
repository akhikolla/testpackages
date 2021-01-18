/*
 * Copyright (C) 2014 Quanli Wang
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 */ 
#include <set>
#include "CData.h"
#include "CLcm.h"
#include "CTrace.h"

static const char* ParameterList[] = {"index","alpha","k_star","Nmis","nu","z","ImputedX","psi"};
static const int NParaList = 8;

CTrace::CTrace(CLcm* pm) {
	trace = NULL;
	m = pm;
  mnsize = 0;
}
CTrace::~CTrace(void)
{
	ClearTrace();
}

void CTrace::ClearTrace() {
	if (TracedParameters.size() > 0 && trace != NULL) {
		for (unsigned int i = 0; i < TracedParameters.size(); i++) {
			delete [] trace[i];
		}
		delete trace; trace = NULL;
	}
}

void CTrace::PrepareTrace() {
	mnindex = 0;
	ClearTrace();
	int nparas = (int)TracedParameters.size();
	if (nparas > 0 && mnsize > 0) {
		trace = new double*[nparas];
		for (int i = 0; i < nparas; i++) { //refactor this later
			int size =0;
			if (TracedParameters[i] == "alpha" ||  TracedParameters[i] == "k_star" || TracedParameters[i] == "Nmis" || TracedParameters[i] == "index") {
				size = 1;
			}
			if (TracedParameters[i] == "nu") {
				size = m->par->K;
			}
			if (TracedParameters[i] == "z") {
				size = m->par->n;
			}
			if (TracedParameters[i] == "ImputedX") {
				size = m->par->n*m->par->J;
			}
			if (TracedParameters[i] == "psi") {
				size = m->par->K * m->par->cumLevelsJ[m->par->J];
			}
			if (size > 0) {//always
				trace[i] = new double[mnsize * size];
			}
		}
	}
}

bool  CTrace::Trace(int index,int currentiter) {
	if (index >= mnsize && TracedParameters.size() > 0) {
		return false;
	}
	mnindex = index+1;
	for (int i = 0; i < (int)TracedParameters.size(); i++) { //refactor this later
		if (TracedParameters[i] == "index") {
			trace[i][index] = currentiter;
		}
		if (TracedParameters[i] == "alpha") {
			trace[i][index] = m->par->alpha;
		}
		if (TracedParameters[i] == "k_star") {
			trace[i][index] = m->par->k_star;
		}
		if (TracedParameters[i] == "Nmis") {
			trace[i][index] = m->par->Nmis;
		}
		if (TracedParameters[i] == "nu") {
			for (int k = 0; k < m->par->K; k++) {
				trace[i][index*m->par->K+k] = m->par->nuK[k];
			}
		}
		if (TracedParameters[i] == "z") {
			for (int k = 0; k < m->par->n; k++) {
				trace[i][index*m->par->n+k] = m->par->zI[k];
			}
		}
		if (TracedParameters[i] == "ImputedX") {
			int size = m->par->n * m->par->J;
			std::copy(m->par->xIJ[0],m->par->xIJ[0]+ size,trace[i] + index * size);
		}
		if (TracedParameters[i] == "psi") {
			int size = m->par->K * m->par->cumLevelsJ[m->par->J];
			std::copy(m->par->psiJKL[0],m->par->psiJKL[0]+ size,trace[i] + index * size);
		}
	}
	return true;
}

std::vector<std::string> CTrace::GetTracedList() {
	return TracedParameters;
}

void CTrace::SetTrace(std::vector< std::string > list_, int size) {
	mnsize = size;
	TracedParameters.clear();
	std::set<std::string> traceAble(ParameterList, ParameterList + NParaList);
	for (unsigned int i = 0; i < list_.size(); i++) {
		if (traceAble.find(list_[i]) != traceAble.end()) {
			TracedParameters.push_back(list_[i]);
		}
	}
}

std::vector<std::string> CTrace::GetParameterList() {
	std::vector<std::string> list_(ParameterList, ParameterList + NParaList);
	return list_;
}
