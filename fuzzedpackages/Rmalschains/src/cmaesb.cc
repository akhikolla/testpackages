/**
 * Copyright 2008, Daniel Molina Cabrera <danimolina@gmail.com>
 * 
 * This file is part of software Realea
 * 
 * Realea is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Realea is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "cmaesb.h"

using namespace realea;
using namespace realea::internal;

CMAESBound::CMAESBound(IEvalReal *eval, DomainRealPtr domain) :
	m_eval(eval), m_domain(domain), m_avanzed(true), m_ndim(domain->getDimension()), 
	m_diagC(m_ndim), m_scale(m_ndim), m_weights(m_ndim), m_dfithist() {
    m_validfitval = 0;
    m_numapplied = 0;

    m_isactive = m_domain->isBound();    

    if (!m_isactive) {
	return;
    }

    m_weights = 0; // weights for bound penalty
    Real meanC = mean(m_diagC);
    m_scale = m_diagC/meanC;
}

void CMAESBound::setParam(int lambda, double mueff, ColumnVector &sigma, MyMatrix &C) {
   m_lambda = lambda;
   m_mueff = mueff;
   m_sigma = sigma;
   m_dfithist_size = 20+(3*m_ndim)/lambda;
   m_dfithist.push(1.0);

   DiagonalMatrix DiagC(m_ndim);
   DiagC <<C;
   copyColumn(DiagC, m_diagC);
}


MyReturnMatrix xintobounds(const ColumnVector &mat, DomainRealPtr domain, vector<int> *corrected) {
    Real val;
    int posi, i;
    int n= mat.Nrows();

    ColumnVector mat_corrected(n);

    for (i = 0; i < n; i++) {
       tReal lower, upper;
       domain->getValues(i, &lower, &upper);
       val = mat[i];
       posi = i;

       if (val < lower) {
	  val = lower;
	  if (corrected != NULL) {
	     corrected->push_back(posi);
	  }
       }
       else if (val > upper) {
	  val = upper;

	  if (corrected != NULL) {
	     corrected->push_back(posi);
	  }
       }

       mat_corrected[i] = val;
    }

    mat_corrected.Release();
    return mat_corrected;
}


void CMAESBound::evalSols(ColumnVector &xmean, MyMatrix &arx, MyMatrix &arxvalid, 
					RowVector &fitness_raw, RowVector &fitness_sel) {
   int lambda = arx.Ncols();
   int i;

   vector<int> v_ti;
   v_ti.reserve(m_ndim);

   // Incremento el contador
   m_numapplied += 1;
//   cout <<"CMAESBnd::eval the "  <<m_numapplied <<" solutions" <<endl;

    // If there is active 'cut' the evaluation
    if (!m_isactive) {
       arxvalid = arx;
    }
    else {
       for (i = 1; i <= lambda; ++i) {
	    ColumnVector arx_i = arx.Column(i);
	    arxvalid.Column(i) = xintobounds(arx_i, m_domain, &v_ti);
       }
	
    }

    for (i = 1; i <= lambda; ++i) {
       ColumnVector arxvalid_i = arxvalid.Column(i);
       fitness_raw[i-1] = m_eval->eval(arxvalid_i.Store(), m_ndim);
    }

    // At the beginning fitness_sel = fitness_raw
    fitness_sel = fitness_raw;

    // Cambio arx if there is active
    if (!m_isactive || !m_avanzed) {
	return;
    }

    // There are punished the solutions that are out of the boundaries
    int perc[2] = {25,75};

    // Get delta fitness values
    RowVector r_val = myprctile(fitness_raw, perc, 2);
    ColumnVector val(m_ndim);
    val = (r_val[1]-r_val[0]) / m_ndim / mean(m_diagC);
    ColumnVector pow2_sigma = pow2(m_sigma);
    val = DotDivide(val, pow2_sigma);
    Real maxval = max(val);

    if (maxval == 0) { // happens if all points are out of bounds
       val = min_positive(m_dfithist);
    }
    else if (m_validfitval == 0) { // first sinsible var
       m_dfithist.pop();
       m_validfitval = 1;
    }

    // Reset the histogram
    if (m_dfithist.size() == m_dfithist_size) {
       m_dfithist.pop();
    }

    Real valuemin = min_positive(val);
    m_dfithist.push(valuemin);

    RowVector dfithist(m_dfithist.size());
    copyRow(m_dfithist, dfithist);

    // Update weight and scale
    m_weights = DotVectors(m_scale, m_weights);

    Real mean_d = mean(log(m_diagC));
    m_scale << exp(0.1*mean_d)*pow(m_diagC, 0.9);

    // cout <<"bnd_scale After :\n" <<m_scale <<endl;
    m_scale = m_scale / exp(mean(log(m_scale))); // prod is 1 initially
    // cout <<"bnd_scale Rescale :\n" <<m_scale <<endl;

    m_weights = DivVectors(m_weights, m_scale);
    // cout <<"bnd_weights:\n" <<m_weights <<endl;

    // check that xmean is in the range
    v_ti.clear();
    ColumnVector tx = xintobounds(xmean, m_domain, &v_ti);

    // cout <<"tx : " <<tx <<endl;

    // Set initial weights
    if (m_iniphase) {
       // cout <<"initphase" <<endl;
       // cout <<"v_ti.size : " <<v_ti.size() <<endl;

       // First time the weigths are updated
       if (v_ti.size() ) {
	  // copy(v_ti.begin(), v_ti.end(), ostream_iterator<int> (cout,","));
    
	  // cout <<endl;
	  Real hist_median =  median(dfithist);
	  m_weights = 2.0002* hist_median;
	  // cout <<"bnd_weights init:\n" <<m_weights <<endl;
	  m_weights = DivVectors(m_weights, m_scale);


	  // cout <<"initphase: bnd_weights:\n" <<m_weights <<endl;
	  if (m_validfitval && m_numapplied > 2) {
	     m_iniphase = false;
	  }
       }
    }

    // Increase/decrease weights
    if  (v_ti.size() ) { // any coordinate of xmean out of bounds
       // judge distance of xmean to boundary
       tx = xmean - tx;
       Real umbral_out = 3*max(1.0,sqrt((double) m_ndim)/m_mueff);
       Real increm = pow(1.1, max(1.0, m_mueff/10.0/m_ndim));
       vector<int>::iterator item;
       int posi;

       for (item = v_ti.begin(); item != v_ti.end(); item++) {
	  posi = *item;

	  if (tx[posi] > umbral_out*m_sigma[posi]*sqrt(m_diagC[posi])) {
	     m_weights *= increm;
	  }
       }
    }

    // cout <<"calculando dif" <<endl;
    Matrix dif = arxvalid - arx;
    // cout <<"dif : \n" <<dif <<endl;

    // cout <<"pow2_m(dif): \n" <<pow2_m(dif) <<endl;
    RowVector arpenalty = m_weights.t() * pow2_m(dif);
    // cout << "arpenalty : " << arpenalty <<endl;
    // Aplica penalización
    // TODO: It is only valide in minimization
    fitness_sel = fitness_raw + arpenalty;
}

/**
 * Esta función permite obtener la matriz indicada manteniendo todas las variables en orden
 */
MyReturnMatrix xintobounds(const ColumnVector   &mat, Real *lower, Real *upper, vector<int> *corrected) {
    Real val;
    int posi, i;
    int n= mat.Nrows();

    ColumnVector mat_corrected(n);

    for (i = 0; i < n; i++) {
       val = mat[i];
       posi = i;

       if (val < lower[i]) {
	  val = lower[i];
	  if (corrected != NULL) {
	     corrected->push_back(posi);
	  }
       }
       else if (val > upper[i]) {
	  val = upper[i];

	  if (corrected != NULL) {
	     corrected->push_back(posi);
	  }
       }

       mat_corrected[i] = val;
    }

    mat_corrected.Release();
    return mat_corrected;
}



