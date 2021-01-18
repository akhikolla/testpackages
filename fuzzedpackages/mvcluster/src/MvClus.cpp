/*
 * MvClus.cpp
 *
 *  Created on: Aug 20, 2015
 *      Author: Javon
 */

#include "MvClus.h"

MvClus::MvClus(const vector<mat> & datasets)
:m_datasets(datasets) {
	m_nView = m_datasets.size();
	m_nSample = m_datasets[0].n_rows;

	m_nFeat.set_size(m_nView);
	for (uint8_t i = 0; i < m_nView; i++) {
		m_nFeat[i] = m_datasets[i].n_cols;
	}
}

MvClus::~MvClus() {

}
