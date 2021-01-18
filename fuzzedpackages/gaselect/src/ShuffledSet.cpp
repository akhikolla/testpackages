//
//  ShuffledSet.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 12.05.2013.
//
//

#include "config.h"
#include "ShuffledSet.h"

#include <RcppArmadillo.h>
#include <stdexcept>
#include <algorithm>

ShuffledSet::ShuffledSet() {}

ShuffledSet::ShuffledSet(arma::uword size) {
	this->reset(size);
}

void ShuffledSet::reset() {
	this->reset(this->set.size());
}

void ShuffledSet::reset(arma::uword size) {
	this->set.resize(size);
	arma::uword i = 0, j = 1;
	for(; j < size; i += 2, j += 2) {
		this->set[i] = i;
		this->set[j] = j;
	}
	if(i < size) {
		this->set[i] = i;
	}
}

const arma::uvec& ShuffledSet::shuffleAll(RNG &rng) {
	for(arma::uword i = 0; i < this->set.size(); ++i) {
		std::swap(this->set[i], this->set[rng(i, this->set.size())]);
	}
	
	return this->set;
}

ShuffledSet::iterator ShuffledSet::shuffle(RNG &rng) {
	std::swap(this->set[0], this->set[rng(0.0, this->set.size())]);

	return ShuffledSet::iterator(*this, rng);
}

ShuffledSet::iterator ShuffledSet::shuffle(arma::uword size, RNG &rng, bool onlyOne) {
	if(onlyOne) {
		this->set.resize(1);
		this->set[0] = (arma::uword) rng(0.0, size);
		return ShuffledSet::iterator(*this, rng);
	} else {
		this->reset(size);
		return this->shuffle(rng);
	}
}

ShuffledSet::iterator& ShuffledSet::iterator::operator++() {
#ifdef SHUFFLED_SET_CHECK_ITERATOR_STATE
	if(this->shifted == true) {
		throw std::logic_error("The iterator has been shifted and is thus not valid incrementing");
	}
#endif
	++this->pos;
	std::swap(obj.set[this->pos], obj.set[this->rng(this->pos, obj.set.size())]);
	return *this;
}

arma::uword ShuffledSet::iterator::operator*() const {
#ifdef SHUFFLED_SET_CHECK_ITERATOR_STATE
	if(this->shifted == true) {
		throw std::logic_error("The iterator has been shifted and is thus not valid for accessing elements");
	}
#endif
	return this->obj.set[this->pos];
}

ShuffledSet::iterator ShuffledSet::iterator::operator+(const arma::uword &shift) {
	ShuffledSet::iterator newIt(this->obj, this->rng);
	newIt.pos = this->pos + shift;
#ifdef SHUFFLED_SET_CHECK_ITERATOR_STATE
	newIt.shifted = true;
#endif
	return newIt;
}


ShuffledSet::iterator& ShuffledSet::iterator::operator+=(const arma::uword &shift) {
	this->pos += shift;
#ifdef SHUFFLED_SET_CHECK_ITERATOR_STATE
	this->shifted = true;
#endif
	return *this;
}

bool ShuffledSet::iterator::operator==(const ShuffledSet::iterator &it) const {
	return ((it.pos == this->pos) && (&(it.obj) == &(this->obj)));
}

bool ShuffledSet::iterator::operator!=(const ShuffledSet::iterator &it) const {
	return ((it.pos != this->pos) || (&(it.obj) != &(this->obj)));
}
