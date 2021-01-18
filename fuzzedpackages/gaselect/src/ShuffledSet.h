//
//  ShuffledSet.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 12.05.2013.
//
//

#ifndef GenAlgPLS_VariablePositionPopulation_h
#define GenAlgPLS_VariablePositionPopulation_h

#include "config.h"

#include <RcppArmadillo.h>
#include <iterator>
#include <vector>

#include "RNG.h"

class ShuffledSet
{
public:
	ShuffledSet();
	ShuffledSet(arma::uword size);
//	~ShuffledSet();
	
	class iterator : public std::iterator<std::input_iterator_tag, arma::uword> {
	public:
		/**
		 * Attention: The first element in obj.set must already
		 * be a shuffled element!
		 */
		iterator(ShuffledSet &obj, RNG &rng) : obj(obj), rng(rng), pos(0) {
#ifdef SHUFFLED_SET_CHECK_ITERATOR_STATE
			this->shifted = false;
#endif
		};
		
		arma::uword operator*() const;
		iterator& operator++();
		bool operator==(const iterator &it) const;
		bool operator!=(const iterator &it) const;

		/**
		 * These methods invalidate the shuffle-process
		 * Use only for comparisons (all checks for misuse are disabled!)
		 */
		iterator operator+(const arma::uword &shift);
		iterator& operator+=(const arma::uword &shift);
	private:
		ShuffledSet &obj;
		RNG &rng;
		arma::uword pos;

#ifdef SHUFFLED_SET_CHECK_ITERATOR_STATE
		bool shifted;
#endif
	};

	/**
	 * Shuffle the whole set and return a reference to the shuffled set (the reference
	 * will be invalid once the ShuffledSet object is destroyed)
	 */
	const arma::uvec& shuffleAll(RNG &rng);
	
	/**
	 *
	 */
	iterator shuffle(RNG &rng);

	/**
	 * First reset the size of the set to `size`
	 *
	 * @param arma::uword size ... The new size of the set
	 * @param RNG rng ... The random number generator instance
	 * @param bool onlyOne ... If this is TRUE, only one element is guaranteed to be shuffled
	 */
	iterator shuffle(arma::uword size, RNG &rng, bool onlyOne = false);
	
	/**
	 * Reset the set to a sorted stated with `size` elements
	 * (i.e. have values 0, 1, 2, 3, ..., size - 1)
	 *
	 * @param arma::uword size ... The new size of the set
	 */
	void reset(arma::uword size);

	/**
	 * Reset the set to a sorted state with the same number of elements
	 * as before
	 */
	void reset();

private:
	arma::uvec set;
};

#endif
