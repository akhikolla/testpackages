#ifndef GenAlgPLS_Chromosome_h
#define GenAlgPLS_Chromosome_h

#include "config.h"

#ifdef HAVE_CLIMITS
#include <climits>
#elif HAVE_LIMITS_H
#include <limits.h>
#endif

#include <vector>
#include <exception>
#include <iostream>
#include <RcppArmadillo.h>

#include "Control.h"
#include "TruncatedGeomGenerator.h"
#include "ShuffledSet.h"
#include "RNG.h"

class InvalidCopulationException : public Rcpp::exception {

public:
	InvalidCopulationException(const char *file, const int line) : Rcpp::exception("The two chromosomes are not compatible for mating", file, line) {};
};

class Chromosome {

public:
	Chromosome(const Control &ctrl, ShuffledSet &shuffledSet, RNG& rng, bool randomInit = true);
	Chromosome(const Chromosome &other, bool copyChromosomeParts = true);
//	~Chromosome();

	/**
	 * Re-initialize the chromosome randomly
	 */
	void randomlyReset(RNG& rng, ShuffledSet &shuffledSet);

	/**
	 * @return bool Returns true if mutation occurred, false otherwise
	 */
	bool mutate(RNG& rng);
	void mateWith(const Chromosome &other, RNG& rng, Chromosome& child1, Chromosome& child2);

	void setFitness(double fitness) { this->fitness = fitness; };
	double getFitness() const { return this->fitness; };

	Rcpp::LogicalVector toLogicalVector() const;
	arma::uvec toColumnSubset() const;

	bool isFitterThan(const Chromosome &ch) const;

	bool operator==(const Chromosome &ch) const;
	bool operator!=(const Chromosome &ch) const;
	Chromosome& operator=(const Chromosome &ch);
	
	uint16_t getVariableCount() const { return this->currentlySetBits; };

	friend std::ostream& operator<<(std::ostream &os, const Chromosome &ch);
private:
	static const uint8_t BITS_PER_PART = sizeof(IntChromosome) * BITS_PER_BYTE;

#ifdef INT_CHROMOSOME_MAX_VAL
	static const IntChromosome INT_CHROMOSOME_MAX = INT_CHROMOSOME_MAX_VAL;
#else
	static IntChromosome INT_CHROMOSOME_MAX;
	static IntChromosome getIntChromosomeMax();
#endif
	const Control &ctrl;
	const TruncatedGeomGenerator rtgeom;

	uint16_t numParts;
	uint16_t unusedBits;
	uint16_t currentlySetBits;
	
	/*
	 * Array with the chromosome parts
	 * If not all bits are used, the k least significant bits of the 1st(!) part are not used (k = this->unusedBits)
	 */
	std::vector<IntChromosome> chromosomeParts;
	double fitness;

//	std::vector<uint16_t> shuffledSet(uint16_t setSize, uint16_t shuffleSize, RNG& rng) const;
	
	/*
	 * Init the internal used chromosome parts completely random taking
	 * the minimum and maximum number of set bits specified by the
	 * control object into account
	 */
	void initChromosomeParts(RNG& rng, ShuffledSet &shuffledSet);

	/*
	 * The RNG only returns 32 random bits
	 * so two random numbers must be "glued" together to form
	 * a 64bit random number
	 */
	IntChromosome rbits(RNG& rng) const;

	std::ostream& printBits(std::ostream &os, IntChromosome bits, uint16_t leaveOut = 0) const;

	inline void updateCurrentlySetBits();

	/*
	 * count trailing zeros
	 */
	inline uint16_t ctz(IntChromosome mask) const;

	void copyFrom(const Chromosome& ch, bool copyChromosomeParts);

#if !(defined HAVE_BUILTIN_POPCOUNTLL | defined HAVE_BUILTIN_POPCOUNTL)
	static const IntChromosome M1 = 0x5555555555555555; // binary: 010101010101... (1 zero, 1 one)
	static const IntChromosome M2 = 0x3333333333333333; // binary: 001100110011... (2 zeros, 2 ones)
	static const IntChromosome M4 = 0x0f0f0f0f0f0f0f0f; // binary: 000011110000... (4 zeros, 4 ones)
	static const IntChromosome H01 = 0x0101010101010101; // the sum of 256 to the power of 0,1,2,3...

	static uint16_t popcount(IntChromosome x);
#endif

};

#endif
