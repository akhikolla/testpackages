//
//  Chromosome.cpp
//
//

#include "config.h"

#include <exception>
#include <algorithm>
#include <vector>
#include <iostream>
#include <set>
#include <RcppArmadillo.h>

#include "Logger.h"
#include "Chromosome.h"
#include "GenAlg.h"

#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#endif

#ifdef TIMING_BENCHMARK
#include <sys/time.h>
#endif

using namespace Rcpp;

#ifdef ENABLE_DEBUG_VERBOSITY
#define IF_DEBUG(expr) if(this->ctrl.verbosity == DEBUG_GA || this->ctrl.verbosity == DEBUG_ALL) { expr; }
#else
#define IF_DEBUG(expr)
#endif

#ifndef INT_CHROMOSOME_MAX_VAL
IntChromosome Chromosome::getIntChromosomeMax() {
	IntChromosome max = 0xFF; // has at least 1 byte!
	size_t remainingBytes = sizeof(IntChromosome) - 1;

	for(size_t i = 0; i < remainingBytes; ++i) {
		max <<= 8; // shift 8 bits (1 byte) to the left
		max += 0xFF; // set first byte
	}
	
	return max;
}

IntChromosome Chromosome::INT_CHROMOSOME_MAX = Chromosome::getIntChromosomeMax();
#endif

Chromosome::Chromosome(const Control &ctrl, ShuffledSet &shuffledSet, RNG& rng, bool randomInit) : ctrl(ctrl), rtgeom(1 - ctrl.mutationProbability), currentlySetBits(0) {
	// Determine the number of IntChromosome bit values that are
	// needed to represent all genes
	
	this->numParts = (uint16_t) this->ctrl.chromosomeSize / Chromosome::BITS_PER_PART;

	int rest = this->ctrl.chromosomeSize % Chromosome::BITS_PER_PART;

	if(rest > 0) {
		this->numParts += 1;
		this->unusedBits = Chromosome::BITS_PER_PART - rest;
	} else {
		this->unusedBits = 0;
	}

	this->fitness = 0.0;

	// Initialize chromosome parts randomly
	this->chromosomeParts.resize(this->numParts, 0);
	if(randomInit == true) {
		this->initChromosomeParts(rng, shuffledSet);
	}
}

Chromosome::Chromosome(const Chromosome &other, bool copyChromosomeParts) : ctrl(other.ctrl), rtgeom(other.rtgeom) {
	this->copyFrom(other, copyChromosomeParts);
}


void Chromosome::copyFrom(const Chromosome &other, bool copyChromosomeParts) {
	this->fitness = other.fitness;
	this->numParts = other.numParts;
	this->unusedBits = other.unusedBits;
	
	// Copy chromosome parts
	if(copyChromosomeParts) {
		this->chromosomeParts = other.chromosomeParts;
		this->currentlySetBits = other.currentlySetBits;
	} else {
		this->chromosomeParts.resize(this->numParts, 0);
		this->currentlySetBits = 0;
	}
}

void Chromosome::randomlyReset(RNG& rng, ShuffledSet &shuffledSet) {
	this->fitness = 0.0;
	this->initChromosomeParts(rng, shuffledSet);
}

inline void Chromosome::initChromosomeParts(RNG& rng, ShuffledSet &shuffledSet) {
	this->currentlySetBits = rng(this->ctrl.minVariables, this->ctrl.maxVariables + 1);

	ShuffledSet::iterator setPosIter = shuffledSet.shuffle(rng);
	ShuffledSet::iterator end = setPosIter + this->currentlySetBits;

	uint16_t randPos;
	uint16_t part;
	uint16_t offset;

	/*
	 * Reset chromosome to 0
	 */
	std::fill(this->chromosomeParts.begin(), this->chromosomeParts.end(), 0);
	
	for(; setPosIter != end; ++setPosIter) {
		randPos = this->unusedBits + (*setPosIter);
		part = randPos / Chromosome::BITS_PER_PART;
		offset = randPos % Chromosome::BITS_PER_PART;
		this->chromosomeParts[part] |= (((IntChromosome) 1) << offset);
	}
	IF_DEBUG(
		GAout << "Initialized chromosome with " << this->currentlySetBits << " bits set" << std::endl;
	)
}

void Chromosome::mateWith(const Chromosome &other, RNG& rng, Chromosome& child1, Chromosome& child2) {
	if(other.ctrl.chromosomeSize != this->ctrl.chromosomeSize) {
		throw InvalidCopulationException(__FILE__, __LINE__);
	}
	
	if(child1.chromosomeParts.size() != this->numParts) {
		child1.chromosomeParts.resize(this->numParts, 0);
	}

	if(child2.chromosomeParts.size() != this->numParts) {
		child2.chromosomeParts.resize(this->numParts, 0);
	}

	switch(this->ctrl.crossover) {
		case RANDOM: {
			IntChromosome randomMask = 0;
			IntChromosome negRandomMask = 0;
			for(uint16_t i = 0; i < this->numParts; ++i) {
				if(this->chromosomeParts[i] == other.chromosomeParts[i]) {
					/*
					 * This part is the same in both parents -- just copy it to the children
					 */
					child1.chromosomeParts[i] = child2.chromosomeParts[i] = this->chromosomeParts[i];
					IF_DEBUG(
						GAout << GAout.lock() << "Chromosome part is the same for both parents -- copy part to both children" << std::endl << GAout.unlock();
					)
				} else {
					/*
					 * Randomly pick some bits from one chromosome and some bits from the other
					 * chromosome
					 */
					do {
						randomMask = this->rbits(rng);
						if(i == 0) {
							randomMask <<= this->unusedBits;
						}
					} while(randomMask == 0);
					negRandomMask = ~randomMask;
					
					child1.chromosomeParts[i] = (this->chromosomeParts[i] & randomMask) | (other.chromosomeParts[i] & negRandomMask);
					child2.chromosomeParts[i] = (this->chromosomeParts[i] & negRandomMask) | (other.chromosomeParts[i] & randomMask);
					
					IF_DEBUG(
						GAout << GAout.lock() << "Mask for part " << i << ": ";
						this->printBits(GAout, randomMask, (i == 0) ? this->unusedBits : 0) << std::endl << GAout.unlock()
					)
				}
			}
			break;
		}
		case SINGLE:
		default: {
			/*
			 * Single crossover - draw a random position
			 */
			uint16_t randPos = (uint16_t) rng(1.0, this->ctrl.chromosomeSize);
			uint16_t chosenPart = randPos / Chromosome::BITS_PER_PART;
			uint16_t crossoverBit = randPos % Chromosome::BITS_PER_PART;
			IntChromosome coMask = (INT_CHROMOSOME_MAX >> crossoverBit);

			IF_DEBUG(
				GAout << GAout.lock() << "Crossover at position " << randPos << " (= part " << chosenPart << ", bit " << crossoverBit << ")" << std::endl << GAout.unlock()
			)
					
			std::copy(this->chromosomeParts.begin(), this->chromosomeParts.begin() + chosenPart, child1.chromosomeParts.begin());
			std::copy(other.chromosomeParts.begin(), other.chromosomeParts.begin() + chosenPart, child2.chromosomeParts.begin());
			
			child1.chromosomeParts[chosenPart] = ((this->chromosomeParts[chosenPart] & (~coMask)) | (other.chromosomeParts[chosenPart] & coMask));
			child2.chromosomeParts[chosenPart] = ((other.chromosomeParts[chosenPart] & (~coMask)) | (this->chromosomeParts[chosenPart] & coMask));
			
			std::copy(other.chromosomeParts.begin() + chosenPart + 1, other.chromosomeParts.end(), child1.chromosomeParts.begin() + chosenPart + 1);
			std::copy(this->chromosomeParts.begin() + chosenPart + 1, this->chromosomeParts.end(), child2.chromosomeParts.begin() + chosenPart + 1);
			
			break;
		}
	}

	IF_DEBUG(
		GAout << GAout.lock() << "1st child: " << child1 << std::endl
		<< "2nd child: " << child2 << std::endl << GAout.unlock();
	)

	child1.updateCurrentlySetBits();
	child2.updateCurrentlySetBits();
}

bool Chromosome::mutate(RNG& rng) {
	if(this->ctrl.mutationProbability == 0.0) {
		return false;
	}
	
	uint16_t currentlyUnsetBits = this->ctrl.chromosomeSize - this->currentlySetBits;
	ShuffledSet shuffledSet;

	int32_t numChangeBits = 0; // a negative value means unsetting bits and a positive value means setting new bits

	if((this->ctrl.minVariables - this->currentlySetBits) > 0) { // Too few bits are set -- we MUST add some
		numChangeBits = (this->ctrl.minVariables - this->currentlySetBits);		
	} else if((this->ctrl.maxVariables - this->currentlySetBits) < 0) { // Too many bits are set -- we MUST unset some
		numChangeBits = (this->ctrl.maxVariables - this->currentlySetBits);
	}

	IF_DEBUG(
		if(numChangeBits != 0) {
			GAout << GAout.lock() << "The number of set bits (" << this->currentlySetBits << ") is not within the valid range. MUST change " << numChangeBits << std::endl << GAout.unlock();
		}
	)

	if(this->ctrl.maxVariables - this->currentlySetBits > 0) {
		/*
		 * We may set some additional bits
		 * numChangeBits is 0 if the number of currently set bits is in the predefined range
		 * or *positive* if too few bits are set
		 */
		numChangeBits += this->rtgeom(this->ctrl.maxVariables - this->currentlySetBits - numChangeBits, rng);
	}

	if(this->currentlySetBits - this->ctrl.minVariables > 0) {
		/*
		 * We may unset some bits
		 * numChangeBits is 0 if the number of currently set bits is in the predefined range
		 * or *negative* if too many bits are set
		 */
		numChangeBits -= this->rtgeom(this->currentlySetBits - this->ctrl.minVariables + numChangeBits, rng);
	}

	IF_DEBUG(
		GAout << GAout.lock();
		if(numChangeBits != 0) {
			GAout << "Changing " << numChangeBits << " bits" << std::endl;
		} else {
			GAout << "No mutation" << std::endl;
		}
		GAout << GAout.unlock();
	)

	if(numChangeBits == 0) {
		return false;
	} else if(numChangeBits < 0) {
		/*
		 * We have to unset -numChangeBits bits, i.e. flip from 1 to 0
		 */
		ShuffledSet::iterator shuffledIt = shuffledSet.shuffle(this->currentlySetBits, rng, (numChangeBits == -1));
		std::vector<uint16_t> removeBitsPos(shuffledIt, shuffledIt + ((uint16_t) -numChangeBits));
		std::sort(removeBitsPos.begin(), removeBitsPos.end());
		std::vector<uint16_t>::iterator removeBitsPosIt = removeBitsPos.begin();

		IntChromosome mask = 0;
		IntChromosome tmp;
		int8_t totalShift = this->unusedBits - 1; // 0 based -- position 0 is the least significant bit
		int8_t trailingZeros = totalShift; // 1 based -- the value 1 means 1 trailing zero
		
		uint16_t onesCount = 0;
		
		for(uint16_t i = 0; (i < this->numParts) && (removeBitsPosIt != removeBitsPos.end()); ++i) {
			tmp = this->chromosomeParts[i];
			do {
				tmp >>= trailingZeros + 1;
				
				trailingZeros = this->ctz(tmp); // get number of 0's before the first 1 (starting with least significant bit)
				
				if(trailingZeros >= Chromosome::BITS_PER_PART - totalShift - 1) {
					trailingZeros = Chromosome::BITS_PER_PART - totalShift - 1;
				}
				
				totalShift += trailingZeros + 1;
				
				if(totalShift < Chromosome::BITS_PER_PART) {
					if(onesCount == (*removeBitsPosIt)) {
						//unset bit
						mask |= ((IntChromosome) 1) << totalShift;
						++removeBitsPosIt;
						
						if(removeBitsPosIt == removeBitsPos.end()) {
							break; // All bits have been removed -- exit while
						}
					}
					++onesCount; // increment after comparision because the position is 0 based
				}
			} while(totalShift < Chromosome::BITS_PER_PART);
			
			this->chromosomeParts[i] ^= mask; // toggle all bits set in the mask
			
			totalShift = -1;
			trailingZeros = -1;
			mask = 0;
		}		
	} else { // (numChangeBits > 0)
		/*
		 * We have to set numChangeBits bits, i.e. flip from 0 to 1
		 */
		ShuffledSet::iterator shuffledIt = shuffledSet.shuffle(currentlyUnsetBits, rng, (numChangeBits == 1));
		std::vector<uint16_t> addBitsPos(shuffledIt, shuffledIt + numChangeBits);
		std::sort(addBitsPos.begin(), addBitsPos.end());
		std::vector<uint16_t>::iterator addBitsPosIt = addBitsPos.begin();

		IntChromosome mask = 0;
		IntChromosome tmp;
		int8_t totalShift = this->unusedBits - 1; // 0 based -- position 0 is the least significant bit
		int8_t trailingZeros = totalShift; // 1 based -- the value 1 means 1 trailing zero
		
		uint16_t onesCount = 0, zerosCount = 0;
		
		for(uint16_t i = 0; (i < this->numParts) && (addBitsPosIt != addBitsPos.end()); ++i) {
			tmp = this->chromosomeParts[i];
			do {
				tmp >>= trailingZeros + 1;
				
				trailingZeros = this->ctz(tmp); // get number of 0's before the first 1 (starting with least significant bit)
				
				if(trailingZeros >= Chromosome::BITS_PER_PART - totalShift - 1) {
					trailingZeros = Chromosome::BITS_PER_PART - totalShift - 1;
				}
				
				totalShift += trailingZeros + 1;
				
				// There are trailingZeros 0-bits between the last found 1 and this 1
				// Some of them may be toggled
				zerosCount += trailingZeros;
				while(addBitsPosIt != addBitsPos.end() && zerosCount > (*addBitsPosIt)) {
					mask |= (((IntChromosome) 1) << (((*addBitsPosIt) + onesCount + this->unusedBits) - i * Chromosome::BITS_PER_PART));
					++addBitsPosIt;
				}
				
				if(totalShift < Chromosome::BITS_PER_PART) {
					++onesCount;
				}
			} while(totalShift < Chromosome::BITS_PER_PART);
			
			this->chromosomeParts[i] ^= mask; // toggle all bits set in the mask
			
			totalShift = -1;
			trailingZeros = -1;
			mask = 0;
		}
	}
	
	this->currentlySetBits += numChangeBits;
	
	IF_DEBUG(
		GAout << GAout.lock() << "After mutation:\n" << *this << std::endl << GAout.unlock();
	)
	
	return true;
}


std::ostream& operator<<(std::ostream &os, const Chromosome &ch) {
	ch.printBits(os, ch.chromosomeParts[0], ch.unusedBits);

	for (uint16_t i = 1; i < ch.numParts; ++i) {
		os << ' ';
		ch.printBits(os, ch.chromosomeParts[i]);
	}

	return os;
}

inline std::ostream& Chromosome::printBits(std::ostream &os, IntChromosome bits, uint16_t leaveOut) const {
	IntChromosome mask = ((IntChromosome) 1) << leaveOut;
	uint8_t delCount = 0;

	do {
		os << (((bits & mask) > 0) ? '1' : '0');
		if((++delCount % DELIMITER_POSITION) == 0) {
			delCount = 0;
			os << ' ';
		}
		mask <<= 1;
	} while(mask > 0);

	return os;
}

Rcpp::LogicalVector Chromosome::toLogicalVector() const {
	LogicalVector varVector;
	IntChromosome mask = ((IntChromosome) 1) << this->unusedBits;

	for (uint16_t i = 0; i < this->numParts; ++i) {
		do {
			varVector.push_back((this->chromosomeParts[i] & mask) > 0);
			mask <<= 1;
		} while(mask > 0);
		mask = (IntChromosome) 1;
	}

	return varVector;
}

arma::uvec Chromosome::toColumnSubset() const {
	arma::uvec columnSubset(this->currentlySetBits);
	IntChromosome mask = ((IntChromosome) 1) << this->unusedBits;
	uint16_t csIndex = 0;
	arma::uword truePos = 0;

	for (uint16_t i = 0; i < this->numParts && csIndex < columnSubset.n_elem; ++i) {
		do {
			if((this->chromosomeParts[i] & mask) > 0) {
				columnSubset[csIndex++] = truePos;
			}
			++truePos;
			mask <<= 1;
		} while(mask > 0 && csIndex < columnSubset.n_elem);
		mask = (IntChromosome) 1;
	}

	return columnSubset;

}

bool Chromosome::operator==(const Chromosome &ch) const {
	return (this->chromosomeParts == ch.chromosomeParts);
}

bool Chromosome::operator!=(const Chromosome &ch) const {
	return !(*this == ch);
}

bool Chromosome::isFitterThan(const Chromosome &ch) const {
	if(this->fitness > ch.fitness) {
		return true;
	} else if(this->fitness < ch.fitness) {
		return false;
	}

	// Both chromosomes have (almost) same fitness
	// so check if this chromosome has less bits set than the other chromosome
	return (this->currentlySetBits < ch.currentlySetBits);
}

Chromosome& Chromosome::operator=(const Chromosome &ch) {
	if(this != &ch) {
		this->copyFrom(ch, true);
	}
	return *this;
}

inline void Chromosome::updateCurrentlySetBits() {
	this->currentlySetBits = 0;
	
#ifdef HAVE_BUILTIN_POPCOUNTLL
	for(uint16_t i = 0; i < this->numParts; ++i) {
		this->currentlySetBits += __builtin_popcountll(this->chromosomeParts[i]);
	}
#elif (defined HAVE_BUILTIN_POPCOUNTL && !(defined HAVE_UNSIGNED_LONG_LONG))
	for(uint16_t i = 0; i < this->numParts; ++i) {
		this->currentlySetBits += __builtin_popcountl(this->chromosomeParts[i]);
	}
#else
	if(sizeof(IntChromosome) > 8) {
		throw std::overflow_error("The 'popcount' fallback algorithm can not handle integers with more than 8 bytes");
	}

	for(uint16_t i = 0; i < this->numParts; ++i) {
		this->currentlySetBits += this->popcount(this->chromosomeParts[i]);
	};
#endif
}

#if !(defined HAVE_BUILTIN_POPCOUNTLL | defined HAVE_BUILTIN_POPCOUNTL)
uint16_t Chromosome::popcount(IntChromosome x) {
	if(x == 0) {
		return 0;
	}
	
	x -= (x >> 1) & Chromosome::M1;							// put count of each 2 bits into those 2 bits
	x = (x & Chromosome::M2) + ((x >> 2) & Chromosome::M2);	// put count of each 4 bits into those 4 bits
	x = (x + (x >> 4)) & Chromosome::M4;					// put count of each 8 bits into those 8 bits
	return ((x * Chromosome::H01) >> 56);					// adds left 8 bits of tmp + (tmp << 8) + (tmp << 16) + (tmp << 24) + ...
}
#endif

/**
 * returns BITS_PER_PART if no bit is set in mask
 * positions start at 0
 */
inline uint16_t Chromosome::ctz(IntChromosome mask) const {
	if(mask == 0) {
		return Chromosome::BITS_PER_PART;
	}

#if defined HAVE_GCC_CTZLL
	return __builtin_ctzll(mask);
#elif defined HAVE_FFSLL
	return  ffsll(mask) - 1;
#elif (defined HAVE_GCC_CTZL && !(defined HAVE_UNSIGNED_LONG_LONG))
	return __builtin_ctzl(mask);
#elif (defined HAVE_FFSL && !(defined HAVE_UNSIGNED_LONG_LONG))
	return  ffsl(mask) - 1;
#else
	// Simple implementation of the "find first set" problem without loop
	// Only valid for IntChromosome with a maximum of 64 bits!
	if(sizeof(mask) > 8) {
		throw std::overflow_error("The 'count trailing zeros' fallback algorithm can not handle integers with more than 8 bytes");
	}
	
	uint16_t index = 0;

	if (sizeof(mask) == 8) {
		if((mask & 0x00000000FFFFFFFF) == 0) {
			index += 32;
			mask >>= 32;
		}
	}

	if((mask & 0x0000FFFF) == 0) {
		index += 16;
		mask >>= 16;
	}
	if((mask & 0x000000FF) == 0) {
		index += 8;
		mask >>= 8;
	}
	if((mask & 0x0000000F) == 0) {
		index += 4;
		mask >>= 4;
	}
	if((mask & 0x00000003) == 0) {
		index += 2;
		mask >>= 2;
	}
	if((mask & 0x00000001) == 0) {
		index += 1;
	}
	
	return index;
#endif
}


inline IntChromosome Chromosome::rbits(RNG& rng) const {
	// Calculate the number of random numbers to combine for one IntChromosome
	static int8_t partsPerRand = (sizeof(IntChromosome) * BITS_PER_BYTE / RNG::RANDOM_BITS);
	
	if(partsPerRand > 1) {
		IntChromosome rand = 0;

		for(int8_t i = partsPerRand; i >= 0; --i) {
			rand |= ((IntChromosome) (rng() << (i * RNG::RANDOM_BITS)));
		}

		return rand;
	} else {
		return (IntChromosome) rng();
	}
}
