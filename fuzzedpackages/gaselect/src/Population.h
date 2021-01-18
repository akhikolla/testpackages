//
//  Population.h
//
//

#ifndef GenAlgPLS_Population_h
#define GenAlgPLS_Population_h

#include "config.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <utility>
#include <algorithm>

#include "Logger.h"
#include "RNG.h"
#include "Chromosome.h"
#include "Evaluator.h"
#include "Control.h"
#include "OnlineStddev.h"

#ifdef ENABLE_DEBUG_VERBOSITY
#define IF_DEBUG(expr) if(this->ctrl.verbosity == DEBUG_GA || this->ctrl.verbosity == DEBUG_ALL) { expr; }
#else
#define IF_DEBUG(expr)
#endif

class Population {
public:
	class ChromosomeComparator : public std::binary_function<Chromosome, Chromosome, bool> {
	public:
		bool operator() (const Chromosome &lhs, const Chromosome &rhs) const {
			return rhs.isFitterThan(lhs);
		}
	};

	typedef std::set<Chromosome, Population::ChromosomeComparator> SortedChromosomes;

protected:
	typedef std::vector<Chromosome*> ChVec;
	typedef std::vector<Chromosome*>::iterator ChVecIt;
	typedef std::vector<Chromosome*>::reverse_iterator ChVecRIt;

	static const uint8_t MAX_DISCARDED_SOLUTIONS_RATIO = 20;
	static const int16_t DEFAULT_SCALING_MEAN = 6;

	const Control& ctrl;
	::Evaluator& evaluator;
	const std::vector<uint32_t> &seed;

	SortedChromosomes elite;
	std::vector<double> currentGenFitnessMap;
	double minEliteFitness;
	bool interrupted;

private:
	OnlineStddev fitStats;
	ChVec currentGeneration;
	std::vector<double> fitnessHistory;

public:
	Population(const Control &ctrl, ::Evaluator &evaluator, const std::vector<uint32_t> &seed) :
		ctrl(ctrl), evaluator(evaluator), seed(seed), currentGenFitnessMap(ctrl.populationSize + ctrl.elitism, 0.0),
		interrupted(false) {
		this->currentGeneration.reserve(this->ctrl.populationSize + this->ctrl.elitism);

		this->minEliteFitness = 0.0;

		this->fitnessHistory.reserve(3 * this->ctrl.numGenerations);

		switch (this->ctrl.fitnessScaling) {
			case EXP:
				this->transformFitness = &Population::transformFitnessExp;
				break;
			default:
				this->transformFitness = &Population::transformFitnessNone;
				break;
		}
	}
	
	virtual ~Population() {
		for(ChVecIt it = this->currentGeneration.begin(); it != this->currentGeneration.end(); ++it) {
			delete *it;
		}
	}

	virtual void run() = 0;

	inline bool wasInterrupted() const {
		return this->interrupted;
	}

	inline const std::vector<double>& getFitnessEvolution() {
		return this->fitnessHistory;
	};

	inline SortedChromosomes getResult() const {
		SortedChromosomes result(this->elite);
		
		/*
		 * result.insert(begin(), end()) can not be used because we need to
		 * insert the VALUE and not the pointer to a chromosome
		 */
		for(ChVec::const_iterator it = this->currentGeneration.begin(); it != this->currentGeneration.end(); ++it) {
			result.insert(**it);
		}
		
		return result;
	}

private:
	double (Population::*transformFitness)(double&) const;

	inline double transformFitnessExp(double& fitness) const {
		return std::exp(fitness);
	}

	inline double transformFitnessNone(double& fitness) const {
		return fitness;
	}

protected:
	inline void initCurrentGeneration(ShuffledSet &shuffledSet, RNG &rng) {
		for(uint32_t i = this->ctrl.elitism + this->ctrl.populationSize; i > 0; --i) {
			this->currentGeneration.push_back(new Chromosome(this->ctrl, shuffledSet, rng, false));
		}
	}

	/**
	 * Update the current generation as well as the fitness map of the current generation
	 *
	 * @param const ChVec &newGeneration The newly generated generation that will be copied to the currentGeneration
	 * @param double minFitness The minimum fitness of the new generation
	 * @param bool updateElite Set to true if the elite should be updated as well
	 */
	inline double updateCurrentGeneration(const ChVec &newGeneration, double minFitness, bool first = false, bool updateElite = false) {
		uint16_t i = 0;
		double sumFitness = 0.0, fitt, fitMean, fitSD;

		if (first) {
			for(; i < this->ctrl.populationSize; ++i) {
				this->fitStats.update(newGeneration[i]->getFitness());
			}
		}

		fitMean = this->fitStats.mean();
		fitSD = this->fitStats.stddev();

		this->fitStats.reset();

		if(this->ctrl.elitism > 0 && this->elite.size() > 0 && minFitness > this->elite.rbegin()->getFitness()) {
			minFitness = this->elite.rbegin()->getFitness();
		}

		minFitness = (minFitness - fitMean) / fitSD;
		minFitness = (this->*transformFitness)(minFitness);

		IF_DEBUG(GAout << "Fitness map:\n")

		/*
		 * Add newly generated chromosomes to the current generation
		 * this may introduce some duplicates but the number is usually
		 * neglectable, as the number of chromosomes is usally much greater than
		 * the number of elite chromosomes
		 */
		for(i = 0; i < this->ctrl.populationSize; ++i) {
			if(updateElite == true) {
				this->addChromosomeToElite(*newGeneration[i]);
			}

			*(this->currentGeneration[i]) = *(newGeneration[i]);

			this->fitStats.update(this->currentGeneration[i]->getFitness());

			fitt = (this->currentGeneration[i]->getFitness() - fitMean) / fitSD;
			sumFitness += (this->*transformFitness)(fitt) - minFitness;

			this->currentGenFitnessMap[i] = sumFitness;

			IF_DEBUG(
				GAout << (std::stringstream() << std::fixed << std::setw(4) << i).rdbuf()
				<< TAB_DELIMITER << sumFitness << "\n";
			)
		}

		/*
		 * Add all chromosomes from the elite
		 */
		for(SortedChromosomes::iterator eliteIt = this->elite.begin(); eliteIt != this->elite.end(); ++eliteIt, ++i) {
			*(this->currentGeneration[i]) = *(eliteIt);

			this->fitStats.update(this->currentGeneration[i]->getFitness());

			fitt = (eliteIt->getFitness() - fitMean) / fitSD;
			sumFitness += (this->*transformFitness)(fitt) - minFitness;

			this->currentGenFitnessMap[i] = sumFitness;

			IF_DEBUG(
				GAout << (std::stringstream() << std::fixed << std::setw(4) << i).rdbuf()
				<< TAB_DELIMITER << sumFitness << "\n";
			)
		}

		IF_DEBUG(GAout << std::endl)

		this->fitnessHistory.push_back(this->elite.rbegin()->getFitness());
		this->fitnessHistory.push_back(this->fitStats.mean());
		this->fitnessHistory.push_back(this->fitStats.stddev());

		return sumFitness;
	}
	
	/**
	 * Pick a chromosome from the current generation at random
	 * where the probability to pick a chromosome is taken from
	 * the currentGenFitnessMap
	 */
	inline Chromosome* drawChromosomeFromCurrentGeneration(double rand) const {
		int imin = 0, imax = static_cast<int>(this->currentGenFitnessMap.size());
		int imid = 0;
		
		/*
		 * Search for the index in the currentGenFitnessMap
		 * that holds a value that is GREATER or EQUAL than
		 * the given random number
		 */
		while(imin < imax) {
			imid = (imax + imin) / 2;
			
			if(this->currentGenFitnessMap[imid] < rand) {
				imin = imid + 1;
			} else {
				imax = imid;
			}
		}
		
		IF_DEBUG(
			GAout << GAout.lock() << "Selected chromosome " << imin << " for mating (rand = " << rand << ")" << std::endl << GAout.unlock();
		);

		return this->currentGeneration[imin];
	};

#ifdef ENABLE_DEBUG_VERBOSITY
	static bool compEqual(Chromosome* c1, Chromosome* c2) {
		return ((*c1) == (*c2));
	}
	
	static bool compLT(const Chromosome* const c1, const Chromosome* const c2) {
		return c1->isFitterThan(*c2);
	}
	
	inline uint16_t countUniques() const {
		std::vector<Chromosome*> gen = this->currentGeneration;
		std::sort(gen.begin(), gen.end(), Population::compLT);
		return std::distance(gen.begin(), std::unique(gen.begin(), gen.end(), Population::compEqual));
	};
#endif

	inline void printCurrentGeneration() {
		int i = 0;
		for(ChVecIt it = this->currentGeneration.begin(); it != this->currentGeneration.end(); ++it) {
			GAout << (std::stringstream() << std::fixed << std::setw(4) << i++ << ": ").rdbuf();
			this->printChromosomeFitness(GAout, **it);
		}
		GAout << "\n" << std::endl;
	}

	inline std::ostream& printChromosomeFitness(std::ostream &os, Chromosome &ch) {
		os << (std::stringstream() << std::fixed << std::setw(WIDTH) << std::setprecision(PRECISION) << ch.getFitness()).rdbuf()
		<< TAB_DELIMITER << ch << std::endl;
		
		return os;
	};
	
	inline void addChromosomeToElite(Chromosome &ch) {
		/* Add chromosome to elite if better than the worst elite-chromosome */
		if(this->ctrl.elitism > 0 && (ch.getFitness() > this->minEliteFitness || this->elite.size() < this->ctrl.elitism)) {
			/* Insert a copy of the chromosome. If the chromosome is a duplicate, it is not inserted */
			this->elite.insert(ch);
			
			/*
			 * If the chromosome was inserted and the elite-set was already `full`
			 * the last worst chromosome is removed
			 */
			if(this->elite.size() > this->ctrl.elitism) {
				this->elite.erase(this->elite.begin());
			}
			
			this->minEliteFitness = this->elite.begin()->getFitness();
			
			IF_DEBUG(
				GAout << "Adding chromosome to elite. New minimum fitness for elite is " << this->minEliteFitness << std::endl;
			)
		}
	};
	
	inline static std::pair<bool, bool> checkDuplicated(ChVecIt begin, ChVecRIt rbegin, const ChVecIt &child1It, const ChVecRIt &child2It) {
		std::pair<bool, bool> duplicated(false, **child1It == **child2It);
		
		while(begin != child1It && (duplicated.first == false || duplicated.second == false)) {
			if(duplicated.first == false && (**begin == **child1It)) {
				duplicated.first = true;
				(**child1It).setFitness((**begin).getFitness());
			}
			if(duplicated.second == false && (**begin == **child2It)) {
				duplicated.second = true;
				(**child2It).setFitness((**begin).getFitness());
			}
			++begin;
		}
		while(rbegin != child2It && (duplicated.first == false || duplicated.second == false)) {
			if(duplicated.first == false && (**rbegin == **child1It)) {
				duplicated.first = true;
				(**child1It).setFitness((**rbegin).getFitness());
			}
			if(duplicated.second == false && (**rbegin == **child2It)) {
				duplicated.second = true;
				(**child2It).setFitness((**rbegin).getFitness());
			}
			++rbegin;
		}
		
		return duplicated;
	};

	class CompChromsomePtr : public std::unary_function<Chromosome*, bool> {
	public:
		CompChromsomePtr(const Chromosome* const ch) : ch(ch) {}
		bool operator()(const Chromosome* const ch) const {
			return ((*this->ch) == (*ch));
		}
		
	private:
		const Chromosome* const ch;
	};
};
#endif
