//
//  OnlineStddev.h
//  gaselect
//
//  Created by David Kepplinger on 09.08.2014.
//
//

#ifndef gaselect_OnlineStddev_h
#define gaselect_OnlineStddev_h

#include <vector>
#include <algorithm>


/**
 * Helper class for online calculation of the standard deviation (and optionally the sum)
 */
class OnlineStddev {
public:
	OnlineStddev(uint16_t dim = 1) : dim(dim),
		meanVec(dim, 0.0), M2(dim, 0.0), counter(dim, 0) {}

	inline void reset() {
		std::memset(&(this->counter[0]), 0, this->dim * sizeof(this->counter[0]));
		std::memset(&(this->meanVec[0]), 0.0, this->dim * sizeof(this->meanVec[0]));
		std::memset(&(this->M2[0]), 0.0, this->dim * sizeof(this->M2[0]));
	};

	inline void update(arma::vec samples, uint16_t dim = 0) {
		for(uint16_t i = 0; i < samples.n_elem; ++i) {
			this->update(samples[i], dim);
		}
	};

	inline void update(double sample, uint16_t dim = 0) {
		double delta = sample - this->meanVec[dim];
		this->meanVec[dim] += delta / (++this->counter[dim]);
		this->M2[dim] += delta * (sample - this->meanVec[dim]);
	};

	inline double stddev(uint16_t dim = 0) const {
		return std::sqrt(this->M2[dim] / (this->counter[dim] - 1));
	};

	inline double var(uint16_t dim = 0) const {
		return this->M2[dim] / (this->counter[dim] - 1);
	};

	inline uint16_t N(uint16_t dim = 0) const {
		return this->counter[dim];
	}

	inline double mean(uint16_t dim = 0) const {
		return this->meanVec[dim];
	};

private:
	const uint16_t dim;
	std::vector<double> meanVec;
	std::vector<double> M2;
	std::vector<uint16_t> counter;
};
#endif
