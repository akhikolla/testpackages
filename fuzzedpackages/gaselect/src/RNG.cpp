/******************************
 *
 * C++ implementation of the WELL19937a random number generator
 * as introduced by:
 * F. Panneton, P. L'Ecuyer and M. Matsumoto. "Improved Long-Period Generators Based on Linear Recurrences Modulo 2"
 * (http://www.iro.umontreal.ca/~lecuyer/myftp/papers/wellrng.pdf)
 *
 * Used the original C code from http://www3.ocn.ne.jp/~harase/megenerators.html
 *
 ******************************/
#include "RNG.h"
#include <vector>
#include <stdexcept>

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define MAT1(v) v
#define MAT3POS(t,v) (v>>t)

#define V0            this->STATE[this->stateIndex]
#define VM1Over       this->STATE[this->stateIndex + RNG::M1 - RNG::R]
#define VM1           this->STATE[this->stateIndex + RNG::M1]
#define VM2Over       this->STATE[this->stateIndex + RNG::M2 - RNG::R]
#define VM2           this->STATE[this->stateIndex + RNG::M2]
#define VM3Over       this->STATE[this->stateIndex + RNG::M3 - RNG::R]
#define VM3           this->STATE[this->stateIndex + RNG::M3]
#define VRm1          this->STATE[this->stateIndex - 1]
#define VRm1Under     this->STATE[this->stateIndex + RNG::M3 - 1]
#define VRm2          this->STATE[this->stateIndex - 2]
#define VRm2Under     this->STATE[this->stateIndex + RNG::M3 - 2]

#define newV0         this->STATE[this->stateIndex - 1]
#define newV0Under    this->STATE[this->stateIndex - 1 + RNG::R]
#define newV1         this->STATE[this->stateIndex ]
#define newVRm1       this->STATE[this->stateIndex - 2]
#define newVRm1Under  this->STATE[this->stateIndex - 2 + RNG::R]

#define newVM2Over    this->STATE[this->stateIndex + RNG::M2 - RNG::R + 1]
#define newVM2        this->STATE[this->stateIndex + RNG::M2 + 1]

#define SEED_CONSTANT 1812433253U

/* tempering paramater */
#define BITMASK 0x41180000

const double RNG::RANDOM_MAX = 4294967296; // 2^W = 2^32

RNG::RNG(void) {
	this->stateIndex = 0;
	this->genFun = &RNG::case1;
}

RNG::RNG(uint32_t seed) {
	this->seed(seed);
}

RNG::RNG(const std::vector<uint32_t> &seed) {
	this->seed(seed);
}

void RNG::seed(const std::vector<uint32_t> &seed) {
	if(seed.size() < RNG::SEED_SIZE) {
		throw std::invalid_argument("The size of the seed must not be smaller than the RNG's seed size");
	}

	std::copy(seed.begin(), seed.begin() + RNG::SEED_SIZE, this->STATE);
	this->stateIndex = 0;
	this->genFun = &RNG::case1;

	for(int16_t i = RNG::BURNIN; i > 0; --i) {
		(this->*genFun)();
	}
}

void RNG::seed(uint32_t seed) {
	this->STATE[0] = seed;
	this->stateIndex = 0;

	// Same generator used to seed Mersenne twister
	uint32_t i = 1;
	for (; i < RNG::R; ++i) {
		this->STATE[i] = (SEED_CONSTANT * (this->STATE[i - 1] ^ (this->STATE[i - 1] >> (RNG::W - 2))) + i);
	}

	this->genFun = &RNG::case1;

	for(int16_t i = RNG::BURNIN; i > 0; --i) {
		(this->*genFun)();
	}
}

uint32_t RNG::case1(void) { // stateIndex = 0
	this->z0 = (VRm1Under & RNG::MASKL) | (VRm2Under & RNG::MASKU);
	this->z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
	this->z2 = MAT3POS (9, VM2) ^ MAT0POS (1, VM3);
	newV1 = this->z1 ^ this->z2;
	newV0Under = MAT1 (z0) ^ MAT0NEG (-9, this->z1) ^ MAT0NEG (-21, this->z2) ^ MAT0POS (21, newV1);
	this->stateIndex = RNG::R - 1;

	this->genFun = &RNG::case3;

	return (this->STATE[this->stateIndex] ^ (newVM2Over & BITMASK));
}

uint32_t RNG::case2(void) { // stateIndex = 1
	this->z0 = (VRm1 & RNG::MASKL) | (VRm2Under & RNG::MASKU);
	this->z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
	this->z2 = MAT3POS (9, VM2) ^ MAT0POS (1, VM3);
	newV1 = this->z1 ^ this->z2;
	newV0 = MAT1 (this->z0) ^ MAT0NEG (-9, this->z1) ^ MAT0NEG (-21, this->z2) ^ MAT0POS (21, newV1);
	this->stateIndex = 0;

	this->genFun = &RNG::case1;
	return (this->STATE[this->stateIndex] ^ (newVM2 & BITMASK));
}

uint32_t RNG::case3(void) { // stateIndex + RNG::M1 >= RNG::R
	this->z0 = (VRm1 & RNG::MASKL) | (VRm2 & RNG::MASKU);
	this->z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1Over);
	this->z2 = MAT3POS (9, VM2Over) ^ MAT0POS (1, VM3Over);
	newV1 = this->z1 ^ this->z2;
	newV0 = MAT1 (this->z0) ^ MAT0NEG (-9, this->z1) ^ MAT0NEG (-21, this->z2) ^ MAT0POS (21, newV1);
	this->stateIndex--;
	if (this->stateIndex + RNG::M1 < RNG::R) {
		this->genFun = &RNG::case5;
	}
	return (this->STATE[this->stateIndex] ^ (newVM2Over & BITMASK));
}

uint32_t RNG::case4(void) { // this->stateIndex + RNG::M3 >= RNG::R
	this->z0 = (VRm1 & RNG::MASKL) | (VRm2 & RNG::MASKU);
	this->z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
	this->z2 = MAT3POS (9, VM2) ^ MAT0POS (1, VM3Over);
	newV1 = this->z1 ^ this->z2;
	newV0 = MAT1 (this->z0) ^ MAT0NEG (-9, this->z1) ^ MAT0NEG (-21, this->z2) ^ MAT0POS (21, newV1);
	this->stateIndex--;
	if (this->stateIndex + RNG::M3 < RNG::R) {
		this->genFun = &RNG::case6;
	}
	return (this->STATE[this->stateIndex] ^ (newVM2 & BITMASK));
}

uint32_t RNG::case5(void) { // this->stateIndex + RNG::M2 >= RNG::R
	this->z0 = (VRm1 & RNG::MASKL) | (VRm2 & RNG::MASKU);
	this->z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
	this->z2 = MAT3POS (9, VM2Over) ^ MAT0POS (1, VM3Over);
	newV1 = this->z1 ^ this->z2;
	newV0 = MAT1 (this->z0) ^ MAT0NEG (-9, this->z1) ^ MAT0NEG (-21, this->z2) ^ MAT0POS (21, newV1);
	this->stateIndex--;
	if (this->stateIndex + RNG::M2 < RNG::R) {
		this->genFun = &RNG::case4;
	}
	return (this->STATE[this->stateIndex] ^ (newVM2Over & BITMASK));
}

uint32_t RNG::case6(void) { // 2 <= this->stateIndex <= (RNG::M3 - RNG::R - 1)
	this->z0 = (VRm1 & RNG::MASKL) | (VRm2 & RNG::MASKU);
	this->z1 = MAT0NEG (-25, V0) ^ MAT0POS (27, VM1);
	this->z2 = MAT3POS (9, VM2) ^ MAT0POS (1, VM3);
	newV1 = this->z1 ^ this->z2;
	newV0 = MAT1 (this->z0) ^ MAT0NEG (-9, this->z1) ^ MAT0NEG (-21, this->z2) ^ MAT0POS (21, newV1);
	this->stateIndex--;
	if (this->stateIndex == 1) {
		this->genFun = &RNG::case2;
	}
	return (this->STATE[this->stateIndex] ^ (newVM2 & BITMASK));
}
