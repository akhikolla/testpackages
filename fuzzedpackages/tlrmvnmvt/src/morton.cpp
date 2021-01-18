#include <RcppEigen.h>
#include "morton.h"

using namespace std;
using namespace Eigen;

// "insert" a 0 bit after each of the 16 low bits of x
uint32_t Part1By1(uint32_t x)
{
  x &= 0x0000ffff;                  // x = ---- ---- ---- ---- fedc ba98 7654 3210
  x = (x ^ (x <<  8)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
  x = (x ^ (x <<  4)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
  x = (x ^ (x <<  2)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
  x = (x ^ (x <<  1)) & 0x55555555; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
  return x;
}

// inverse of Part1By1 - "delete" all odd-indexed bits
uint32_t Compact1By1(uint32_t x)
{
  x &= 0x55555555;                  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
  x = (x ^ (x >>  1)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
  x = (x ^ (x >>  2)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
  x = (x ^ (x >>  4)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
  x = (x ^ (x >>  8)) & 0x0000ffff; // x = ---- ---- ---- ---- fedc ba98 7654 3210
  return x;
}

uint32_t EncodeMorton2(uint32_t x, uint32_t y) {
        return (Part1By1(y) << 1) + Part1By1(x); }
uint32_t DecodeMorton2X(uint32_t code) { return Compact1By1(code >> 0); }
uint32_t DecodeMorton2Y(uint32_t code) { return Compact1By1(code >> 1); }

/*
	Sort in Morton order
	geom should be contained in the unit sq
*/
std::vector<int> zsort(const Eigen::MatrixXd &geom)
{
        int n = geom.rows();
        vector<uint32_t> z(n);
        for (int i = 0; i < n; i++) {
                uint16_t x = static_cast<uint16_t>  (geom(i, 0) * (double)
                        UINT16_MAX + .5);  //  (1 << 16) ); 
                uint16_t y = static_cast<uint16_t>  (geom(i, 1) * (double)
                        UINT16_MAX + .5);  // (1 << 16) ); 
                z[i] = EncodeMorton2(x, y);
        }
        std::vector<int> idx(n);
        std::iota(idx.begin() , idx.end() , 0);
        std::sort(idx.begin() , idx.end() , [&z](int i , int j)
                {return z[i] < z[j];});
        return idx;
}

// [[Rcpp::export]]
Eigen::VectorXi zorder(const Eigen::MatrixXd &geom)
{
	if(geom.rows() == 0) Rcpp::stop("The number of rows of geom cannot be"
		" zero\n");
        if(geom.cols() != 2) Rcpp::stop("The geometry for `mvn.tlr` should be 2D, "
                "each row represents one location\n");
        if(geom.maxCoeff() > 1.0 || geom.minCoeff() < 0.0) Rcpp::stop("The geometry "
                "for `mvn.tlr` should be contained in the unit square\n");

        int n = geom.rows();
        vector<uint32_t> z(n);
        for (int i = 0; i < n; i++) {
                uint16_t x = static_cast<uint16_t>  (geom(i, 0) * (double)
                        UINT16_MAX + .5);  //  (1 << 16) ); 
                uint16_t y = static_cast<uint16_t>  (geom(i, 1) * (double)
                        UINT16_MAX + .5);  // (1 << 16) ); 
                z[i] = EncodeMorton2(x, y);
        }
        VectorXi idx(n);
        std::iota(idx.data(), idx.data()+n, 0);
        std::sort(idx.data(), idx.data()+n , [&z](int i , int j)
                {return z[i] < z[j];});
	std::for_each(idx.data(), idx.data()+n, [](int &x){x++;});
        return idx;
}

