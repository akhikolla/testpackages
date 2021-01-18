#pragma once
#ifndef SOBOL_POINT_H
#define SOBOL_POINT_H

#include <iostream>
#include <inttypes.h>
#if defined(IN_RCPP)
#include <Rcpp.h>
#endif

// [[Rcpp::plugins(cpp11)]]

namespace {
    inline int get_sobol_s_max() {
        return 21201;
    }
    inline int get_sobol_s_min() {
        return 2;
    }
    inline int get_sobol_m_max() {
        return 31;
    }
    inline int get_sobol_m_min() {
        return 8;
    }
}

namespace DigitalNetNS {
#if defined(USE_FILE)
    bool get_sobol_base(std::istream& is,
                        uint32_t s, uint32_t m,  uint64_t base[]);
#endif
#if defined(USE_DF)
    bool read_sobol_base(Rcpp::DataFrame df,
                         uint32_t s, uint32_t m,  uint64_t base[]);
#endif
#if defined(USE_SQL)
    bool select_sobol_base(const std::string& path,
                           uint32_t s, uint32_t m,  uint64_t base[]);
#endif
}
#endif //SOBOL_POINT_H
