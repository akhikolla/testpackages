#ifndef GRAYINDEX_H
#define GRAYINDEX_H
/**
 * @file grayindex.h
 *
 * @brief gray code class
 *
 * @author Shinsuke Mori (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 * @author Mutsuo Saito
 *
 * Copyright (C) 2017 Shinsuke Mori, Makoto Matsumoto, Mutsuo Saito
 * and Hiroshima University.
 * All rights reserved.
 *
 * The GPL ver.3 is applied to this software, see
 * COPYING
 */
#include "bit_operator.h"
#include <stdint.h>

// [[Rcpp::plugins(cpp11)]]

class GrayIndex {
public:
    GrayIndex() {
        count = 1;
    }
    void clear() {
        count = 1;
    }
    void next() {
        count++;
    }
    int index() {
        return tailingZeroBit(count);
    }
private:
    uint64_t count;
};
#endif // GRAYINDEX_HPP
