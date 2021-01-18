#include <testthat.h>
#include "intrapolate.h"

using namespace PSQN;

context("testing intrapolate class") {
  test_that("intrapolate gives the correct result with a 2nd order poly") {
    constexpr double const f0 = 0,
                           d0 = -1,
                           x0 = 0,
                            f = 3.75,
                            x = 2.5;

    intrapolate inter(f0, d0, x, f);
    double const val = inter.get_value(-2, 3);

    expect_true(std::abs((val - .5) / .5) < 1e-8);
  }

  test_that("intrapolate gives the correct result with a 3rd order poly") {
    constexpr double const f0 = 0,
                           d0 = -1,
                           x0 = 0,
                         fold = 5.3125,
                         xold = 2.5,
                         fnew = -0.2336,
                         xnew = .4,
                        truth = 0.467251416997127;

    intrapolate inter(f0, d0, xold, fold);
    inter.update(xnew, fnew);
    double const val = inter.get_value(xnew, xold);

    expect_true(std::abs((val - truth) / truth) < 1e-8);
  }
}
