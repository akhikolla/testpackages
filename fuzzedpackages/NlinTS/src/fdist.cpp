/**
 * @authors Hmamouche Youssef
 **/

#include "cmath"
#include <tgmath.h>

#include "../inst/include/fdist.h"

#define PI 3.14159265358

// Beta function
double beta(double a, double b) {
    if (a <= 0 || b <= 0) {
        return 0;
    }
    return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
}

double fi(int N, double x, double a, double b) {
    int n = N / 2;
    double f = 0.0, f1, s1, s2, tU, tV;
    int i;
    for (i = n; i >= 1; i--) {
        tU = (a + 2.0 * i - 1.0) * (a + 2.0 * i);
        s2 = i * (b - i) * x / tU;
        f1 = s2 / (1.0 + f);
        tV = (a + 2.0 * i - 2.0) * (a + 2.0 * i - 1.0);
        s1 = -(a + i - 1.0) * (b + a + i - 1.0) * x / tV;
        f = s1 / (1.0 + f1);
    }
    return 1.0 / (1.0 + f);
}

double inBeta(double x, double a, double b) {
    if (a <= 0.0 || b <= 0.0) {
        return 0.0;
    }
    if (fabs(x - 0.0) < 1.0e-15 || fabs(x - 1.0) < 1.0e-15) {
        return 0.0;
    }

    double c1, c2, c3, f1, f2;
    int n;
    c1 = pow(x, a);
    c2 = pow(1.0 - x, b);
    c3 = beta(a, b);
    if (x < (a + 1.0) / (a + b + 2.0)) {
        n = 1;
        while (1) {
            f1 = fi(2 * n, x, a, b);
            f2 = fi(2 * n + 2, x, a, b);
            if (fabs(f2 - f1) < 1.0e-15)
                return f2 * c1 * c2 / a / c3;
            else
                n++;
        }
    } else {
        if (fabs(x - 0.5) < 1.0e-15 && fabs(a - b) < 1.0e-15)
            return 0.5;
        else {
            n = 1;
            while (1) {
                f1 = fi(2 * n, 1.0 - x, b, a);
                f2 = fi(2 * n + 2, 1.0 - x, b, a);
                if (fabs(f2 - f1) < 1.0e-15)
                    return 1.0 - f2 * c1 * c2 / b / c3;
                else
                    n++;
            }
        }
    }
    return 0;
}
 // p-value of the f-test
double getPvalue(double f, double n1, double n2) {
    if (f <= 0.0)
        return 1.0;
    else
      return inBeta(n2 / (n2 + n1 * f), n2 / 2.0, n1 / 2.0);
}

// Student distribution
double getStudent (double t, double n)
{
    if ( t < 0)
        t = -t;
    double a = 1 / sqrt (PI * n);
    double b = tgamma ((n + 1) / 2) / tgamma (n / 2);
    double c = pow ((1 +  (t * t) / n), - (n + 1) / 2);
    return a * b * c;
}
