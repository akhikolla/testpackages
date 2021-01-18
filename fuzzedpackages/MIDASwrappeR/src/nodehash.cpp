#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

#include <iostream>
#include <limits>
#include "nodehash.h"
#include <Rcpp.h>

using namespace Rcpp;

Nodehash::Nodehash(int r, int b)
{
    num_rows = r;
    num_buckets = b;
    hash_a.resize(num_rows);
    hash_b.resize(num_rows);
    for (int i = 0; i < r; i++) {
        // a is in [1, p-1]; b is in [0, p-1]
        hash_a[i] = floor(unif_rand() * (num_buckets - 1) + 1);
        hash_b[i] = floor(unif_rand() * (num_buckets) + 1);
    }
    this->clear();
}

Nodehash::~Nodehash()
{
}

int Nodehash::hash(int a, int i)
{
    int resid = (a * hash_a[i] + hash_b[i]) % num_buckets;
    return resid + (resid < 0 ? num_buckets : 0);
}

void Nodehash::insert(int a, double weight)
{
    int bucket;
    for (int i = 0; i < num_rows; i++) {
        bucket = hash(a, i);
        count[i][bucket] += weight;
    }
}

double Nodehash::get_count(int a)
{
    double min_count = numeric_limits<double>::max();
    int bucket;
    for (int i = 0; i < num_rows; i++) {
        bucket = hash(a, i);
        min_count = MIN(min_count, count[i][bucket]);
    }
    return min_count;
}

void Nodehash::clear()
{
    count = vector<vector<double> >(num_rows, vector<double>(num_buckets, 0.0));
}

void Nodehash::lower(double factor)
{
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_buckets; j++) {
            count[i][j] = count[i][j] * factor;
        }
    }
}
