#include "core/utils/random.hpp"
#include "core/exceptions/OperationNotSupportedException.hpp"
#include <algorithm>
#include <iostream>

namespace uu {
namespace core {

std::mt19937 &
get_random_engine(
)
{
    static std::mt19937 engine;
    static bool seed = true;

    if (seed)
    {
        engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        seed = false;
    }

    return engine;
}

size_t
getRandomInt(
    size_t max
)
{
    std::uniform_int_distribution<int> distribution(0,max-1);
    return distribution(get_random_engine());

}

long
getRandomLong(
    long max
)
{
    std::uniform_int_distribution<long> distribution(0,max-1);
    return distribution(get_random_engine());

}

size_t
get_binomial(
    size_t tests,
    double p
)
{
    std::binomial_distribution<size_t> distribution(tests, p);
    return distribution(get_random_engine());

}



double
drand()
{
    std::uniform_real_distribution<double> distribution(0,1);
    return distribution(get_random_engine()); // [0,1[
}

size_t
random_level(
    size_t MAX_LEVEL,
    double P
)
{
    double r = drand();

    if (r==0)
    {
        r=1;    // avoid taking logarithm of 0
    }

    double num = std::log(r);
    double denum = std::log(1.0-P);
    size_t lvl = (size_t)(num/denum);
    return lvl < MAX_LEVEL ? lvl : MAX_LEVEL;
}


std::set<size_t>
getKRandom(
    size_t max,
    size_t k
)
{
    if (max<k)
    {
        throw OperationNotSupportedException("Only " + std::to_string(max) + " values available, requested " + std::to_string(k));
    }

    std::set<size_t> res;

    while (res.size()<k)
    {
        res.insert(getRandomInt(max));
    }

    return res;
}


std::vector<size_t>
get_k_uniform(
    size_t max,
    size_t k
)
{

    std::vector<size_t> res(k, 0);

    size_t rand = getRandomInt(max);

    size_t last_position = 1;

    res[0] = rand;

    for (size_t i = 1; i < k; i++)
    {
        rand = getRandomInt(max-i);

        size_t pos = 0;

        while (pos < last_position && res[pos] <= rand)
        {
            rand++;
            pos++;
        }

        last_position++;

        for (size_t idx = last_position - 1; idx > pos; idx--)
        {
            res[idx] = res[idx-1];
        }

        res[pos] = rand;

    }

    return res;
}

bool
test(
    double probability
)
{
    std::bernoulli_distribution distribution(probability);
    return distribution(get_random_engine());
}

size_t
test(
    const std::vector<double>& options
)
{
    // For efficiency reasons, we do not check if the values sum to 1
    double prob_failing_previous_tests=1;

    for (size_t idx=0; idx<options.size()-1; idx++)
    {
        double adjusted_prob = options.at(idx)/prob_failing_previous_tests;

        if (test(adjusted_prob))
        {
            return idx;
        }

        prob_failing_previous_tests *= (1-adjusted_prob);
    }

    // In practice, the last value of the input is assumed to be 1 minus the sum of the previous values
    return options.size()-1;
}

size_t
test(
    const std::vector<std::vector<double> >& options,
    size_t row_num
)
{
    // For efficiency reasons, we do not check if the values sum to 1
    double prob_failing_previous_tests=1;

    for (size_t idx=0; idx<options.at(row_num).size()-1; idx++)
    {
        double adjusted_prob = options.at(row_num).at(idx)/prob_failing_previous_tests;

        if (test(adjusted_prob))
        {
            return idx;
        }

        prob_failing_previous_tests *= (1-adjusted_prob);
    }

    // In practice, the last value of the input is assumed to be 1 minus the sum of the previous values
    return options.at(row_num).size()-1;
}

} // namespace core
} // namespace uu
