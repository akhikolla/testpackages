/**
 * This header defines:
 * basic mathematical/statistical functions (mean, standard deviation, intersection, union).
 *
 * History:
 * - 2018.01.01 file imported from version 1.0 of the multinet library
 */

#ifndef UU_CORE_UTILS_MATH_H_
#define UU_CORE_UTILS_MATH_H_

#include <set>
#include <unordered_set>
#include <cmath>
#include <vector>
#include <type_traits>
#include <numeric>

namespace uu {
namespace core {

/**
 * Mathematical mean (sum of elements divided by number of elements).
 * @param vec input set of values
 * @return the mean of the input values
 */
template <typename InputIterator>
double
mean(
    InputIterator begin,
    InputIterator end
)
{
    double sum = 0;
    int n = 0;

    for (auto it=begin; it!=end; ++it)
    {
        sum += *it;
        n++;
    }

    return sum / n;
}

/**
 * Standard deviation of a population.
 * @param vec input set of values
 * @return the standard deviation of the input values
 */
template <typename InputIterator>
double
stdev(
    InputIterator begin,
    InputIterator end
)
{

    // mean
    double m = mean(begin,end);

    // variance
    double variance = 0.0;

    int n = 0;

    for (auto it=begin; it!=end; ++it)
    {
        variance += (*it - m) * (*it - m);
        n++;
    }

    // The POPULATION stdev is computed: divide by n
    variance /= n;

    // standard deviation
    return std::sqrt(variance);
}

/**
 * Tha Jaccard similarity of a set of sets is the size of their intersection divided by the size of their union.
 * @param sets a vector of sets
 * @return Jaccard similiarity of the input sets. If the sets are all empty, 0 is returned.
 */
template <class T>
double
jaccard_similarity(
    const std::vector<std::unordered_set<T> >& sets
)
{
    long union_size = s_union(sets).size();

    if (union_size==0)
    {
        return 0;
    }

    long intersection_size = s_intersection(sets).size();
    return (double)intersection_size/union_size;
}

/**
 * Set-based intersection, for a combination of sorted and unordered sets.
 * @param set1 a set of values
 * @param set2 a set of values
 * @return the intersection of the two input sets
 */
template <class T>
int
intersection_size(
    const std::unordered_set<T>& set1,
    const std::unordered_set<T>& set2
)
{
    int common_elements = 0;

    if (set1.size()<set2.size())
    {
        for (T el: set1)
        {
            if (set2.count(el)>0)
            {
                common_elements++;
            }
        }
    }

    else
    {
        for (T el: set2)
        {
            if (set1.count(el)>0)
            {
                common_elements++;
            }
        }
    }

    return common_elements;
}

/**
 * Set-based intersection, for unordered sets.
 * @param sets a vector of sets
 * @return the intersection of the input sets
 */
template <class T>
std::unordered_set<T>
s_intersection(
    const std::vector<std::unordered_set<T> >& sets
)
{
    std::unordered_set<T> result;
    size_t idx = 0; // index of the smallest set

    for (size_t i=1; i<sets.size(); i++)
    {
        if (sets.at(i).size() < sets.at(idx).size())
        {
            idx=i;
        }
    }

    for (T element: sets.at(idx))
    {
        bool in_intersection = true;

        for (size_t i=0; i<sets.size(); i++)
        {
            if (i==idx)
            {
                continue;
            }

            if (sets.at(i).count(element)==0)
            {
                in_intersection = false;
                break;
            }
        }

        if (in_intersection)
        {
            result.insert(element);
        }
    }

    return result;
}

/**
 * Set-based intersection, for sorted sets.
 * @param sets a vector of sets
 * @return the intersection of the input sets
 */
template <class T>
std::unordered_set<T>
s_intersection(
    const std::vector<std::set<T> >& sets
)
{
    // NOTE: it can be made more efficient exploiting sorting
    std::unordered_set<T> result;
    size_t idx = 0; // index of the smallest set

    for (size_t i=1; i<sets.size(); i++)
    {
        if (sets.at(i).size() < sets.at(idx).size())
        {
            idx=i;
        }
    }

    for (T element: sets.at(idx))
    {
        bool in_intersection = true;

        for (size_t i=0; i<sets.size(); i++)
        {
            if (i==idx)
            {
                continue;
            }

            if (sets.at(i).count(element)==0)
            {
                in_intersection = false;
                break;
            }
        }

        if (in_intersection)
        {
            result.insert(element);
        }
    }

    return result;
}


/**
 * Set-based intersection, for unordered sets.
 * @param set1 a set of values
 * @param set2 a set of values
 * @return the intersection of the two input sets
 */
template <class T>
std::unordered_set<T>
s_intersection(
    const std::unordered_set<T>& set1,
    const std::unordered_set<T>& set2
)
{
    std::vector<std::unordered_set<T> > sets({set1,set2});
    return s_intersection(sets);
}

/**
 * Set-based intersection, for sorted sets.
 * @param set1 a set of values
 * @param set2 a set of values
 * @return the intersection of the two input sets
 */
template <class T>
std::unordered_set<T>
s_intersection(
    const std::set<T>& set1,
    const std::set<T>& set2
)
{
    std::vector<std::set<T> > sets({set1,set2});
    return s_intersection(sets);
}


/**
 * Set-based intersection, for a combination of sorted and unordered sets.
 * @param set1 a set of values
 * @param set2 a set of values
 * @return the intersection of the two input sets
 */
template <class T>
std::unordered_set<T>
s_intersection(
    const std::set<T>& set1,
    const std::unordered_set<T>& set2
)
{
    std::unordered_set<T> result;

    for (T element: set1)
    {
        if (set2.count(element)>0)
        {
            result.insert(element);
        }
    }

    return result;
}

/**
 * Set-based intersection, for a combination of sorted and unordered sets.
 * @param set1 a set of values
 * @param set2 a set of values
 * @return the intersection of the two input sets
 */
template <class T>
std::unordered_set<T>
s_intersection(
    const std::unordered_set<T>& set1,
    const std::set<T>& set2
)
{
    std::unordered_set<T> result;

    for (T element: set2)
    {
        if (set1.count(element)>0)
        {
            result.insert(element);
        }
    }

    return result;
}


/**
 * Set-based union, for unordered sets.
 * @param sets a vector of sets
 * @return the union of the input sets
 */
template <class T>
std::unordered_set<T>
s_union(
    const std::vector<std::unordered_set<T> >& sets
)
{
    std::unordered_set<T> result;

    for (std::unordered_set<T> S: sets)
    {
        result.insert(S.begin(), S.end());
    }

    return result;
}

/**
 * Set-based union, for sorted sets.
 * @param sets a vector of sets
 * @return the union of the input sets
 */
template <class T>
std::unordered_set<T>
s_union(
    const std::vector<std::set<T> >& sets
)
{
    std::unordered_set<T> result;

    for (std::set<T> S: sets)
    {
        result.insert(S.begin(), S.end());
    }

    return result;
}

/**
 * Set-based union, for two sets.
 * @param set1 a set of values
 * @param set2 a set of values
 * @return the union of the two input sets
 */
template <class T>
std::unordered_set<T>
s_union(
    const std::unordered_set<T>& set1,
    const std::unordered_set<T>& set2
)
{
    std::unordered_set<T> result;
    result.insert(set1.begin(), set1.end());
    result.insert(set2.begin(), set2.end());
    return result;
}

/**
 * Set-based union, for two sorted sets.
 * @param set1 a set of values
 * @param set2 a set of values
 * @return the union of the two input sets
 */
template <class T>
std::unordered_set<T>
s_union(
    const std::set<T>& set1,
    const std::set<T>& set2
)
{
    std::unordered_set<T> result;
    result.insert(set1.begin(), set1.end());
    result.insert(set2.begin(), set2.end());
    return result;
}


}
}

#endif
