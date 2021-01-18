/**
 * Functions to summarize a property matrix.
 *
 * Several generic summarization functions can be computed on a property matrix, e.g., to
 * obtain the amount of overlapping or the statistical correlation between different types
 * of structures in different contexts.
 *
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_DATASTRUCTURES_PROPERTYMATRIX_SUMMARIZATION_H_
#define UU_CORE_DATASTRUCTURES_PROPERTYMATRIX_SUMMARIZATION_H_

#include <algorithm>
#include <cmath>
#include <limits>
#include "core/utils/math.hpp"
#include "core/utils/Counter.hpp"
#include "core/datastructures/propertymatrix/PropertyMatrix.hpp"
#include "core/datastructures/propertymatrix/StructureComparisonFunction.hpp"


namespace uu {
namespace core {

/**
 * The result of the comparison of two binary vectors, comparing their elements at each coordinate
 * and counting the number of occurrences for all possible configurations
 * (true-true, true-false, false-true, false-false)
 */
struct binary_vector_comparison
{
  public:
    /** number of coordinates where both vectors are true */
    long yy = 0;
    /** number of coordinates where the first vector is true and the second is false */
    long yn = 0;
    /** number of coordinates where the first vector is false and the second is true */
    long ny = 0;
    /** number of coordinates where both vectors are false */
    long nn = 0;
};

/**
 * Basic context summaries
 */

template <class STRUCTURE, class CONTEXT, class VALUE>
double
min(
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c
)
{
    double min = std::numeric_limits<double>::infinity();

    long checked_columns = 0;

    for (STRUCTURE s: P.structures())
    {
        Value<VALUE> v = P.get(s,c);

        if (!v.null && min>v.value)
        {
            min = v.value;
        }

        checked_columns++;
    }

    if ((P.num_structures>checked_columns) && min>P.get_default())
    {
        min = P.get_default();
    }

    return min;
}


template <class STRUCTURE, class CONTEXT, class VALUE>
double
max(
    const  PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c
)
{
    double max = -std::numeric_limits<double>::infinity();

    long checked_columns = 0;

    for (STRUCTURE s: P.structures())
    {
        Value<VALUE> v = P.get(s,c);

        if (!v.null && max<v.value)
        {
            max = v.value;
        }

        checked_columns++;
    }

    if ((P.num_structures>checked_columns) && max<P.get_default())
    {
        max = P.get_default();
    }

    return max;
}

template <class STRUCTURE, class CONTEXT, class VALUE>
double
sum(
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c
)
{
    double sum = 0.0;

    long checked_columns = 0;

    for (STRUCTURE s: P.structures())
    {
        Value<VALUE> v = P.get(s,c);

        if (!v.null)
        {
            sum += (double)v.value;
        }

        checked_columns++;
    }

    sum += P.get_default()*(P.num_structures-checked_columns);
    return sum;
}

template <class STRUCTURE, class CONTEXT, class VALUE>
double
mean(
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c
)
{
    double _mean = 0.0;

    long checked_columns = 0;

    for (STRUCTURE s: P.structures())
    {
        Value<VALUE> v = P.get(s,c);

        if (!v.null)
        {
            _mean += (double)v.value;
        }

        checked_columns++;
    }

    _mean += P.get_default()*(P.num_structures-checked_columns);
    return _mean/(P.num_structures-P.num_na(c));
}

template <class STRUCTURE, class CONTEXT, class VALUE>
double
sd(
    const  PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c
)
{
    double _mean = mean(P,c);

    double sd = 0.0;

    long checked_columns = 0;

    for (STRUCTURE s: P.structures())
    {
        Value<VALUE> v = P.get(s,c);

        if (!v.null)
        {
            sd += (double)pow(v.value-_mean,2);
        }

        checked_columns++;
    }

    sd += (double)std::pow(P.get_default()-_mean,2)*(P.num_structures-checked_columns);
    return std::sqrt(sd/(P.num_structures-P.num_na(c)));
}

template <class STRUCTURE, class CONTEXT, class VALUE>
double
skew(
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c
)
{
    double _mean = mean(P,c);

    double sd = 0.0;
    double skew = 0.0;

    long checked_columns = 0;

    for (STRUCTURE s: P.structures())
    {
        Value<VALUE> v = P.get(s,c);

        if (!v.null)
        {
            sd += (double)pow(v.value-_mean,2);
            skew += (double)pow(v.value-_mean,3);
        }

        checked_columns++;
    }

    sd += (double)std::pow(P.get_default()-_mean,2)*(P.num_structures-checked_columns);
    sd = std::sqrt(sd/(P.num_structures-P.num_na(c)));

    skew += (double)std::pow(P.get_default()-_mean,3)*(P.num_structures-checked_columns);
    return skew/std::pow(sd,3)/(P.num_structures-P.num_na(c));
}


template <class STRUCTURE, class CONTEXT, class VALUE>
double
kurt(
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c
)
{
    double _mean = mean(P,c);

    double sd = 0.0;
    double kurt = 0.0;

    long checked_columns = 0;

    for (STRUCTURE s: P.structures())
    {
        Value<VALUE> v = P.get(s,c);

        if (!v.null)
        {
            sd += (double)pow(v.value-_mean,2);
            kurt += (double)pow(v.value-_mean,4);
        }

        checked_columns++;
    }

    sd += (double)std::pow(P.get_default()-_mean,2)*(P.num_structures-checked_columns);
    sd = std::sqrt(sd/(P.num_structures-P.num_na(c)));

    kurt += (double)std::pow(P.get_default()-_mean,4)*(P.num_structures-checked_columns);

    return kurt/std::pow(sd,4)/(P.num_structures-P.num_na(c));
}

template <class STRUCTURE, class CONTEXT, class VALUE>
double
entropy(
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c
)
{
    Counter<double> count;

    double entropy = 0.0;

    long checked_columns = 0;


    for (STRUCTURE s: P.structures())
    {
        Value<VALUE> v = P.get(s,c);

        if (!v.null)
        {
            count.inc(v.value);
        }

        checked_columns++;
    }

    count.set(P.get_default(),count.count(P.get_default())+(P.num_structures-checked_columns));

    for (auto pair: count.map())
    {
        double fr = (double)pair.second/(P.num_structures-P.num_na(c));

        if (fr!=0)
        {
            entropy += -fr*std::log(fr);
        }
    }

    return entropy;
}

template <class STRUCTURE, class CONTEXT, class VALUE>
double
CV(
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c
)
{
    double _mean=mean(P,c);
    return sd(P,c)/_mean;
}

template <class STRUCTURE, class CONTEXT, class VALUE>
double
jarque_bera(
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c
)
{
    return (P.num_structures-P.num_na(c))/6.0*(std::pow(skew(P,c),2)+std::pow(kurt(P,c)-3,2)/4);
}

template <class STRUCTURE, class CONTEXT, class VALUE>
double
relative_difference(
    double (*f)(const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
                const CONTEXT& c),
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c1,
    const CONTEXT& c2
)
{
    double f1 = (*f)(P,c1);
    double f2 = (*f)(P,c2);

    if (std::abs(f1)+std::abs(f2)==0)
    {
        return 0;
    }

    return std::abs(f1-f2)*2/(std::abs(f1)+std::abs(f2));
}


/*
 * K number of bins -
 * assumes only non-negative -
 * check if empty bins
 */
template <class STRUCTURE, class CONTEXT, class VALUE>
std::pair<Counter<int>,Counter<int> >
histograms(
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c1,
    const CONTEXT& c2,
    int K
)
{
    Counter<int> hist1;
    Counter<int> hist2;

    // build histograms
    double _min = std::min(min(P,c1),min(P,c2));
    double _max = std::max(max(P,c1),max(P,c2));
    long checked_columns = 0;

    for (STRUCTURE s: P.structures())
    {
        Value<VALUE> val1 = P.get(s,c1);

        if (!val1.null)
        {
            int v1 = std::floor(((double)val1.value-_min)*K/(_max-_min));

            if (v1==K)
            {
                v1=K-1;
            }

            hist1.inc(v1);
        }

        Value<VALUE> val2 = P.get(s,c2);

        if (!val2.null)
        {
            int v2 = std::floor(((double)val2.value-_min)*K/(_max-_min));

            if (v2==K)
            {
                v2=K-1;
            }

            hist2.inc(v2);
        }

        checked_columns++;
    }

    hist1.set(P.get_default(),hist1.count(P.get_default())+(P.num_structures-checked_columns));
    hist2.set(P.get_default(),hist2.count(P.get_default())+(P.num_structures-checked_columns));

    return std::pair<Counter<int>,Counter<int> >(hist1,hist2);
}


template <class STRUCTURE, class CONTEXT, class VALUE>
double
dissimilarity_index(
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c1,
    const CONTEXT& c2,
    int K
)
{

    std::pair<Counter<int>,Counter<int> > h = histograms(P,c1,c2,K);

    // compare histograms
    double diss = 0;

    for (int i=0; i<K; i++)
    {
        double fr1 = (double)h.first.count(i)/(P.num_structures-P.num_na(c1));
        double fr2 = (double)h.second.count(i)/(P.num_structures-P.num_na(c2));
        diss += std::abs(fr1-fr2)*.5;
    }

    return diss;
}

// epsilon correction to avoid division by 0
template <class STRUCTURE, class CONTEXT, class VALUE>
double
KL_divergence(
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c1,
    const CONTEXT& c2,
    int K
)
{

    std::pair<Counter<int>,Counter<int> > h = histograms(P,c1,c2,K);

    // compare histograms
    double diss = 0;

    for (int i=0; i<K; i++)
    {
        int num_elements_h1 = K+P.num_structures-P.num_na(c1);
        int num_elements_h2 = K+P.num_structures-P.num_na(c2);
        double fr1 = ((double)h.first.count(i)+1)/num_elements_h1;
        double fr2 = ((double)h.second.count(i)+1)/num_elements_h2;

        if (fr1!=0)
        {
            diss += fr1*std::log(fr1/fr2);
        }
    }

    return diss;
}


// epsilon correction to avoid division by 0
template <class STRUCTURE, class CONTEXT, class VALUE>
double
JS_divergence(
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c1,
    const CONTEXT& c2,
    int K
)
{

    std::pair<Counter<int>,Counter<int> > h = histograms(P,c1,c2,K);

    // compare histograms
    double diss = 0;

    for (int i=0; i<K; i++)
    {
        int num_elements_h1 = K+P.num_structures-P.num_na(c1);
        int num_elements_h2 = K+P.num_structures-P.num_na(c2);

        double fr1 = ((double)h.first.count(i)+1)/num_elements_h1;
        double fr2 = ((double)h.second.count(i)+1)/num_elements_h2;

        double fr_joint = .5*fr1+.5*fr2;

        if (fr_joint!=0)
        {
            diss += .5*fr1*std::log(fr1/fr_joint) + .5*fr2*std::log(fr2/fr_joint);
        }
    }

    return diss;
}

template <class STRUCTURE, class CONTEXT, class VALUE>
double
jeffrey_divergence(
    const PropertyMatrix<STRUCTURE,CONTEXT,VALUE>& P,
    const CONTEXT& c1,
    const CONTEXT& c2,
    int K
)
{

    std::pair<Counter<int>,Counter<int> > h = histograms(P,c1,c2,K);

    // compare histograms
    double diss = 0;

    for (int i=0; i<K; i++)
    {
        double fr1 = (double)h.first.count(i)/(P.num_structures-P.num_na(c1));
        double fr2 = (double)h.second.count(i)/(P.num_structures-P.num_na(c2));

        if (fr1!=0 && fr2!=0)
        {
            diss += fr1*std::log(fr1/fr2) + fr2*std::log(fr2/fr1);
        }
    }

    return diss;
}


/**
 * Compares two binary vectors, comparing their elements at each coordinate
 * and counting the number of occurrences for all possible configurations (true-true, true-false, false-true, false-false)
 */
template <class STRUCTURE, class CONTEXT>
binary_vector_comparison
compare(
    const PropertyMatrix<STRUCTURE,CONTEXT,bool>& P,
    const CONTEXT& c1,
    const CONTEXT& c2
)
{
    binary_vector_comparison res;
    res.yy = 0;
    res.yn = 0;
    res.ny = 0;
    res.nn = 0;

    long checked_columns = 0;

    for (STRUCTURE s: P.structures())
    {
        Value<bool> v1 = P.get(s,c1);
        bool p1 = v1.value && !v1.null;
        Value<bool> v2 = P.get(s,c2);
        bool p2 = v2.value && !v2.null;

        if (p1 && p2)
        {
            res.yy++;
        }

        else if (p1)
        {
            res.yn++;
        }

        else if (p2)
        {
            res.ny++;
        }

        else
        {
            res.nn++;
        }

        checked_columns++;
    }

    if (P.get_default())
    {
        res.yy += P.num_structures - checked_columns;
    }

    else
    {
        res.nn += P.num_structures - checked_columns;
    }

    return res;
}

template <class STRUCTURE, class CONTEXT>
double
russell_rao(
    const PropertyMatrix<STRUCTURE,CONTEXT,bool>& P,
    const CONTEXT& c1,
    const CONTEXT& c2
)
{
    binary_vector_comparison comp = compare(P, c1, c2);
    return (double)(comp.yy)/(comp.yy+comp.ny+comp.yn+comp.nn);
}

template <class STRUCTURE, class CONTEXT>
double
jaccard(
    const PropertyMatrix<STRUCTURE,CONTEXT,bool>& P,
    const CONTEXT& c1,
    const CONTEXT& c2
)
{
    binary_vector_comparison comp = compare(P, c1, c2);
    return (double)comp.yy/(comp.yy+comp.yn+comp.ny);
}

template <class STRUCTURE, class CONTEXT>
double
coverage(
    const PropertyMatrix<STRUCTURE,CONTEXT,bool>& P,
    const CONTEXT& c1,
    const CONTEXT& c2
)
{
    binary_vector_comparison comp = compare(P, c1, c2);
    return (double)comp.yy/(comp.yy+comp.ny);
}

template <class STRUCTURE, class CONTEXT>
double
kulczynski2(
    const PropertyMatrix<STRUCTURE,CONTEXT,bool>& P,
    const CONTEXT& c1,
    const CONTEXT& c2
)
{
    binary_vector_comparison comp = compare(P, c1, c2);
    return ((double)comp.yy/(comp.yy+comp.yn)+(double)comp.yy/(comp.yy+comp.ny))/2;
}

template <class STRUCTURE, class CONTEXT>
double
simple_matching(
    const PropertyMatrix<STRUCTURE,CONTEXT,bool>& P,
    const CONTEXT& c1,
    const CONTEXT& c2
)
{
    binary_vector_comparison comp = compare(P, c1, c2);
    return (double)(comp.yy+comp.nn)/(comp.yy+comp.ny+comp.yn+comp.nn);
}

template <class STRUCTURE, class CONTEXT>
double
hamann(
    const PropertyMatrix<STRUCTURE,CONTEXT,bool>& P,
    const CONTEXT& c1,
    const CONTEXT& c2
)
{
    binary_vector_comparison comp = compare(P, c1, c2);
    return (double)(comp.yy+comp.nn-comp.yn-comp.ny)/(comp.yy+comp.ny+comp.yn+comp.nn);
}


template <class STRUCTURE, class CONTEXT>
double
L2(
    const PropertyMatrix<STRUCTURE,CONTEXT,double>& P,
    const CONTEXT& c1,
    const CONTEXT& c2
)
{
    double dist = 0;

    for (STRUCTURE s: P.structures())
    {
        Value<double> v1 = P.get(s,c1);
        Value<double> v2 = P.get(s,c2);

        if (!v1.null && !v2.null)
        {
            dist += (v1.value-v2.value)*(v1.value-v2.value);
        }
    }

    return std::sqrt(dist);
}

template <class STRUCTURE, class CONTEXT>
double
cosine(
    const PropertyMatrix<STRUCTURE,CONTEXT,double>& P,
    const CONTEXT& c1,
    const CONTEXT& c2
)
{
    double product = 0;
    double norm1 = 0;
    double norm2 = 0;

    long checked_columns = 0;
    double default_val = P.get_default();

    for (STRUCTURE s: P.structures())
    {
        Value<double> v1 = P.get(s,c1);
        Value<double> v2 = P.get(s,c2);

        if (!v1.null && !v2.null)
        {
            product += v1.value*v2.value;
            norm1 += v1.value*v1.value;
            norm2 += v2.value*v2.value;
        }

        checked_columns++;
    }

    product += default_val*default_val*(P.num_structures-checked_columns);
    norm1 += default_val*default_val*(P.num_structures-checked_columns);
    norm2 += default_val*default_val*(P.num_structures-checked_columns);

    return product/std::sqrt(norm1)/std::sqrt(norm2);
}

template <class STRUCTURE, class CONTEXT>
double
pearson(
    const PropertyMatrix<STRUCTURE,CONTEXT,double>& P,
    const CONTEXT& c1,
    const CONTEXT& c2
)
{
    double cov = 0;
    double mean1 = 0;
    double mean2 = 0;
    double std1 = 0;
    double std2 = 0;

    long checked_columns = 0;
    long num_incomplete = 0;
    double default_val = P.get_default();

    for (STRUCTURE s: P.structures())
    {
        Value<double> v1 = P.get(s,c1);
        Value<double> v2 = P.get(s,c2);

        if (!v1.null && !v2.null)
        {
            mean1 += v1.value;
            mean2 += v2.value;
        }

        else
        {
            num_incomplete++;
        }

        checked_columns++;
    }

    mean1 += default_val*(P.num_structures-checked_columns);
    mean2 += default_val*(P.num_structures-checked_columns);
    mean1 /= (P.num_structures-num_incomplete);
    mean2 /= (P.num_structures-num_incomplete);

    checked_columns = 0;

    for (STRUCTURE s: P.structures())
    {
        Value<double> v1 = P.get(s,c1);
        Value<double> v2 = P.get(s,c2);

        if (!v1.null && !v2.null)
        {
            cov += (v1.value-mean1)*(v2.value-mean2);
            std1 += (v1.value-mean1)*(v1.value-mean1);
            std2 += (v2.value-mean2)*(v2.value-mean2);
        }

        //if (!v1.null) std1 += (v1.value-_mean1)*(v1.value-_mean1);
        //if (!v2.null) std2 += (v2.value-_mean2)*(v2.value-_mean2);
        checked_columns++;
        //std::cout << val2 << " " << mean2 << " " << ((val2-mean2)*(val2-mean2)) << std::endl;
    }

    cov += (default_val-mean1)*(default_val-mean2)*(P.num_structures-checked_columns);
    std1 += (default_val-mean1)*(default_val-mean1)*(P.num_structures-checked_columns);
    std2 += (default_val-mean2)*(default_val-mean2)*(P.num_structures-checked_columns);
    //std::cout << val2 << " " << mean2 << " " << ((val2-mean2)*(val2-mean2)) << std::endl;

    return cov/std::sqrt(std1)/std::sqrt(std2);
}

template <class STRUCTURE, class CONTEXT, class NUMBER>
void
PropertyMatrix<STRUCTURE,CONTEXT,NUMBER>::
rankify(
)
{
    for (CONTEXT c: contexts())
    {
        std::vector<STRUCTURE> ranks(_structures.begin(),_structures.end());

        StructureComparisonFunction<STRUCTURE,CONTEXT,NUMBER> f(this,&c);
        std::sort(ranks.begin(), ranks.end(), f);

        size_t i=0;

        while (i<ranks.size())
        {
            Value<NUMBER> v1 = get(ranks[i],c);

            if (v1.null)
            {
                i++;
                continue;
            }

            size_t last_tie = i;

            while (i+1<ranks.size())
            {
                Value<NUMBER> v2 = get(ranks[i+1],c);

                if (v1.null || v2.value>v1.value)
                {
                    break;
                }

                i++;
            }

            for (size_t j=last_tie; j<=i; j++)
            {
                set(ranks[j],c,((double)last_tie+i)/2+1);
            }

            i++;
        }
    }
}

}
}

#endif
