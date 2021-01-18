/**
 * A property matrix associates values to pairs (structure, context).
 *
 * Examples of structures are edges, or triangles, in a network, and examples of contexts are layers
 * in a multilayer network.
 *
 * Several generic summarization functions can be computed on a property matrix, e.g., to
 * obtain the amount of overlapping or the statistical correlation between different types
 * of structures in different contexts.
 *
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_DATASTRUCTURES_PROPERTYMATRIX_PROPERTYMATRIX_H_
#define UU_CORE_DATASTRUCTURES_PROPERTYMATRIX_PROPERTYMATRIX_H_

#include <unordered_map>
#include <unordered_set>
#include "core/attributes/Value.hpp"
#include "core/utils/Counter.hpp"

namespace uu {
namespace core {

/**
 * A property matrix is a view that associates a value to each structure
 * (e.g., node, pair of nodes, ...) in each context (e.g., layer).
 */
template <class STRUCTURE, class CONTEXT, class VALUE>
class
    PropertyMatrix
{

  public:

    typedef STRUCTURE struct_type;
    typedef CONTEXT context_type;

    /** number of structures in this matrix (e.g., actors, or edges) */
    const long num_structures;

    /** number of observation contexts, typically layers */
    const long num_contexts;

    /**
     * Creates a property matrix with a given number of columns and rows.
     * @param num_contexts
     * @param num_structures
     * @param default_value
     */
    PropertyMatrix(
        long num_structures,
        long num_contexts,
        VALUE default_value
    );

    /**
     * @param s a structure
     * @param c a context
     * @return the value associated to the input structure in the input context
     */
    Value<VALUE>
    get(
        const STRUCTURE& s,
        const CONTEXT& c
    ) const;

    /**
     * @param s a structure
     * @param c a context
     * @param v the value to be associated to the input structure in the input context
     */
    void
    set(
        const STRUCTURE& s,
        const CONTEXT& c,
        VALUE v
    );

    /**
     * Sets s in c as not available (NA)
     * @param s a structure
     * @param c a context
     */
    void
    set_na(
        const STRUCTURE& s,
        const CONTEXT& c
    );

    /**
     * number of "not available" values in c (NA)
     * @param c a context
     */
    long
    num_na(
        const CONTEXT& c
    ) const;


    /**
     * Changes the fields of the matrix replacing absolute values with their ranking for each context.
     * For example, a context (1.2, 5.4, 5) would be replaced by (1, 3, 2)
     */
    void
    rankify(
    );

    /**
     * @return a set of all contexts in the matrix
     */
    const
    std::unordered_set<CONTEXT>&
    contexts(
    ) const;

    /**
     * @return a set of all structures in the matrix
     */
    const
    std::unordered_set<STRUCTURE>&
    structures(
    ) const;

    /**
     * @return the default value assigned to a structure for a context where no explicit value has been specified.
     */
    VALUE
    get_default(
    ) const;


  private:

    std::unordered_set<STRUCTURE> _structures;

    std::unordered_set<CONTEXT> _contexts;

    std::unordered_map<CONTEXT, std::unordered_map<STRUCTURE, Value<VALUE> > > data;

    VALUE default_value;

    Counter<CONTEXT> na;

};

/* TEMPLATE CODE */

template <class STRUCTURE, class CONTEXT, class VALUE>
PropertyMatrix<STRUCTURE,CONTEXT,VALUE>::
PropertyMatrix(
    long num_structures,
    long num_contexts,
    VALUE default_value
) :
    num_structures(num_structures),
    num_contexts(num_contexts),
    default_value(default_value)
{

}

template <class STRUCTURE, class CONTEXT, class VALUE>
Value<VALUE>
PropertyMatrix<STRUCTURE,CONTEXT,VALUE>::
get(
    const STRUCTURE& s,
    const CONTEXT& c
) const
{
    if (data.count(c)==0)
    {
        return  Value<VALUE>(default_value,false);
    }

    if (data.at(c).count(s)==0)
    {
        return  Value<VALUE>(default_value,false);
    }

    return data.at(c).at(s);
}

template <class STRUCTURE, class CONTEXT, class VALUE>
void
PropertyMatrix<STRUCTURE,CONTEXT,VALUE>::
set(
    const STRUCTURE& s,
    const CONTEXT& c,
    VALUE v
)
{
    data[c][s] = Value<VALUE>(v,false);
    _contexts.insert(c); // @todo this might slow down the function significantly - check
    _structures.insert(s);
}

template <class STRUCTURE, class CONTEXT, class VALUE>
void
PropertyMatrix<STRUCTURE,CONTEXT,VALUE>::
set_na(
    const STRUCTURE& s,
    const CONTEXT& c
)
{
    Value<VALUE> v = get(s,c);

    if (!v.null)
    {
        na.inc(c);
    }

    data[c][s] =  Value<VALUE>(v.value,true);
    _contexts.insert(c); // @todo this might slow down the function significantly - check
    _structures.insert(s);
}

template <class STRUCTURE, class CONTEXT, class VALUE>
long
PropertyMatrix<STRUCTURE,CONTEXT,VALUE>::
num_na(
    const CONTEXT& c
) const
{
    return na.count(c);
}

template <class STRUCTURE, class CONTEXT, class VALUE>
const
std::unordered_set<CONTEXT>&
PropertyMatrix<STRUCTURE,CONTEXT,VALUE>::
contexts(
) const
{
    return _contexts;
}

template <class STRUCTURE, class CONTEXT, class VALUE>
const
std::unordered_set<STRUCTURE>&
PropertyMatrix<STRUCTURE,CONTEXT,VALUE>::
structures(
) const
{
    return _structures;
}


template <class STRUCTURE, class CONTEXT, class VALUE>
VALUE
PropertyMatrix<STRUCTURE,CONTEXT,VALUE>::
get_default(
) const
{
    return default_value;
}

// rankify is defined in summarization.h

} // namespace core
} // namespace uu

#endif /* UU_CORE_DATASTRUCTURES_PROPERTYMATRIX_H_ */
