/**
 * History:
 * - 2018.01.01 file created, adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_ATTRIBUTES_VALUE_H_
#define UU_CORE_ATTRIBUTES_VALUE_H_

#include <ostream>

namespace uu {
namespace core {

/**
 * This structure stores an attribute value associated to an object, or a flag
 * indicating that the value is missing (NA/NULL).
 */
template <class T>
struct Value
{
    /** The value, which is only valid if member `null` is false. */
    T value;

    /** Indicates whether member `value` is valid or not available (NA/NULL). */
    bool null;


    /** Constructs a NULL value. */
    Value<T>
    (
    ) : value(T()),
        null(true)
    {
    };


    /**
     * Constructs an object with an available value.
     * The «null« flag is automatically set to false.
     * @param v the value of this object.
     */
    Value<T>
    (
        const T& v
    ) :
        value(v),
        null(false)
    {
    };


    /**
     * Constructs an object allowing to specify a value and also if the object is available or NULL.
     * @param v the value of this object, to be used if member «null« is false
     * @param is_null a flag indicating if the value is available or not (NULL)
     */
    Value<T>
    (
        const T& v,
        bool is_null
    ) :
        value(v),
        null(is_null)
    {
    };


    template <class U>
    friend std::ostream&
    operator<<(std::ostream&, const Value<U>&);

};


/**
 * Prints an object of type Value to an output stream.
 * If `null` is true, "NA" (Not Available) is written to the stream
 * independently of the value stored in member `value`.
 */
template <class T>
std::ostream&
operator<<(
    std::ostream& os,
    const Value<T>& v
)
{
    if (v.null)
    {
        os << "NA";
    }

    else
    {
        os << v.value;
    }

    return os;
}

} // namespace core
} // namespace uu

#endif /* UU_CORE_ATTRIBUTES_VALUE_H_ */
