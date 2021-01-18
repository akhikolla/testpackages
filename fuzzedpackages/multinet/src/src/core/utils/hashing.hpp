/**
 * Allows more data types to be used as keys of hash-based indexes.
 */

#ifndef UU_CORE_UTILS_HASHING_H_
#define UU_CORE_UTILS_HASHING_H_

#include <set>
#include <functional>

namespace uu {
namespace core {

// From StackOverflow - itself from boost.
template <class T>
inline
void
hash_combine(
    std::size_t& seed,
    const T& v
)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

}
}

namespace std {

template <class T1, class T2>
struct hash<std::pair<T1,T2>>
{
    std::size_t
    operator () (const std::pair<T1,T2> &p) const
    {
        size_t res;
        uu::core::hash_combine(res, p.first);
        uu::core::hash_combine(res, p.second);
        return res;
    }
};

template <class T>
struct hash<std::set<T> >
{
    std::size_t
    operator () (const std::set<T> &s) const
    {
        size_t res;

        for (T el: s)
        {
            uu::core::hash_combine(res, el);
        }

        return res;
    }
};
}


#endif /* UU_CORE_UTILS_HASHING_H_ */
