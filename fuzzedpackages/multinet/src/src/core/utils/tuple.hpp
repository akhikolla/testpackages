/**
 * Utility functions to handle tuples, waiting for C++17 adoption.
 */

#ifndef UU_CORE_UTILS_TUPLE_H_
#define UU_CORE_UTILS_TUPLE_H_

#include <type_traits>
#include <tuple>
#include <cstddef>
#include <iostream> // cout, endl

namespace uu {
namespace core {


template <size_t N1, size_t N2, typename Func, size_t I = 0>
typename std::enable_if<I == N1*N2>::type
for_each_in_seq(Func) {}

template <size_t N1, size_t N2, typename Func, size_t I = 0>
typename std::enable_if<I < N1*N2>::type
for_each_in_seq(Func func)
{
    const int i = I / N2;
    const int j = I % N2;
    func.template go<N1, N2, i, j>();
    for_each_in_seq<N1, N2, Func, I + 1>(func);
}

// an example
struct print_pair
{
    template <size_t N1, size_t N2, size_t I1, size_t I2>
    void
    go ()
    {
        std::cout << "(" << I1 << "," << I2 << ") ";

        if (I2==N2-1)
        {
            std::cout << std::endl;
        }
    }
};



// From StackOverflow
// https://stackoverflow.com/questions/16387354/template-tuple-calling-a-function-on-each-element

template<size_t I = 0, typename Func, typename ...Ts>
typename std::enable_if<I == sizeof...(Ts)>::type
for_each_in_tuple(std::tuple<Ts...> &, Func) {}

template<size_t I = 0, typename Func, typename ...Ts>
typename std::enable_if<I < sizeof...(Ts)>::type
for_each_in_tuple(std::tuple<Ts...> & tpl, Func func)
{
    func(std::get<I>(tpl));
    for_each_in_tuple<I + 1>(tpl,func);
}


// an example
struct print_tuple_element
{
    template<typename T>
    void
    operator () (T&& t)
    {
        std::cout << t << std::endl;
    }
};

}
}



#endif /* UU_CORE_UTILS_HASHING_H_ */
