// Copyright Timothy H. Keitt 2017
// By use of this header, you agree to the following terms:
// 1. Your use conforms to either
//     a) the Boost Software License http://www.boost.org/users/license.html or
//     b) the license dictated by the R package https://github.com/thk686/strider and
// 2. You acknowledge use of this software in your product by either
//     a) place a visible notice somewhere in your product or
//     b) cite the software as you would any R package and
// 3. You do not remove these notices from this file

#ifndef __STRIDER_H__
#define __STRIDER_H__

#include <vector>
#include <functional>

#include <boost/iterator/iterator_adaptor.hpp>

namespace strider {
namespace detail {

using boost::iterator_adaptor;
using boost::enable_if_convertible;
using boost::iterator_core_access;

using std::next;
using std::size_t;
using std::vector;
using std::iterator_traits;

template <class Iterator>
class strided_iterator
  : public iterator_adaptor<strided_iterator<Iterator>, Iterator>
{
  using super_t = iterator_adaptor<strided_iterator<Iterator>, Iterator>;
  
  friend class iterator_core_access;
  
public:
  using difference_type = typename super_t::difference_type;
  using value_type = typename super_t::value_type;
  using pointer = typename super_t::pointer;
  using reference = typename super_t::reference;
  using iterator_category = typename super_t::iterator_category;
  
  strided_iterator() {}
  
  explicit strided_iterator(Iterator x, difference_type stride)
    : super_t(x), m_stride(stride) {}
  
  template<class OtherIterator>
  strided_iterator(
    strided_iterator<OtherIterator> const& r,
    typename enable_if_convertible<OtherIterator, Iterator>::type* = 0) 
    : super_t(r.base()), m_stride(r.m_stride) {}
  
private:
  reference dereference() const
  {
    return *this->base_reference();
  }
  
  void increment() { std::advance(this->base_reference(),  m_stride); }

  void decrement() { std::advance(this->base_reference(), -m_stride); }
  
  void advance(difference_type n)
  {
    std::advance(this->base_reference(), n * m_stride);
  }
  
  template <class OtherIterator>
  difference_type
  distance_to(strided_iterator<OtherIterator> const& y) const
  {
    return std::distance(this->base_reference(), y.base()) / m_stride;
  }
  
  template <class OtherIterator>
  bool equal(strided_iterator<OtherIterator> const& y) const
  {
    return this->base_reference() == y.base();
  }
  
  difference_type m_stride;
};

template<typename T>
inline strided_iterator<T>
make_strided(T iter,
             typename iterator_traits<T>::difference_type stride = 0,
             typename iterator_traits<T>::difference_type strides = 0)
{
  return strided_iterator<T>(next(iter, stride * strides), stride);
}

template<typename T>
class strided_range
{
public:
  using difference_type = typename iterator_traits<T>::difference_type;
  strided_range(T iter, difference_type stride, difference_type strides)
    : m_iter(iter), m_stride(stride), m_strides(strides) {}
  strided_iterator<T> begin() const
  {
    return make_strided(m_iter, m_stride);
  }
  strided_iterator<T> end() const
  {
    return make_strided(m_iter, m_stride, m_strides);
  }
private:
  T m_iter;
  difference_type m_stride, m_strides;
};

template<typename T>
inline strided_range<T>
make_strided_range(T iter,
                   typename iterator_traits<T>::difference_type stride,
                   typename iterator_traits<T>::difference_type strides)
{
  return strided_range<T>(iter, stride, strides);
}

}; // namespace detail

using detail::make_strided;
using detail::strided_range;
using detail::make_strided_range;
using detail::strided_iterator;

}; // namespace strider

#endif // __STRIDER_H__
