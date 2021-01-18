#ifndef ALGO_IT_H
#define ALGO_IT_H

template <class InputIterator, class Predicate>
InputIterator find_if_it(InputIterator it, InputIterator end,
                         Predicate pred) {
  while (it != end && !pred(it)) ++it;
  return it;
}

template <class InputIterator, class OutputIterator, class Predicate>
OutputIterator remove_copy_if_it(InputIterator it, InputIterator end,
                                 OutputIterator result, Predicate pred) {
  while (it != end) {
    if (!pred(it)) {
      *result = *it;
      ++result;
    }
    ++it;
  }
  return result;
}

template <class ForwardIterator, class Predicate>
ForwardIterator remove_if_it(ForwardIterator it, ForwardIterator end,
                             Predicate pred) {
  it = find_if_it(it, end, pred);
  ForwardIterator next = it;
  return it == end ? it : remove_copy_if_it(++next, end, it, pred);
}

#endif
