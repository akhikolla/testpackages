#ifndef FSELECTOR_ENTROPY_H
#define FSELECTOR_ENTROPY_H

#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <utility>
#include "support/table.h"
#include <numeric>

namespace fselector
{

namespace entropy
{


template<class InputIterator> double freq_entropy(InputIterator first, InputIterator last)
{

  const double sum = std::accumulate(first, last, 0.0,
                              [](const double& a, typename std::iterator_traits<InputIterator>::value_type iter)
                              {
                                return a + iter.second;
                              });

  const double result = std::accumulate(first, last, 0.0, [sum](const double& a, typename std::iterator_traits<InputIterator>::value_type iter)
                            {
                              if(iter.second > 0)
                              {
                                const double res = iter.second/sum;
                                return a + res * std::log(res);
                              } else
                              {
                                return a;
                              }

                            });



  return -result;
}

//// numeric entropy
template<class InputIterator> double numeric_entropy(InputIterator first, InputIterator last)
{

  const double sum = std::accumulate(first, last, 0.0);

  const double result = std::accumulate(first, last, 0.0, [sum](const double& a, typename std::iterator_traits<InputIterator>::value_type iter)
  {
    if(iter > 0)
    {
      const double res = iter/sum;
      return a + res * std::log(res);
    } else
    {
      return a;
    }

  });



  return -result;
}

template<class InputIterator> double entropy1d(InputIterator first, InputIterator last)
{
  const auto table = fselector::support::table1d(first, last);
  return freq_entropy(table.begin(), table.end());
}

template<typename T> class RollEntropy
{
  std::unordered_map<T, std::pair<int, double> > _map;
  int _size;
  double _sizeLog;

  public:

    RollEntropy() : _size(0) {};
    template<class InputIterator> RollEntropy(InputIterator first, InputIterator last) :
      _size(0)
    {
      for(; first != last; first++)
      {
        auto mit = _map.find(*first);

        if(mit != _map.end())
        {
          mit->second.first++;

        }
        else
        {
          _map[*first] = std::make_pair(1, 0.0);
        }

        _size++;
      }

      for(auto& mit : _map)
      {
        mit.second.second = std::log(mit.second.first);
      }

      _sizeLog = std::log(_size);

    }


    void add_sample(const T& val)
    {
      _size++;
      _sizeLog = std::log(_size);

      auto mit = _map.find(val);

      if(mit != _map.end())
      {
        mit->second.first++;
        mit->second.second = std::log(mit->second.first);
      }
      else
      {
        _map[val] = std::make_pair(1, 0.0);
      }
    }

    void remove_sample(const T& val)
    {
      auto mit = _map.find(val);
      if(mit != _map.end())
      {
        mit->second.first--;
        mit->second.second = std::log(mit->second.first);

        _size--;
        _sizeLog = _size > 0 ? std::log(_size) : 0.0;
      }
      else
      {
        _map[val] = std::make_pair(0, 0.0);
      }
    }

    double get_entropy()
    {
      double total = 0.0;
      for(const auto& it : _map)
      {
        if(it.second.first != 0)
        {
          const double res = it.second.second - _sizeLog;
          total += double(it.second.first)/double(_size) * res;
        }
      }

      return -total;
    }

};


class RollIntEntropy
{
  int _size;
  std::vector<int> _values;

public:

  RollIntEntropy(size_t levelsNo) : _size(0),
  _values(levelsNo + 1) {};
  template<class InputIterator> RollIntEntropy(InputIterator first, InputIterator last, size_t levelsNo) :
    _size(0),
    _values(levelsNo + 1)
  {
    for(; first != last; first++)
    {
      add_sample(*first);
      _size++;
    }
  }


  void add_sample(const int& val)
  {
    _size++;
    _values[val]++;
  }

  void remove_sample(const int& val)
  {
    _size--;
    _values[val]--;
  }

  double get_entropy()
  {
    double total = 0.0;
    double size = _size;
    for(const auto& it : _values)
    {
      if(it > 0)
      {
        const double res = it/size;
        total += res * std::log(res);
      }
    }

    return total;
  }

};



}
}


#endif
