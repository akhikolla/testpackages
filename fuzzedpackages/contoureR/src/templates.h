#include <vector>

using namespace std;

#ifndef __templates__
#define __templates__

//Template Functions
template<typename T>
void uniqueOnly(vector<T>& vec){
  typedef typename std::vector<T>::iterator It;
  It it = std::unique(vec.begin(), vec.end());
  //vec.resize( std::distance(vec.begin(),it) );
  vec.erase(it,vec.end());
}

template<typename T>
bool isIn(const vector<T>& x,const T v){
  return std::find(x.begin(), x.end(), v) != x.end();
}

#endif /* __templates__ */
