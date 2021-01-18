#include <unordered_map> 

typedef struct Dico Dico;
struct Dico {
  std::unordered_map<int,int> dict;
  int last;
};
