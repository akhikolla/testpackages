#include <Rcpp.h>
#include "concatenate_lists.h"

using namespace Rcpp;

List concatenate_lists(
    List a,
    List b){

  int length_a = a.size();
  int length_b = b.size();

  List return_list(length_a + length_b);

  CharacterVector names_a = a.names();
  CharacterVector names_b = b.names();

  CharacterVector return_list_names(length_a + length_b);
  return_list.attr("names") = return_list_names;
  for(int i = 0; i < length_a; i++) {
    return_list[i] = a[i];
    return_list_names[i] = names_a[i];
  }

  for(int i = 0; i < length_b; i++) {
    return_list[length_a+i] = b[i];
    return_list_names[length_a+i] = names_b[i];
  }

  return(return_list);
}
