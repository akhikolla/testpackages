// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <stdio.h>
#include "convert.h"

using namespace Rcpp;
using namespace std;

/*  Swedish stemmer tring to remove inflectional suffixes */

wstring swedish_stemming(wstring& word) {
  int len = word.size() - 1;

    if (len > 3) {   /*  -s  genitive form */
if (word[len]=='s') {
  word.erase(len);
  len--;
}
    }

    if (len > 6) {   /*  -elser  -heten  */
if ((word[len]=='r') && (word[len-1]=='e') && (word[len-2]=='s') &&
    (word[len-3]=='l') && (word[len-4]=='e')) {
  word.erase(len-4);
  return(word);
}
if ((word[len]=='n') && (word[len-1]=='e') && (word[len-2]=='t') &&
    (word[len-3]=='e') && (word[len-4]=='h')) {
  word.erase(len-4);
  return(word);
}
    }  /* len > 6 */


if (len > 5) {   /*  -arne  -erna  -ande  -else  -aste  -orna  -aren  */
if ((word[len]=='e') && (word[len-1]=='n') && (word[len-2]=='r') &&
    (word[len-3]=='a')) {
  word.erase(len-3);
  return(word);
}
if ((word[len]=='a') && (word[len-1]=='n') && (word[len-2]=='r') &&
    (word[len-3]=='e')) {
  word.erase(len-3);
  return(word);
}
if ((word[len]=='e') && (word[len-1]=='d') && (word[len-2]=='n') &&
    (word[len-3]=='a')) {
  word.erase(len-3);
  return(word);
}
if ((word[len]=='e') && (word[len-1]=='s') && (word[len-2]=='l') &&
    (word[len-3]=='e')) {
  word.erase(len-3);
  return(word);
}
if ((word[len]=='e') && (word[len-1]=='t') && (word[len-2]=='s') &&
    (word[len-3]=='a')) {
  word.erase(len-3);
  return(word);
}
if ((word[len]=='a') && (word[len-1]=='n') && (word[len-2]=='r') &&
    (word[len-3]=='o')) {
  word.erase(len-3);
  return(word);
}
if ((word[len]=='n') && (word[len-1]=='e') && (word[len-2]=='r') &&
    (word[len-3]=='a')) {
  word.erase(len-3);
  return(word);
}
}  /* len > 5 */


if (len > 4) {   /*  -are  comparative form */
if ((word[len]=='e') && (word[len-1]=='r') && (word[len-2]=='a')) {
  word.erase(len-2);
  return(word);
}
/* -ast  superlative form */
if ((word[len]=='t') && (word[len-1]=='s') && (word[len-2]=='a')) {
  word.erase(len-2);
  return(word);
}
/* -het  form */
if ((word[len]=='t') && (word[len-1]=='e') && (word[len-2]=='h')) {
  word.erase(len-2);
  return(word);
}
} /* if len > 4 */


if (len > 3) {
  /* -{aeo}r    */
  if ((word[len]=='r') &&
  ((word[len-1]=='a') || (word[len-1]=='o') || (word[len-1]=='e'))) {
    word.erase(len-1);
    return(word);
  }
  /* -en   */
  if ((word[len-1]=='e') && (word[len]=='n')) {
    word.erase(len-1);
    return(word);
  }
  /* -at  */
  if ((word[len-1]=='a') && (word[len]=='t')) {
    word.erase(len-1);
    return(word);
  }
  /* -te  */
  if ((word[len-1]=='t') && (word[len]=='e')) {
    word.erase(len-1);
    return(word);
  }
  /* -et  */
  if ((word[len-1]=='t') && (word[len]=='e')) {
    word.erase(len-1);
    return(word);
  }
}  /* end if len > 3 */


  if (len > 2) { /* -{taen}  */
  if ((word[len]=='t') || (word[len]=='a') ||
  (word[len]=='e') || (word[len]=='n')) {
    word.erase(len);
    return(word);
  }
  }  /* end if len > 2 */


  return(word);
}

//' Stem Swedish words
//'
//' Stemmer for Swedish words
//'
//' @param words a [character] containing the original words.
//' @return [character] with stemmed words.
//' @examples
//' swedish_stemmer(c("stiga"))
//' @export
// [[Rcpp::export]]
CharacterVector swedish_stemmer(Rcpp::StringVector words) {

  CharacterVector result(words.size());

  for (int i = 0; i < words.size(); ++i) {
    string s1 = static_cast<string>(words[i]);
    wstring str2 = utf8_to_utf16(s1);
    result[i] = swedish_stemming(str2);
    Rcpp::checkUserInterrupt();
  }

  return result;
}
