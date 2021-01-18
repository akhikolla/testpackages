// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <stdio.h>
#include "convert.h"
using namespace Rcpp;
using namespace std;

/*  Spanish stemmer tring to remove inflectional suffixes */

wstring removeSpanishAccent (wstring& word) {
  int len = word.size() - 1;
  int i;

  for(i=len; i>=0; i--) {
    if ((word[i]==u'à') || (word[i]==u'á') || (word[i]==u'â') || (word[i]==u'ä')) {
      word[i] = 'a';
    }
    if ((word[i]==u'ò') || (word[i]==u'ó') || (word[i]==u'ô') || (word[i]==u'ö')) {
      word[i] = 'o';
    }
    if ((word[i]==u'è') || (word[i]==u'é') || (word[i]==u'ê') || (word[i]==u'ë')) {
      word[i] = 'e';
    }
    if ((word[i]==u'ù') || (word[i]==u'ú') || (word[i]==u'û') || (word[i]==u'ü')) {
      word[i] = 'u';
    }
    if ((word[i]==u'ì') || (word[i]==u'í') || (word[i]==u'î') || (word[i]==u'ï')) {
      word[i] = 'i';
    }
  }
  return(word);
}

wstring spanish_word_stemmer(wstring word) {
  int len = word.size() - 1;

  if (len > 3) {
    removeSpanishAccent(word);
    if ((word[len]==u's') && (word[len-1]==u'e') && (word[len-2]==u's') && (word[len-3]==u'e')) {
      /*  corteses -> cortés  */
      word.erase(len-1);
      return(word);
    }
    if ((word[len]==u's') && (word[len-1]==u'e') && (word[len-2]==u'c')) {
      word[len-2]='z';        /*  dos veces -> una vez  */
      word.erase(len-1);
      return(word);
    }
    if (word[len]==u's') {  /*  ending with -os, -as  or -es */
      if (word[len-1]==u'o' || word[len-1]==u'a' || word[len-1]==u'e' ) {
        word.erase(len-1);  /*  remove -os, -as  or -es */
      return (word);
      }
    }
    if (word[len]==u'o') {   /*  ending with  -o  */
      word.erase(len);
      return(word);
    }
    if (word[len]==u'a') {   /*  ending with  -a  */
      word.erase(len);
      return(word);
    }
    if (word[len]==u'e') {   /*  ending with  -e  */
      word.erase(len);
      return(word);
    }
  } /* end if (len > 3) */
      return(word);
}

//' Stem Spanish words
//'
//' Stemmer for Spanish words
//'
//' @param words a [character] containing the original words.
//' @return [character] with stemmed words.
//' @examples
//' spanish_stemmer(c("perros"))
//' @export
// [[Rcpp::export]]
CharacterVector spanish_stemmer(Rcpp::StringVector words) {
  CharacterVector result(words.size());

  for (int i = 0; i < words.size(); ++i) {
    string s1 = static_cast<string>(words[i]);
    wstring str2 = utf8_to_utf16(s1);
    result[i] = spanish_word_stemmer(str2);
    Rcpp::checkUserInterrupt();
  }
  return result;
}
