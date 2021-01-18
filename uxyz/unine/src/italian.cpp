// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <stdio.h>
#include "convert.h"
using namespace Rcpp;
using namespace std;

/*  Italian stemmer tring to remove inflectional suffixes */


static void removeItalianAccent(wstring& word) {
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
}

static wstring italian_stemming(wstring word) {
  int len = word.size() - 1;

  if (len > 4) {
    removeItalianAccent(word);
    if (word[len]==u'e') {  /*  ending with -ie or -he  */
if (word[len-1]==u'i' || word[len-1]==u'h') {
  word.erase(len-1);
  return (word);
}
word.erase(len);  /*  ending with -e  */
return(word);
    }
    if (word[len]==u'i') {  /*  ending with -hi or -ii */
if ((word[len-1]==u'h') || (word[len-1]==u'i')) {
  word.erase(len-1);
  return (word);
}
word.erase(len);  /*  ending with -i  */
return(word);
    }
    if (word[len]==u'a') {  /*  ending with -ia  */
if (word[len-1]==u'i') {
  word.erase(len-1);
  return (word);
}
word.erase(len);  /*  ending with -a  */
return(word);
    }
    if (word[len]==u'o') {  /*  ending with -io  */
if (word[len-1]==u'i') {
  word.erase(len-1);
  return (word);
}
word.erase(len);  /*  ending with -o  */
return(word);
    }

  } /* end if (len > 4) */
return(word);
}

//' Stem Italian words
//'
//' Stemmer for Italian words
//'
//' @param words a [character] containing the original words.
//' @return [character] with stemmed words.
//' @examples
//' italian_stemmer(c("arrivederci"))
//' @export
// [[Rcpp::export]]
CharacterVector italian_stemmer(Rcpp::StringVector words) {
  CharacterVector result(words.size());

  for (int i = 0; i < words.size(); ++i) {
    string s1 = static_cast<string>(words[i]);
    wstring str2 = utf8_to_utf16(s1);
    result[i] = italian_stemming(str2);
    Rcpp::checkUserInterrupt();
  }

  return result;
}
