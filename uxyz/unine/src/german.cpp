// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <stdio.h>
#include "convert.h"
using namespace Rcpp;
using namespace std;


/*  German stemmer tring to remove inflectional suffixes */

static wstring removeGermanAccent(wstring& word) {
  int len = word.size() - 1;
  int i;

  for(i=len; i>=0; i--) {
    if ((word[i]==u'ä') || (word[i]==u'à') || (word[i]==u'á') || (word[i]==u'â')) {
      word[i] = 'a';
    }
    if ((word[i]==u'ö') || (word[i]==u'ò') || (word[i]==u'ó') || (word[i]==u'ô')) {
      word[i] = 'o';
    }
    if ((word[i]==u'ï') || (word[i]==u'ì') || (word[i]==u'í') || (word[i]==u'î')) {
      word[i] = 'i';
    }
    if ((word[i]==u'ü') || (word[i]==u'ù') || (word[i]==u'ú') || (word[i]==u'û')) {
      word[i] = 'u';
    }
  }
  return(word);
}


static int STEnding (wchar_t aLetter) {
  if (aLetter==u'b' || aLetter==u'd' || aLetter==u'f' ||
      aLetter==u'g' || aLetter==u'h' || aLetter==u'k' ||
      aLetter==u'l' || aLetter==u'm' || aLetter==u'n' ||
      aLetter==u't')
    return(1);
  return(0);
}


wstring remove_Step1 (wstring& word) {
  int len = word.size() - 1;

  if (len > 4) {
    if (word[len]==u'n' && word[len-1]==u'r' && word[len-2]==u'e') {
      word.erase(len-2);  /*  ending with -ern ->   */
return(word);
    }
  }
  if (len > 3) {
    if (word[len]==u'm' && word[len-1]==u'e') {
      word.erase(len-1);  /*  ending with -em ->  */
return(word);
    }
    if (word[len]==u'n' && word[len-1]==u'e') {
      word.erase(len-1);  /*  ending with -en ->  */
return(word);
    }
    if (word[len]==u'r' && word[len-1]==u'e') {
      word.erase(len-1);  /*  ending with -er ->  */
return(word);
    }
    if (word[len]==u's' && word[len-1]==u'e') {
      word.erase(len-1);  /*  ending with -es ->  */
return(word);
    }
  }
  if (len > 2) {
    if (word[len]==u'e') {
      word.erase(len);  /*  ending with -e ->  */
return(word);
    }
    if (word[len]==u's' && STEnding(word[len-1])) {
      word.erase(len);  /*  ending with -s ->  */
return(word);
    }
  }
  return(word);
}

wstring remove_Step2 (wstring& word) {
  int len = word.size() - 1;

  if (len > 4) {
    if (word[len]==u't' && word[len-1]==u's' && word[len-2]==u'e') {
      word.erase(len-2);  /*  ending with -est ->   */
return(word);
    }
  }
  if (len > 3) {
    if (word[len]==u'r' && word[len-1]==u'e') {
      word.erase(len-1);  /*  ending with -er ->  */
return(word);
    }
    if (word[len]==u'n' && word[len-1]==u'e') {
      word.erase(len-1);  /*  ending with -en ->  */
return(word);
    }
    if (word[len]==u't' && word[len-1]==u's' && STEnding(word[len-2])) {
      word.erase(len-1);  /*  ending with -st ->  */
return(word);
    }
  }
  return(word);
}



wstring german_stemming (wstring word) {
  removeGermanAccent(word);
  remove_Step1(word);
  remove_Step2(word);
  return(word);
}

//' Stem German words
//'
//' Stemmer for German words
//'
//' @param words a [character] containing the original words.
//' @return [character] with stemmed words.
//' @examples
//' german_stemmer(c("kinder"))
//' @export
// [[Rcpp::export]]
CharacterVector german_stemmer(Rcpp::StringVector words) {
  CharacterVector result(words.size());

  for (int i = 0; i < words.size(); ++i) {
    string s1 = static_cast<string>(words[i]);
    wstring str2 = utf8_to_utf16(s1);
    str2 = removeGermanAccent(str2);
    str2 = remove_Step1(str2);
    result[i] = remove_Step2(str2);

    Rcpp::checkUserInterrupt();
  }

  return result;
}
