// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <stdio.h>
#include "convert.h"
using namespace Rcpp;
using namespace std;

/*  Portuguese stemmer tring to remove inflectional suffixes for nouns and adjectives */

static wstring removeAllPTAccent(wstring& word) {
  int len = word.size() - 1;
  int i;

  for(i=len; i>=0; i--) {
    if (word[i] == u'ä') {
      word[i] = u'a';
    }
    if (word[i] == u'â') {
      word[i] = u'a';
    }
    if (word[i] == u'à') {
      word[i] = u'a';
    }
    if (word[i] == u'á') {
      word[i] = u'a';
    }
    if (word[i] == u'ã') {
      word[i] = u'a';
    }
    if (word[i] == u'ê') {
      word[i] = u'e';
    }
    if (word[i] == u'é') {
      word[i] = u'e';
    }
    if (word[i] == u'è') {
      word[i] = u'e';
    }
    if (word[i] == u'ë') {
      word[i] = u'e';
    }
    if (word[i] == u'ï') {
      word[i] = u'i';
    }
    if (word[i] == u'î') {
      word[i] = u'i';
    }
    if (word[i] == u'ì') {
      word[i] = u'i';
    }
    if (word[i] == u'í') {
      word[i] = u'i';
    }
    if (word[i] == u'ü') {
      word[i] = u'u';
    }
    if (word[i] == u'ú') {
      word[i] = u'u';
    }
    if (word[i] == u'ù') {
      word[i] = u'u';
    }
    if (word[i] == u'û') {
      word[i] = u'u';
    }
    if (word[i] == u'ô') {
      word[i] = u'o';
    }
    if (word[i] == u'ö') {
      word[i] = u'o';
    }
    if (word[i] == u'ó') {
      word[i] = u'o';
    }
    if (word[i] == u'ò') {
      word[i] = u'o';
    }
    if (word[i] == u'õ') {
      word[i] = u'o';
    }
    if (word[i] == u'ç') {
      word[i] = u'c';
    }
  }
  return(word);
}

static wstring finalVowelPortuguese(wstring& word) {
  int len = word.size() - 1;

  if (len > 3) {
    if ((word[len]== u'e') || (word[len]== u'a') || (word[len]== u'o')) {
      word.erase(len);  /* remove final -e or -a or -o */
return(word);
    }
  }
  return(word);
}


/* Remove plural and feminine form of Portuguese words */

static wstring remove_PTsuffix (wstring& word) {
  int len = word.size() - 1;

  if (len > 3) {   /* plural in -es when sing form ends with -r, -s, -l or -z*/
if ((word[len]== u's') && (word[len-1]== u'e') &&
    ((word[len-2]== u'r') || (word[len-2]== u's') ||
    (word[len-2]== u'z') || (word[len-2]== u'l'))) {
  word.erase(len-1);  /* doutores (plur) --> doutor (sing) */
return(word);
}
  }  /* len > 3 */

if (len > 2) {   /* when sing form ends with -em, change -m in -n in plur */
if ((word[len]== u's') && (word[len-1]== u'n')) {
  word[len-1]= u'm';     /* homens (plur) --> homem (sing) */
word.erase(len);
return(word);
}
} /* len > 2 */

if (len > 3) {   /* when sing form ends with -el, change -el in -eis  */
if (((word[len]== u's') && (word[len-1]== u'i')) &&
    ((word[len-2]== u'e') || (word[len-2]== u'é'))) {
  word[len-2]= u'e';     /* papéis (plur) --> papel (sing) */
word[len-1]= u'l';     /* error:  faceis (plur) --> facil (sing) */
word.erase(len);
return(word);
}
} /* len > 3 */

if (len > 3) {   /* when sing form ends with -ais, change -ais in -al in plur */
if ((word[len]== u's') && (word[len-1]== u'i') && (word[len-2]== u'a')) {
  word[len-1]= u'l';     /* normais (plur) --> normal (sing) */
word.erase(len);
return(word);
}
} /* len > 3 */

if (len > 3) {   /* when sing form ends with -'ois, change -ais in -al in plur */
if ((word[len]== u's') && (word[len-1]== u'i') && (word[len-2]== u'ó')) {
  word[len-2]= u'o';     /* lencois (plur) --> lencol (sing) */
word[len-1]= u'l';
word.erase(len);
return(word);
}
} /* len > 3 */

if (len > 3) {   /* when sing form ends with -is, change -is in -il in plur */
if ((word[len]== u's') && (word[len-1]== u'i')) {
  word[len]= u'l';     /* barris (plur) --> barril (sing) */
return(word);
}
} /* len > 3 */

if (len > 2) {   /* when plur form ends with -ões, change -ões in -ão  */
if ((word[len]== u's') && (word[len-1]== u'e') && (word[len-2]== u'õ')) {
  word[len-2]= u'ã';     /* botões (plur) --> botão (sing) */
word[len-1]= u'o';
word.erase(len);
return(word);
}         /* when plur form ends with -ães, change -ães in -ão  */
if ((word[len]== u's') && (word[len-1]== u'e') && (word[len-2]== u'ã')) {
  word[len-1]= u'o';     /* caes (plur) --> cao (sing) */
word.erase(len);
return(word);
}
} /* len > 2 */

if (len > 5) {   /* for adverb -mente */
if ((word[len]== u'e') && (word[len-1]== u't') && (word[len-2]== u'n') &&
    (word[len-3]== u'e') && (word[len-4]== u'm')) {
  word.erase(len-4);
  return(word);
}
} /* len > 5 */

if (len > 2) {   /* usually plural in -s */
if (word[len]== u's') {
  word.erase(len);   /* of course, invariable word, pires->pires */
len--;
}
} /* len > 2 */

return(word);
}


static wstring normFemininPortuguese(wstring& word) {
  int len = word.size() - 1;

  if ((len < 3) || (word[len]!= u'a')) {
    return(word);
  }
  if (len > 6) {
    /*  -inha  -> inho */
    if ((word[len-1]== u'h') && (word[len-2]== u'n') && (word[len-3]== u'i')) {
      word[len]= u'o';
      return(word);
    }
    /*  -iaca  -> iaco */
    if ((word[len-1]== u'c') && (word[len-2]== u'a') && (word[len-3]== u'i')) {
      word[len]= u'o';
      return(word);
    }
    /*  -eira  -> eiro */
    if ((word[len-1]== u'r') && (word[len-2]== u'i') && (word[len-3]== u'e')) {
      word[len]= u'o';
      return(word);
    }
  } /* len > 6 */

    if (len > 5) {
      /*  -ona  -> ão */
      if ((word[len-1]== u'n') && (word[len-2]== u'o')) {
        word[len-2]= u'ã';
        word[len-1]= u'o';
        word.erase(len);
        return(word);
      }
      /*  -ora  -> or */
      if ((word[len-1]== u'r') && (word[len-2]== u'o')) {
        word.erase(len);
        return(word);
      }
      /*  -osa  -> oso */
      if ((word[len-1]== u's') && (word[len-2]== u'o')) {
        word[len]= u'o';
        return(word);
      }
      /*  -esa  -> ês */
      if ((word[len-1]== u's') && (word[len-2]== u'e')) {
        word[len-2]= u'ê';
        word.erase(len);
        return(word);
      }
      /*  -ica  -> ico */
      if ((word[len-1]== u'c') && (word[len-2]== u'i')) {
        word[len]= u'o';
        return(word);
      }
      /*  -ida  -> ido */
      if ((word[len-1]== u'd') && (word[len-2]== u'i')) {
        word[len]= u'o';
        return(word);
      }
      /*  -ada  -> ido */
      if ((word[len-1]== u'd') && (word[len-2]== u'a')) {
        word[len]= u'o';
        return(word);
      }
      /*  -iva  -> ivo */
      if ((word[len-1]== u'v') && (word[len-2]== u'i')) {
        word[len]= u'o';
        return(word);
      }
      /*  -ama  -> imo */
      if ((word[len-1]== u'm') && (word[len-2]== u'a')) {
        word[len]= u'o';
        return(word);
      }
      /*  -na  -> no */
      if (word[len-1]== u'n') {
        word[len]= u'o';
        return(word);
      }
    } /* len > 5 */

      return(word);
}

wstring portuguese_stemming (wstring& word) {
  int len = word.size() - 1;

  if (len > 2) {
    remove_PTsuffix(word);
    normFemininPortuguese(word);
    finalVowelPortuguese(word);
    removeAllPTAccent(word);
  }
  return(word);
}

//' Stem Portuguese words
//'
//' Stemmer for Portuguese words
//'
//' @param words a [character] containing the original words.
//' @return [character] with stemmed words.
//' @examples
//' portuguese_stemmer(c("adeus"))
//' @export
// [[Rcpp::export]]
CharacterVector portuguese_stemmer(Rcpp::StringVector words) {
  CharacterVector result(words.size());

  for (int i = 0; i < words.size(); ++i) {
    string s1 = static_cast<string>(words[i]);
    wstring str2 = utf8_to_utf16(s1);
    result[i] = portuguese_stemming(str2);
    Rcpp::checkUserInterrupt();
  }

  return result;
}
