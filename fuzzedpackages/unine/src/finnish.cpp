// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <stdio.h>
#include "convert.h"
using namespace Rcpp;
using namespace std;

/*  Finnish stemmer tring to remove inflectional suffixes */


bool IsVowel(const wchar_t& c) {
  bool result = (u'a'==(c)||u'e'==(c)||u'i'==(c)||u'o'==(c)||u'u'==(c)||\
                 u'y'==(c));
  return(result);
}


static void removeFinnishAccent(wstring& word) {
  int len = word.size() - 1;
  int i;

  for(i=len; i>=0; i--) {
    if (word[i] == u'ä') {
      word[i] = u'a';
    }
    if (word[i] == u'å') {
      word[i] = u'a';
    }
    if (word[i] == u'ö') {
      word[i] = u'o';
    }
  }
}


static wstring norm_finnish(wstring& word) {
  int len = word.size() - 1;

  if (len > 4) {   /* -hde  -> -ksi  */
if ((word[len]== u'e') && (word[len-1]== u'd') && (word[len-2]== u'h')) {
  word[len-2]= u'k';
  word[len-1]= u's';
  word[len]= u'i';
}
  }  /* end if len > 4 */
if (len > 3) {   /* -ei  -> -  */
if ((word[len]== u'i') && (word[len-1]== u'e')) {
  word.erase(len-1);
  return(word);
}
if ((word[len]== u't') && (word[len-1]== u'a')) {
  word.erase(len-1);
  return(word);
}
}  /* end if len > 3 */
if (len > 2) {   /* plural    -t  or -(aeiouy)i */
if ((word[len]== u't') || (word[len]== u's') || (word[len]== u'j')
      || (word[len]== u'e') || (word[len]== u'a')) {
  word.erase(len);
}
else {
  /*      if ((word[len]== u'i') && (IsVowel(word[len-1]))) { */
  if ((word[len]== u'i')) {
    word.erase(len);
  }
}
} /* end if (len > 2) */
  return(word);
}

static void removeDoubleKPT(wstring& word) {
  int len = word.size() - 1;
  int i, position;
  char currentChar;

  if (len > 3) {
    currentChar = word[0];
    position = 1;
    while (word[position]) {  /*  remove double kk pp tt  */
  if ((currentChar==word[position]) && ((currentChar== u'k') ||
      (currentChar== u'p') || (currentChar== u't'))) {
    i = position-1;
    while (word[i] != u'\0') {
      word[i] = word[i+1];
      i++;
    }
  }  /* end if */
  else {
    currentChar = word[position];
    position++;
  }
    }  /* end while */
  } /* end if len */
}

static wstring norm2_finnish (wstring& word) {
  int len = word.size() - 1;

  if (len > 7) {   /* -e, -o,  -u */
  if ((word[len]== u'e') || (word[len]== u'o') || (word[len]== u'u')) {
    word.erase(len);
    len--;
  }
  }
  if (len > 3) {   /* plural    -i  */
  if (word[len]== u'i') {
    word.erase(len);
  }
  removeDoubleKPT(word);
  } /* end if (len > 3) */
  return(word);
}

static wstring finnishStep1 (wstring&  word) {
  int len = word.size() - 1;

  if (len > 7) {
    /*    -kin  */
    if ((word[len]== u'n') && (word[len-1]== u'i') && (word[len-2]== u'k')) {
      word.erase(len-2);
      return(finnishStep1(word));
    }
    /*    -ko  */
    if ((word[len]== u'o') && (word[len-1]== u'k')) {
      word.erase(len-1);
      return(finnishStep1(word));
    }
  } /* end if (len > 7) */
    if (len > 10) {
      /*    -dellinen  for adjective  */
      if ((word[len]== u'n') && (word[len-1]== u'e') && (word[len-2]== u'n')
            && (word[len-3]== u'i') && (word[len-4]== u'l') && (word[len-5]== u'l')
            && (word[len-6]== u'e') && (word[len-7]== u'd')) {
            word.erase(len-7);
        return(word);
      }
      /*    -dellisuus  for adverb  */
      if ((word[len]== u's') && (word[len-1]== u'u') && (word[len-2]== u'u')
            && (word[len-3]== u's') && (word[len-4]== u'i') && (word[len-5]== u'l')
            && (word[len-6]== u'l') && (word[len-7]== u'e') && (word[len-8]== u'd')) {
            word.erase(len-8);
        return(word);
      }
    } /* end if (len > 10) */
      return(word);
}


static wstring finnishStep2(wstring& word) {
  int len = word.size() - 1;

  if (len > 4) {
    /*    -lla  */
    if ((word[len]== u'a') && (word[len-1]== u'l') && (word[len-2]== u'l')) {
      word.erase(len-2);
      return(word);
    }
    /*    -tse  */
    if ((word[len]== u'e') && (word[len-1]== u's') && (word[len-2]== u't')) {
      word.erase(len-2);
      return(word);
    }
    /*    -sti  */
    if ((word[len]== u'i') && (word[len-1]== u't') && (word[len-2]== u's')) {
      word.erase(len-2);
      return(word);
    }
    /*    -ni  */
    if ((word[len]== u'i') && (word[len-1]== u'n')) {
      word.erase(len-1);
      return(word);
    }
    /*    -a  if -aa  */
    if ((word[len]== u'a') && (word[len-1]== u'a')) {
      word.erase(len);
      return(word);
    }
  } /* end if (len > 4) */
    return(word);
}


static wstring finnishStep3 (wstring& word) {
  int len = word.size() - 1;

  if (len > 7) {
    /* genetive -nnen  -s  */
    if ((word[len]== u'n') && (word[len-1]== u'e') &&
    (word[len-2]== u'n') && (word[len-3]== u'n')) {
      word[len-3]= u's';
      word.erase(len-2);
      return(word);
    }
    /* essive -ntena  -s  */
    if ((word[len]== u'a') && (word[len-1]== u'n') && (word[len-2]== u'e') &&
    (word[len-3]== u't') && (word[len-4]== u'n')) {
      word[len-4]= u's';
      word.erase(len-3);
      return(word);
    }
    /*  -tten  -s  */
    if ((word[len]== u'n') && (word[len-1]== u'e') && (word[len-2]== u't') &&
    (word[len-3]== u't')) {
      word.erase(len-3);
      return(word);
    }
    /* genitive plural   -eiden  -s  */
    if ((word[len]== u'n') && (word[len-1]== u'e') &&
    (word[len-2]== u'd') && (word[len-3]== u'i') && (word[len-4]== u'e')) {
      word.erase(len-4);
      return(word);
    }
  }
  if (len > 5) {
    /* komitatiivi plural   -neen  */
    if ((word[len]== u'n') && (word[len-1]== u'e') &&
    (word[len-2]== u'e') && (word[len-3]== u'n')) {
      word.erase(len-3);
      return(word);
    }
    /* illatiivi   -siin,  etc.  */
    if ((word[len]== u'n') && (word[len-1]== u'i') &&
    (word[len-2]== u'i') && (word[len-3]== u'n')) {
      word.erase(len-3);
      return(word);
    }
    /* illatiivi   -seen,  etc.  */
    if ((word[len]== u'n') && (word[len-1]== u'e') &&
    (word[len-2]== u'e') && (word[len-3]== u's')) {
      word.erase(len-3);
      return(word);
    }
    /* illatiivi   -h?n,  etc.  */
    if ((word[len]== u'n') && (IsVowel(word[len-1])) &&
    (word[len-2]== u'h')) {
      word.erase(len-2);
      return(word);
    }
    /* genitive plural   -teen,  etc.  */
    if ((word[len]== u'n') && (word[len-1]== u'e') &&
    (word[len-2]== u'e') && (word[len-3]== u't')) {
      word.erase(len-3);
      return(word);
    }
    /* genitive plural   -den  -s  */
    if ((word[len]== u'n') && (word[len-1]== u'e') &&
    (word[len-2]== u'd')) {
      word[len-2]= u's';
      word.erase(len-1);
      return(word);
    }
    /* genitive plural   -ksen  -s  */
    if ((word[len]== u'n') && (word[len-1]== u'e') &&
    (word[len-2]== u's') && (word[len-3]== u'k')) {
      word[len-3]= u's';
      word.erase(len-2);
      return(word);
    }
    /* komitatiivi plural   -neni  */
    if ((word[len]== u'n') && (word[len-1]== u'e') &&
    (word[len-2]== u'n') && (word[len-3]== u'i')) {
      word.erase(len-3);
      return(word);
    }
    /* inessiivi   -ssa  */
    if ((word[len]== u'a') && (word[len-1]== u's') && (word[len-2]== u's')) {
      word.erase(len-2);
      return(word);
    }
    /* elatiivi   -sta  */
    if ((word[len]== u'a') && (word[len-1]== u't') && (word[len-2]== u's')) {
      word.erase(len-2);
      return(word);
    }
    /* adessiivi   -lla  */
    if ((word[len]== u'a') && (word[len-1]== u'l') && (word[len-2]== u'l')) {
      word.erase(len-2);
      return(word);
    }
    /* ablatiivi   -lta  */
    if ((word[len]== u'a') && (word[len-1]== u't') && (word[len-2]== u'l')) {
      word.erase(len-2);
      return(word);
    }
    /* abessiivi   -tta  */
    if ((word[len]== u'a') && (word[len-1]== u't') && (word[len-2]== u't')) {
      word.erase(len-2);
      return(word);
    }
    /* translatiivi   -ksi  */
    if ((word[len]== u'i') && (word[len-1]== u's') && (word[len-2]== u'k')) {
      word.erase(len-2);
      return(word);
    }
    /* allatiivi   -lle  */
    if ((word[len]== u'e') && (word[len-1]== u'l') && (word[len-2]== u'l')) {
      word.erase(len-2);
      return(word);
    }
  } /* end if (len > 5) */
    if (len > 4) {
      /* essiivi   -na  */
      if ((word[len]== u'a') && (word[len-1]== u'n')) {
        word.erase(len-1);
        return(word);
      }
      /* komitatiivi   -ne-  */
      if ((word[len]== u'e') && (word[len-1]== u'n')) {
        word.erase(len-1);
        return(word);
      }
      if ((word[len]== u'i') && (word[len-1]== u'e') && (word[len-2]== u'n')) {
        word.erase(len-2);
        return(word);
      }
    } /* end if (len > 4) */
      if (len > 3) {
        /* partitiivi   -(t,j)a  */
        if ((word[len]== u'a') && ((word[len-1]== u't') || (word[len-1]== u'j'))) {
          word.erase(len-1);
          return(word);
        }
        if (word[len]== u'a') {
          word.erase(len);
          return(word);
        }
        /* illatiivi   -an, -en, -on, -in, -un, -yn, etc.  */
        if ((word[len]== u'n') && ((word[len-1]== u'a') || (word[len-1]== u'e')
                                   || (word[len-1]== u'o') || (word[len-1]== u'i')
                                   || (word[len-1]== u'u') || (word[len-1]== u'y'))) {
                                   word.erase(len-1);
          return(word);
        }
        /* genetiivi or instruktiivi   -n  */
        if (word[len]== u'n') {
          word.erase(len);
          return(word);
        }
      } /* end if (len > 3) */
        return(word);
}

wstring finnish_stemming(wstring& word) {
  int len = word.size() - 1;

  if (len > 2) {
    removeFinnishAccent(word);
    word = finnishStep1(word);
    word = finnishStep2(word);
    word = finnishStep3(word);
    word = norm_finnish(word);
    word = norm2_finnish(word);
  }
  return(word);
}

//' Stem Finnish words
//'
//' Stemmer for Finnish words
//'
//' @param words a [character] containing the original words.
//' @return [character] with stemmed words.
//' @examples
//' finnish_stemmer(c("taivas"))
//' @export
// [[Rcpp::export]]
CharacterVector finnish_stemmer(Rcpp::StringVector words) {
  CharacterVector result(words.size());

  for (int i = 0; i < words.size(); ++i) {
    string s1 = static_cast<string>(words[i]);
    wstring str2 = utf8_to_utf16(s1);
    result[i] = finnish_stemming(str2);
    Rcpp::checkUserInterrupt();
  }

  return result;
}
