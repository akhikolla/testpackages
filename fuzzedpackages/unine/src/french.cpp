// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <stdio.h>
#include "convert.h"
using namespace Rcpp;
using namespace std;

static void removeAllFEAccent (wstring& word) {
    int len = word.size() -1;
    int i;

    for(i=len; i>=0; i--) {
      if (word[i] == u'â') {
        word[i] = 'a';
      }
      if (word[i] == u'à') {
        word[i] = 'a';
      }
      if (word[i] == u'á') {
        word[i] = 'a';
      }
      if (word[i] == u'ê') {
        word[i] = 'e';
      }
      if (word[i] == u'é') {
        word[i] = 'e';
      }
      if (word[i] == u'è') {
        word[i] = 'e';
      }
      if (word[i] == u'î') {
        word[i] = 'i';
      }
      if (word[i] == u'ù') {
        word[i] = 'u';
      }
      if (word[i] == u'û') {
        word[i] = 'u';
      }
      if (word[i] == u'ô') {
        word[i] = 'o';
      }
      if (word[i] == u'ç') {
        word[i] = 'c';
      }
    }
}

static void removeDoublet(wstring& word) {
    int len = word.size() - 1;
    int i, position;
    wchar_t currentChar;

    if (len > 3) {
      currentChar = word[0];
      position = 1;
      while (word[position]) {
        if (currentChar == word[position]) {
          i = position-1;
          while (word[i] != '\0') {
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


wstring normfrenchword(wstring& word) {
    int len = word.size() - 1;

    if (len > 3) {
      removeAllFEAccent(word);
      removeDoublet(word);
      len = word.size() - 1;
    }

    if (len > 3) {
      if (word[len]==u'e' && word[len-1]==u'i') {
        word.erase(len-1);
        len = len -2;
      }
    }
    if (len > 3) {
      if (word[len]==u'r') {
        word.erase(len);
        len--;
      }
      if (word[len]==u'e') {
        word.erase(len);
        len--;
      }
      if (word[len]==u'e') {
        word.erase(len);
        len--;
      }
      if (word[len] == word[len-1])
        word.erase(len);
    }
    return(word);
}

wstring french_stemming_word(wstring& word) {
    int len = word.size() - 1;

    if (len > 4) {
      if (word[len]==u'x') {
        if (word[len-1]==u'u' && word[len-2]==u'a' && word[len-3]!=u'e') {
          word[len-1]=u'l';  /*  chevaux -> cheval  */
        }                 /*  error :  travaux -> traval but not travail  */
      word.erase(len);
      /*  anneaux -> anneau,  neveux -> neveu  */
      len--;               /*  error :  aulx -> aul but not ail (rare)  */
      }
    }                       /*  error :  yeux -> yeu but not oeil (rare)  */
      if (len > 2) {
        if (word[len]==u'x') {
          word.erase(len);
          /*  peaux -> peau,  poux -> pou  */
          len--;               /*  error :  affreux -> affreu */
        }
      }

      if (len > 2 && word[len]==u's') {  /*  remove final --s --> -- */
      word.erase(len);
        len--;
      }

      if (len > 8) {  /* --issement  -->   --ir */
      if (word[len]==u't'   && word[len-1]==u'n' && word[len-2]==u'e' &&
          word[len-3]==u'm' && word[len-4]==u'e' && word[len-5]==u's' &&
          word[len-6]==u's' && word[len-7]==u'i') {
        word[len-6]=u'r';       /* investissement --> investir */
      word.erase(len-5);
      return(normfrenchword(word));
      }
      }

      if (len > 7) {  /* ---issant  -->   ---ir */
      if (word[len]==u't'   && word[len-1]==u'n' && word[len-2]==u'a' &&
          word[len-3]==u's' && word[len-4]==u's' && word[len-5]==u'i') {
        word[len-4]=u'r';     /* assourdissant --> assourdir */
      word.erase(len-3);
      return(normfrenchword(word));
      }
      }

      if (len > 5) {    /* --ement  -->   --e */
      if (word[len]==u't'   && word[len-1]==u'n' && word[len-2]==u'e' &&
          word[len-3]==u'm' && word[len-4]==u'e') {
        word.erase(len-3);
        /* pratiquement --> pratique */
      if (word[len-5]==u'v' && word[len-6]==u'i') {
        word[len-5]=u'f';     /* administrativement --> administratif */
      word.erase(len-4);
      }
      return(normfrenchword(word));
      }
      }

      if (len > 10) {    /* ---ficatrice  -->   --fier */
      if (word[len]==u'e'   && word[len-1]==u'c' && word[len-2]==u'i' &&
          word[len-3]==u'r' && word[len-4]==u't' && word[len-5]==u'a' &&
          word[len-6]==u'c' && word[len-7]==u'i' && word[len-8]==u'f') {
        word[len-6]=u'e';
        word[len-5]=u'r';
        word.erase(len-4);
        /* justificatrice --> justifier */
      return(normfrenchword(word));
      }
      }

      if (len > 9) {    /* ---ficateur -->   --fier */
      if (word[len]==u'r'   && word[len-1]==u'u' && word[len-2]==u'e' &&
          word[len-3]==u't' && word[len-4]==u'a' && word[len-5]==u'c' &&
          word[len-6]==u'i' && word[len-7]==u'f') {
        word[len-5]=u'e';
        word[len-4]=u'r';
        word.erase(len-3);
        /* justificateur --> justifier */
      return(normfrenchword(word));
      }
      }

      if (len > 8) {    /* ---catrice  -->   --quer */
      if (word[len]==u'e'   && word[len-1]==u'c' && word[len-2]==u'i' &&
          word[len-3]==u'r' && word[len-4]==u't' && word[len-5]==u'a' &&
          word[len-6]==u'c') {
        word[len-6]=u'q';
        word[len-5]=u'u';
        word[len-4]=u'e';
        word[len-3]=u'r';
        word.erase(len-2);
        /* educatrice--> eduquer */
      return(normfrenchword(word));
      }
      }

      if (len > 7) {    /* ---cateur -->   --quer */
      if (word[len]==u'r'   && word[len-1]==u'u' && word[len-2]==u'e' &&
          word[len-3]==u't' && word[len-4]==u'a' && word[len-5]==u'c') {
        word[len-5]=u'q';
        word[len-4]=u'u';
        word[len-3]=u'e';
        word[len-2]=u'r';
        word.erase(len-1);
        /* communicateur--> communiquer */
      return(normfrenchword(word));
      }
      }

      if (len > 7) {    /* ---atrice  -->   --er */
      if (word[len]==u'e'   && word[len-1]==u'c' && word[len-2]==u'i' &&
          word[len-3]==u'r' && word[len-4]==u't' && word[len-5]==u'a') {
        word[len-5]=u'e';
        word[len-4]=u'r';
        word.erase(len-3);
        /* accompagnatrice--> accompagner */
      return(normfrenchword(word));
      }
      }

      if (len > 6) {    /* ---ateur  -->   --er */
      if (word[len]==u'r'   && word[len-1]==u'u' && word[len-2]==u'e' &&
          word[len-3]==u't' && word[len-4]==u'a') {
        word[len-4]=u'e';
        word[len-3]=u'r';
        word.erase(len-2);
        /* administrateur--> administrer */
      return(normfrenchword(word));
      }
      }

      if (len > 5) {    /* --trice  -->   --teur */
      if (word[len]==u'e'   && word[len-1]==u'c' && word[len-2]==u'i' &&
          word[len-3]==u'r' && word[len-4]==u't') {
        word[len-3]=u'e';
        word[len-2]=u'u';
        word[len-1]=u'r';  /* productrice --> producteur */
      word.erase(len);
      /* matrice --> mateur ? */
      len--;
      }
      }

      if (len > 4) {    /* --ième  -->   -- */
      if (word[len]==u'e' && word[len-1]==u'm' && word[len-2]==u'è' &&
          word[len-3]==u'i') {
        word.erase(len-3);
        return(normfrenchword(word));
      }
      }

      if (len > 6) {    /* ---teuse  -->   ---ter */
      if (word[len]==u'e'   && word[len-1]==u's' && word[len-2]==u'u' &&
          word[len-3]==u'e' && word[len-4]==u't') {
        word[len-2]=u'r';
        word.erase(len-1);
        /* acheteuse --> acheter */
      return(normfrenchword(word));
      }
      }

      if (len > 5) {    /* ---teur  -->   ---ter */
      if (word[len]==u'r'   && word[len-1]==u'u' && word[len-2]==u'e' &&
          word[len-3]==u't') {
        word[len-1]=u'r';
        word.erase(len);
        /* planteur --> planter */
      return(normfrenchword(word));
      }
      }

      if (len > 4) {    /* --euse  -->   --eu- */
      if (word[len]==u'e' && word[len-1]==u's' && word[len-2]==u'u' &&
          word[len-3]==u'e') {
        word.erase(len-1);
        /* poreuse --> poreu-,  plieuse --> plieu- */
      return(normfrenchword(word));
      }
      }

      if (len > 7) {    /* ------ère  -->   ------er */
      if (word[len]==u'e' && word[len-1]==u'r' && word[len-2]==u'è') {
        word[len-2]=u'e';
        word[len-1]=u'r';
        word.erase(len);
        /* bijoutière --> bijoutier,  caissière -> caissier */
      return(normfrenchword(word));
      }
      }

      if (len > 6) {    /* -----ive  -->   -----if */
      if (word[len]==u'e' && word[len-1]==u'v' && word[len-2]==u'i') {
        word[len-1]=u'f';   /* but not convive */
        word.erase(len);
      /* abrasive --> abrasif */
      return(normfrenchword(word));
      }
      }

      if (len > 3) {    /* folle or molle  -->   fou or mou */
      if (word[len]==u'e' && word[len-1]==u'l' && word[len-2]==u'l' &&
          word[len-3]==u'o' && (word[len-4]==u'f' || word[len-4]==u'm')) {
        word[len-2]=u'u';
        word.erase(len-1);
        /* folle --> fou */
      return(normfrenchword(word));
      }
      }

      if (len > 8) {    /* ----nnelle  -->   ----n */
        if (word[len]==u'e'   && word[len-1]==u'l' && word[len-2]==u'l' &&
            word[len-3]==u'e' && word[len-4]==u'n' && word[len-5]==u'n') {
          word.erase(len-4);
          /* personnelle --> person */
          return(normfrenchword(word));
        }
      }

      if (len > 8) {    /* ----nnel  -->   ----n */
        if (word[len]==u'l'   && word[len-1]==u'e' && word[len-2]==u'n' &&
            word[len-3]==u'n') {
          word.erase(len-2);
          /* personnel --> person */
          return(normfrenchword(word));
        }
      }

      if (len > 3) {    /* --ète  -->  et */
        if (word[len]==u'e' && word[len-1]==u't' && word[len-2]==u'è') {
          word[len-2]=u'e';
          word.erase(len);
          /* complète --> complet */
          len--;
        }
      }

      if (len > 7) {    /* -----ique  -->   */
        if (word[len]==u'e' && word[len-1]==u'u' && word[len-2]==u'q' &&
            word[len-3]==u'i') {
          word.erase(len-3);
          /* aromatique --> aromat */
          len = len-4;
        }
      }

      if (len > 7) {    /* -----esse -->    */
        if (word[len]==u'e' && word[len-1]==u's' && word[len-2]==u's' &&
            word[len-3]==u'e') {
          word.erase(len-2);
          /* faiblesse --> faible */
          return(normfrenchword(word));
        }
      }

      if (len > 6) {    /* ---inage -->    */
        if (word[len]==u'e' && word[len-1]==u'g' && word[len-2]==u'a' &&
            word[len-3]==u'n' && word[len-4]==u'i') {
          word.erase(len-2);
          /* patinage --> patin */
          return(normfrenchword(word));
        }
      }

      if (len > 8) {    /* ---isation -->   - */
        if (word[len]==u'n'   && word[len-1]==u'o' && word[len-2]==u'i' &&
            word[len-3]==u't' && word[len-4]==u'a' && word[len-5]==u's' &&
            word[len-6]==u'i') {
          word.erase(len-6);
          /* sonorisation --> sonor */
          if (len > 11 && word[len-7]==u'l' && word[len-8]==u'a' && word[len-9]==u'u')
            word[len-8]=u'e';  /* ritualisation --> rituel */
          return(normfrenchword(word));
        }
      }

      if (len > 8) {    /* ---isateur -->   - */
        if (word[len]==u'r'   && word[len-1]==u'u' && word[len-2]==u'e' && word[len-3]==u't' &&
            word[len-4]==u'a' && word[len-5]==u's' && word[len-6]==u'i') {
          word.erase(len-6);
          /* colonisateur --> colon */
          return(normfrenchword(word));
        }
      }

      if (len > 7) {    /* ----ation -->   - */
        if (word[len]==u'n'   && word[len-1]==u'o' && word[len-2]==u'i' &&
            word[len-3]==u't' && word[len-4]==u'a') {
          word.erase(len-4);
          /* nomination --> nomin */
          return(normfrenchword(word));
        }
      }

      if (len > 7) {    /* ----ition -->   - */
        if (word[len]==u'n'   && word[len-1]==u'o' && word[len-2]==u'i' &&
            word[len-3]==u't' && word[len-4]==u'i') {
          word.erase(len-4);
          /* disposition --> dispos */
          return(normfrenchword(word));
        }
      }

      /* various other suffix */
      return(normfrenchword(word));
}

//' Stem French words
//'
//' Stemmer for French words
//'
//' @param words a [character] containing the original words.
//' @return [character] with stemmed words.
//' @examples
//' french_stemmer(c("tester", "testament", "clients"))
//' @export
// [[Rcpp::export]]
CharacterVector french_stemmer(Rcpp::StringVector words) {
  CharacterVector result(words.size());

  for (int i = 0; i < words.size(); ++i) {
    string s1 = static_cast<string>(words[i]);
    wstring str2 = utf8_to_utf16(s1);
    result[i] = french_stemming_word(str2);
    Rcpp::checkUserInterrupt();
  }

  return result;
}
