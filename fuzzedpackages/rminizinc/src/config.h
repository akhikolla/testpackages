#define PKG_MZN_PARSE 
#define PKG_MZN_EVAL 

#if defined(PKG_MZN_PARSE) && !(~(~PKG_MZN_PARSE + 0) == 0 && ~(~PKG_MZN_PARSE + 1) == 1) 
  #define MZN_PARSE 30
#elif defined(PKG_MZN_EVAL) && !(~(~PKG_MZN_EVAL + 0) == 0 && ~(~PKG_MZN_EVAL + 1) == 1)
    #define MZN_EVAL 30
#else
    #define MZN_ABSENT 10
#endif

