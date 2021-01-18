#if defined(__GNUC__) && __GNUC__ > 3
/* Define to 1 if the system has the `__builtin_popcount' built-in function */
#define HAVE___BUILTIN_POPCOUNT 1

/* Define to 1 if the system has the `__builtin_popcountl' built-in function
   */
#define HAVE___BUILTIN_POPCOUNTL 1

/* Define to 1 if the system has the `__builtin_popcountll' built-in function
   */
#define HAVE___BUILTIN_POPCOUNTLL 1
#else // old
/* Define to 1 if the system has the `__builtin_popcount' built-in function */
#define HAVE___BUILTIN_POPCOUNT 0

/* Define to 1 if the system has the `__builtin_popcountl' built-in function
   */
#define HAVE___BUILTIN_POPCOUNTL 0

/* Define to 1 if the system has the `__builtin_popcountll' built-in function
   */
#define HAVE___BUILTIN_POPCOUNTLL 0
#endif
