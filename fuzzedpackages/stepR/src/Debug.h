#ifndef DEBUG

// globally turn debugging on or off
// changes require recompilation: rm stepR/src/*.o
#define DEBUG
#undef DEBUG // comment out to turn debugging on
// debug specific functions
#ifdef DEBUG
#define DEBUGbounded // for Step::bounded
#endif

#endif
