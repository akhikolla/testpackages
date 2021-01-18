#ifndef CONSTANTS_H
#define CONSTANTS_H

inline constexpr size_t cacheline_size(){
  return 128L;
}

#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#define PSQN_RESTRICT
#else
#define PSQN_RESTRICT __restrict__
#endif

#endif
