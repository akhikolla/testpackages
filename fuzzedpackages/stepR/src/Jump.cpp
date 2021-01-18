
#include "Jump.h"

/***************
* class Jump
* describes a jump in a block signal
* Thomas Hotz, 2007
***************/

Jump::Jump() : number(NA_INTEGER), rightIndex(-1), improve(-1) {}

Jump::Jump(int n, int ri, double im) : number(n), rightIndex(ri), improve(im) {}
