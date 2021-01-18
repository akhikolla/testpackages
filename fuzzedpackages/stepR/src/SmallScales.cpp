#include "SmallScales.h"
# include <algorithm>

/***************
* class small scales
* maintains additional finding on small scales
* Florian Pein, 2017
***************/
std::list<SmallScales> SmallScales::listSmallScales_;
std::list<SmallScales>::iterator SmallScales::it_;

SmallScales::SmallScales() {}

SmallScales::SmallScales(unsigned int left, unsigned int right, unsigned int li, unsigned int ri,
                         double stat, bool noDeconvolution) :
  left_(left), right_(right), li_(li), ri_(ri), stat_(stat), noDeconvolution_(noDeconvolution) {}

unsigned int SmallScales::left() {
  return left_;
}

unsigned int SmallScales::right() {
  return right_;
}

unsigned int SmallScales::li() {
  return li_;
}

unsigned int SmallScales::ri() {
  return ri_;
}

double SmallScales::stat() {
  return stat_;
}

bool SmallScales::noDeconvolution() {
  return noDeconvolution_;
}

void SmallScales::update(unsigned int start, unsigned int len, double stat) {
  // find first segment that touches, intersets or is right of the new one
  while (it_ != listSmallScales_.end() && it_ -> ri() < start) {
    ++it_;
  }
  
  SmallScales newSegment = SmallScales(start + 1L, start + len + 1L, start + 1L, start + len + 1L,
                                       stat, false);
  
  // consider all segments for which li() <= start + len + 2L (all other are right of it with a gap)
  std::list<SmallScales>::iterator it(it_);
  unsigned int number = 0u;
  bool replace = true;
  while (it != listSmallScales_.end() && it -> li() <= start + len + 2L) {
    ++number;
    newSegment.extend(it -> li(), it -> ri());
    if (it -> stat() >= stat) {
      replace = false;
    }
    
    ++it;
  }
  
  if (number == 0u) {
    listSmallScales_.insert(it_, newSegment);
    --it_;
  } else {
    if (replace) {
      if (number > 1u) {
        it_ -> replace(start, len, newSegment.li(), newSegment.ri(), stat, true);
        it = it_;
        ++it;
        while (it != listSmallScales_.end() && it -> left() <= start + len + 2L) {
          it = listSmallScales_.erase(it);
        }
      } else {
        it_ -> replace(start, len, newSegment.li(), newSegment.ri(), stat, false);
      }
    } else {
      it = it_;
      while (it != listSmallScales_.end() && it -> li() <= start + len + 2L) {
        it -> extend(start + 1L, start + len + 1L);
        ++ it;
      }
    }
  }
}

void SmallScales::replace(unsigned int start, unsigned int len, unsigned int li, unsigned int ri,
                          double stat, bool noDe) {
  left_ = start + 1;
  right_ = start + len + 1;
  li_ = li;
  ri_ = ri;
  stat_ = stat;
  if (noDe) {
    noDeconvolution_ = true;
  }
}

void SmallScales::extend(unsigned int li, unsigned int ri) {
  li_ = std::min(li_, li);
  ri_ = std::max(ri_, ri);
}

void SmallScales::cleanUpGlobalVariables() {
  std::list<SmallScales> tmpListSmallScales;
  listSmallScales_.swap(tmpListSmallScales);
}
