#ifndef LOWPASSFILTER_H_FILTER
#define LOWPASSFILTER_H_FILTER

/***************
* class Filter
* abstract class maintaining the filter
* Florian Pein, 2016
***************/
class Filter {
  public:
    virtual ~Filter();
    
    // returns the antiderivative at point t
    virtual double antiderivative(const double &t) const = 0;
};

#endif
