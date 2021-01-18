/*
 * CumulativeTimeMeasurement.h
 *
 *  Created on: 24 Mar 2017
 *      Author: mikolajr
 */

#ifndef LIBSIGINTERVALSEARCH_CUMULATIVETIMEMEASUREMENT_H_
#define LIBSIGINTERVALSEARCH_CUMULATIVETIMEMEASUREMENT_H_

namespace SignificantPattern {

class CumulativeTimeMeasurement {
private:
    double time = 0;
    double start = 0;

    double measureTime();

public:
    CumulativeTimeMeasurement() = default;
    virtual ~CumulativeTimeMeasurement() = default;

    // Cumulative marked time (in seconds)
    inline double getTime() const { return time; }

    inline void markStart() { start = measureTime(); };
    inline void markEnd() { time += measureTime() - start; start = 0; };

    inline void reset() { time = 0; start = 0; }
};

} /* namespace SignificantPattern */

#endif /* LIBSIGINTERVALSEARCH_CUMULATIVETIMEMEASUREMENT_H_ */
