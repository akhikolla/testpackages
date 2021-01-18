/*
 * Profiler.h
 *
 *  Created on: Sep 9, 2016
 *      Author: mabaker
 */

#ifndef SRC_PROFILER_H_
#define SRC_PROFILER_H_

#include "mem_rss.h"
#include "CumulativeTimeMeasurement.h"

#include <stddef.h> // C's size_t

#include <string>

namespace SignificantPattern
{

    class Profiler
    {
    // all measurements are _cumulative_; use reset to start over
    private:
        CumulativeTimeMeasurement executionTime;
        CumulativeTimeMeasurement initialisationTime;
        CumulativeTimeMeasurement postprocessingAndCleanupTime;
        CumulativeTimeMeasurement fileIOTime;
        CumulativeTimeMeasurement significanceThresholdComputeTime;
        CumulativeTimeMeasurement significantIntervalsComputeTime;

        size_t peak_memory;

    protected:
        double measureTime();

    public:
        Profiler ();
        virtual ~Profiler ();

        /// total execution time in secs
        inline double getExecutionTime() const {
            return executionTime.getTime();
        }
        inline void markStartExecution() { executionTime.markStart(); };
        inline void markEndExecution() { executionTime.markEnd(); };

        /// initialisation time in secs
        inline double getInitialisationTime() const {
            return initialisationTime.getTime();
        }
        inline void markStartInitialisation() { initialisationTime.markStart(); };
        inline void markEndInitialisation() { initialisationTime.markEnd(); };

        /// file I/O time in secs
        inline double getFileIOTime() const { return fileIOTime.getTime(); }
        inline void markStartFileIO() { fileIOTime.markStart(); };
        inline void markEndFileIO() { fileIOTime.markEnd(); };

        /// time to compute corrected significance threshold secs
        inline double getSignificanceThresholdComputeTime() const {
            return significanceThresholdComputeTime.getTime();
        }
        inline void markStartSignificantIntervalsCompute() { significantIntervalsComputeTime.markStart(); };
        inline void markEndSignificantIntervalsCompute() { significantIntervalsComputeTime.markEnd(); };

        // time to find significant intervals
        inline double getSignificantIntervalsComputeTime() const {
            return significantIntervalsComputeTime.getTime();
        }
        inline void markStartSignificanceThresholdCompute() { significanceThresholdComputeTime.markStart(); };
        inline void markEndSignificanceThresholdCompute() { significanceThresholdComputeTime.markEnd(); };

        /// post-processing and cleanup time in secs
        inline double getPostprocessingAndCleanupTimeTime() const {
            return postprocessingAndCleanupTime.getTime();
        }
        inline void markStartPostprocessingAndCleanup() { postprocessingAndCleanupTime.markStart(); };
        inline void markEndPostprocessingAndCleanup() { postprocessingAndCleanupTime.markEnd(); };


        // current and peak memory (0 if not marked)
        inline const size_t getCurrMemory() const { return getCurrentRSS(); }
        inline const size_t getPeakMemory() const { return peak_memory; }
        void markPeakMemory();


        void reset(bool resetPeakMemory = false, bool resetIO = false);


        void writeToFile(const std::string& filename) const;

    };

} /* namespace SignificantPattern */

#endif /* SRC_PROFILER_H_ */
