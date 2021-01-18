/*
 * Profiler.cpp
 *
 *  Created on: Sep 9, 2016
 *      Author: mabaker
 */

#include "Profiler.h"
#include "Exception.h"
#include "mem_rss.h"

#include <fstream>
#include <stddef.h> // C's size_t


namespace SignificantPattern
{

Profiler::Profiler (): peak_memory(0)
{
}
Profiler::~Profiler ()
{
}

void Profiler::markPeakMemory()
{
//     peak_memory = std::max(peak_memory, getPeakRSS());
//     NOTE: Needed to remove getPeakRSS
    peak_memory = 0;
}


void Profiler::reset(bool resetPeakMemory, bool resetIO) {
    executionTime.reset();
    initialisationTime.reset();
    significanceThresholdComputeTime.reset();
    significantIntervalsComputeTime.reset();
    postprocessingAndCleanupTime.reset();
    if (resetIO) fileIOTime.reset();
    if (resetPeakMemory) peak_memory=0;
}

void Profiler::writeToFile (const std::string& filename) const
{
    std::ofstream timing_file;
    timing_file.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
    try
    {
            timing_file.open(filename.c_str());
    }
    catch (const std::ios_base::failure& e)
    {
        throw Exception("Failed opening " + filename);
    }

    timing_file << "CODE PROFILING" << std::endl;
    timing_file << "Total Execution time: " << getExecutionTime() << " (s)." << std::endl;
    timing_file << "\tInitialisation time: " << getInitialisationTime() << " (s)." << std::endl;
    timing_file << "\tTime to compute corrected significance threshold: " << getSignificanceThresholdComputeTime() << " (s)." << std::endl;
    timing_file << "\tTime to find significant intervals: " << getSignificantIntervalsComputeTime() << " (s)." << std::endl;
    timing_file << "\tPost-processing and cleanup time: " << getPostprocessingAndCleanupTimeTime() << " (s)." << std::endl;
    timing_file << "File I/O time: " << getFileIOTime() << " (s)." << std::endl;

    timing_file << "Peak memory usage: " << ((size_t) getPeakMemory()/1024) << " (KB)" << std::endl;

    timing_file.close();
}

} /* namespace SignificantPattern */

