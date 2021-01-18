// Copyright (c) 2019 ISciences, LLC.
// All rights reserved.
//
// This software is licensed under the Apache License, Version 2.0 (the "License").
// You may not use this file except in compliance with the License. You may
// obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef EXACTEXTRACT_PROCESSOR_H
#define EXACTEXTRACT_PROCESSOR_H

#include <iostream>
#include <string>

#include "gdal_dataset_wrapper.h"
#include "output_writer.h"
#include "stats_registry.h"


namespace exactextract {
    class Processor {

    public:
        // FIXME add GEOS error/notice handlers
        Processor(GDALDatasetWrapper & ds, OutputWriter & out, const std::vector<Operation> & ops) :
                m_reg{},
                m_geos_context{initGEOS_r(nullptr, nullptr)},
                m_output{out},
                m_shp{ds},
                m_operations{ops}
        {
            m_output.set_registry(&m_reg);
        }

        virtual ~Processor() {
            finishGEOS_r(m_geos_context);
        }

        virtual void process()= 0;

        void set_max_cells_in_memory(size_t n) {
            m_max_cells_in_memory = n;
        }

        void show_progress(bool val) {
            m_show_progress = val;
        }

    protected:

        template<typename T>
        void progress(const T & name) const {
            if (m_show_progress)
                std::cout << std::endl << "Processing " << name << std::flush;
        }

        void progress() const {
            if (m_show_progress)
                std::cout << "." << std::flush;
        }

        StatsRegistry m_reg;

        GEOSContextHandle_t m_geos_context;

        OutputWriter& m_output;

        GDALDatasetWrapper& m_shp;

        bool store_values=false;
        bool m_show_progress=false;

        std::vector<Operation> m_operations;

        size_t m_max_cells_in_memory = 1000000L;
    };
}

#endif //EXACTEXTRACT_PROCESSOR_H
