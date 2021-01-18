// Copyright (c) 2020 ISciences, LLC.
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

#ifndef EXACTEXTRACT_RASTER_SOURCE_H
#define EXACTEXTRACT_RASTER_SOURCE_H

#include "box.h"
#include "grid.h"
#include "raster.h"

namespace exactextract {
    class RasterSource {
    public:

        virtual const Grid<bounded_extent> &grid() const = 0;
        virtual std::unique_ptr<AbstractRaster<double>> read_box(const Box &box) = 0;

        virtual ~RasterSource() = default;

        void set_name(const std::string & name) {
            m_name = name;
        }

        std::string name() const {
            return m_name;
        }

    private:
        std::string m_name;
    };
}

#endif //EXACTEXTRACT_RASTER_SOURCE_H
