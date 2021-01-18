// Copyright (c) 2019-2020 ISciences, LLC.
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

#ifndef EXACTEXTRACT_OPERATION_H
#define EXACTEXTRACT_OPERATION_H

#include <functional>
#include <string>

#include "grid.h"
#include "gdal_raster_wrapper.h"
#include "raster_stats.h"

namespace exactextract {

    class Operation {
    public:
        Operation(std::string p_stat, std::string p_name, RasterSource* p_values, RasterSource* p_weights = nullptr) :
                stat{std::move(p_stat)},
                name{std::move(p_name)},
                values{p_values},
                weights{p_weights} {}

        bool weighted() const {
            return weights != nullptr;
        }

        Grid<bounded_extent> grid() const {
            if (weighted()) {
                return values->grid().common_grid(weights->grid());
            } else {
                return values->grid();
            }
        }

        std::function<nonstd::optional<double>(const RasterStats<double>&)> result_fetcher() const {
            if (stat == "mean") {
                return [](const RasterStats<double> & s) { return s.mean(); };
            } else if (stat == "sum") {
                return [](const RasterStats<double> & s) { return s.sum(); };
            } else if (stat == "count") {
                return [](const RasterStats<double> & s) { return s.count(); };
            } else if (stat == "weighted_mean") {
                return [](const RasterStats<double> & s) { return s.weighted_mean(); };
            } else if (stat == "weighted_sum") {
                return [](const RasterStats<double> & s) { return s.weighted_sum(); };
            } else if (stat == "min") {
                return [](const RasterStats<double> & s) { return s.min(); };
            } else if (stat == "max") {
                return [](const RasterStats<double> & s) { return s.max(); };
            } else if (stat == "majority" || stat == "mode") {
                return [](const RasterStats<double> & s) { return s.mode(); };
            } else if (stat == "minority") {
                return [](const RasterStats<double> & s) { return s.minority(); };
            } else if (stat == "variety") {
                return [](const RasterStats<double> & s) { return s.variety(); };
            } else if (stat == "stdev") {
                return [](const RasterStats<double> & s) { return s.stdev(); };
            } else if (stat == "variance") {
                return [](const RasterStats<double> & s) { return s.variance(); };
            } else if (stat == "coefficient_of_variation") {
                return [](const RasterStats<double> & s) { return s.coefficient_of_variation(); };
            } else {
                throw std::runtime_error("Unknown stat: '" + stat + "'");
            }
        }

        std::string stat;
        std::string name;
        RasterSource* values;
        RasterSource* weights;
    };

}

#endif //EXACTEXTRACT_OPERATION_H
