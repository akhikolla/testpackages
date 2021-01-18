// Copyright (c) 2018-2019 ISciences, LLC.
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

#include "gdal_dataset_wrapper.h"

#include <algorithm>
#include <memory>

namespace exactextract {

    GDALDatasetWrapper::GDALDatasetWrapper(const std::string & filename, const std::string & layer, std::string id_field) :
    m_id_field{std::move(id_field)}
    {
        m_dataset = GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr);
        if (m_dataset == nullptr) {
            throw std::runtime_error("Failed to open " + filename);
        }

        bool numeric = std::all_of(layer.begin(),layer.end(),
                                   [](char c) { return std::isdigit(c); });

        if (numeric) {
            m_layer = GDALDatasetGetLayer(m_dataset, std::stoi(layer));
        } else {
            m_layer = GDALDatasetGetLayerByName(m_dataset, layer.c_str());
        }

        if (m_layer == nullptr) {
            throw std::runtime_error("No layer " + layer + " found in " + filename);
        }

        OGR_L_ResetReading(m_layer);
        m_feature = nullptr;

        auto defn = OGR_L_GetLayerDefn(m_layer);
        auto index = OGR_FD_GetFieldIndex(defn, m_id_field.c_str());

        if (index == -1) {
            throw std::runtime_error("ID field '" + m_id_field + "' not found in " + filename + ".");
        }
    }

    bool GDALDatasetWrapper::next() {
        if (m_feature != nullptr) {
            OGR_F_Destroy(m_feature);
        }
        m_feature = OGR_L_GetNextFeature(m_layer);
        return m_feature != nullptr;
    }

    GEOSGeometry* GDALDatasetWrapper::feature_geometry(const GEOSContextHandle_t &geos_context) const {
        OGRGeometryH geom = OGR_F_GetGeometryRef(m_feature);

        auto sz = static_cast<size_t>(OGR_G_WkbSize(geom));
        auto buff = std::make_unique<unsigned char[]>(sz);
        OGR_G_ExportToWkb(geom, wkbXDR, buff.get());

        return GEOSGeomFromWKB_buf_r(geos_context, buff.get(), sz);
    }

    std::string GDALDatasetWrapper::feature_field(const std::string &field_name) const {
        int index = OGR_F_GetFieldIndex(m_feature, field_name.c_str());
        // TODO check handling of invalid field name
        return OGR_F_GetFieldAsString(m_feature, index);
    }

    void GDALDatasetWrapper::copy_field(const std::string & name, OGRLayerH copy_to) const {
        auto src_layer_defn = OGR_L_GetLayerDefn(m_layer);
        auto src_index = OGR_FD_GetFieldIndex(src_layer_defn, name.c_str());

        if (src_index == -1) {
            throw std::runtime_error("Cannot find field " + name);
        }

        auto src_field_defn = OGR_FD_GetFieldDefn(src_layer_defn, src_index);

        OGR_L_CreateField(copy_to, src_field_defn, true);
    }

    GDALDatasetWrapper::~GDALDatasetWrapper(){
        GDALClose(m_dataset);

        if (m_feature != nullptr) {
            OGR_F_Destroy(m_feature);
        }
    }
}
