/**
 * @file IsSorted.h
 *
 * This file is part of Cyclops
 *
 * Copyright 2020 Observational Health Data Sciences and Informatics
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef __IsSorted_h__
#define __IsSorted_h__

#include <Rcpp.h>
#include <vector>
#include <cstdint>

using namespace Rcpp ;

namespace ohdsi {
	namespace cyclops {

		struct IsSorted {
		public:
			static bool isSorted(const DataFrame& dataFrame,const std::vector<std::string>& indexes,const std::vector<bool>& ascending);
            static bool isSorted(const List& vectorList, const std::vector<bool>& ascending);
		};
	}
}

#endif // __IsSorted_h__
