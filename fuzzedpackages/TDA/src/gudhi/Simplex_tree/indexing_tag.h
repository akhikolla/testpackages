/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SIMPLEX_TREE_INDEXING_TAG_H_
#define SIMPLEX_TREE_INDEXING_TAG_H_

namespace Gudhi {

/** \brief Tag for a linear ordering of simplices. 
 *
 * \implements IndexingTag
 */
struct linear_indexing_tag {
};

/* \brief Tag for a zigzag ordering of simplices. */
//  struct zigzag_indexing_tag {};
}  // namespace Gudhi

#endif  // SIMPLEX_TREE_INDEXING_TAG_H_
