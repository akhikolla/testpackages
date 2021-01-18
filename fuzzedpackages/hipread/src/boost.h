#ifndef HIPREAD_BOOST_H_
#define HIPREAD_BOOST_H_

// Adapted from readr boost.h
// (I belive the purpose is to have a header-only file so that
// we can have the '#pragma' line, which hides warnings from the
// boost library - before doing it this way, there were a bunch
// about how the shared_ptr is deprecated)

#pragma GCC system_header

#include <boost/shared_ptr.hpp>
#include <boost/spirit/include/qi.hpp>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

#endif
