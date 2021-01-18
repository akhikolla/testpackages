# Install script for directory: /Users/hadley/Documents/ingest/feather/cpp/src/feather

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/hadley/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "DEBUG")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/feather" TYPE FILE FILES
    "/Users/hadley/Documents/ingest/feather/cpp/src/feather/api.h"
    "/Users/hadley/Documents/ingest/feather/cpp/src/feather/buffer.h"
    "/Users/hadley/Documents/ingest/feather/cpp/src/feather/common.h"
    "/Users/hadley/Documents/ingest/feather/cpp/src/feather/io.h"
    "/Users/hadley/Documents/ingest/feather/cpp/src/feather/metadata.h"
    "/Users/hadley/Documents/ingest/feather/cpp/src/feather/reader.h"
    "/Users/hadley/Documents/ingest/feather/cpp/src/feather/status.h"
    "/Users/hadley/Documents/ingest/feather/cpp/src/feather/types.h"
    "/Users/hadley/Documents/ingest/feather/cpp/src/feather/writer.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/hadley/Documents/ingest/feather/cpp/build/debug/libfeather.dylib")
endif()

