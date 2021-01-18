# Install script for directory: /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/examples

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
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

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/share/doc/libsoxr/examples/1-single-block.c;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/share/doc/libsoxr/examples/1a-lsr.c;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/share/doc/libsoxr/examples/2-stream.C;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/share/doc/libsoxr/examples/3-options-input-fn.c;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/share/doc/libsoxr/examples/4-split-channels.c;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/share/doc/libsoxr/examples/5-variable-rate.c;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/share/doc/libsoxr/examples/examples-common.h;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/share/doc/libsoxr/examples/README")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/share/doc/libsoxr/examples" TYPE FILE FILES
    "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/examples/1-single-block.c"
    "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/examples/1a-lsr.c"
    "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/examples/2-stream.C"
    "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/examples/3-options-input-fn.c"
    "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/examples/4-split-channels.c"
    "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/examples/5-variable-rate.c"
    "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/examples/examples-common.h"
    "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/examples/README"
    )
endif()

