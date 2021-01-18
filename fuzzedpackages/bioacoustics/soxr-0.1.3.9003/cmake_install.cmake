# Install script for directory: /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003

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
   "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/share/doc/libsoxr/README;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/share/doc/libsoxr/LICENCE;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/share/doc/libsoxr/NEWS")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/share/doc/libsoxr" TYPE FILE FILES
    "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/README"
    "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/LICENCE"
    "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/NEWS"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/src/cmake_install.cmake")
  include("/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/cmake_install.cmake")
  include("/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/examples/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
