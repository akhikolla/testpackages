# CMake generated Testfile for 
# Source directory: /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests
# Build directory: /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(20-bit-perfect-44100-192000 "/usr/bin/cmake" "-Dbits=20" "-DBIN=./" "-DEXAMPLES_BIN=../examples/" "-DlenToSkip=1" "-Dorate=192000" "-Dirate=44100" "-Dlen=16" "-P" "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/cmp-test.cmake")
set_tests_properties(20-bit-perfect-44100-192000 PROPERTIES  _BACKTRACE_TRIPLES "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;32;add_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;53;add_cmp_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;0;")
add_test(20-bit-perfect-192000-44100 "/usr/bin/cmake" "-Dbits=20" "-DBIN=./" "-DEXAMPLES_BIN=../examples/" "-DlenToSkip=1" "-Dorate=44100" "-Dirate=192000" "-Dlen=16" "-P" "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/cmp-test.cmake")
set_tests_properties(20-bit-perfect-192000-44100 PROPERTIES  _BACKTRACE_TRIPLES "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;32;add_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;54;add_cmp_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;0;")
add_test(20-bit-perfect-44100-65537 "/usr/bin/cmake" "-Dbits=20" "-DBIN=./" "-DEXAMPLES_BIN=../examples/" "-DlenToSkip=1" "-Dorate=65537" "-Dirate=44100" "-Dlen=16" "-P" "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/cmp-test.cmake")
set_tests_properties(20-bit-perfect-44100-65537 PROPERTIES  _BACKTRACE_TRIPLES "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;32;add_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;53;add_cmp_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;0;")
add_test(20-bit-perfect-65537-44100 "/usr/bin/cmake" "-Dbits=20" "-DBIN=./" "-DEXAMPLES_BIN=../examples/" "-DlenToSkip=1" "-Dorate=44100" "-Dirate=65537" "-Dlen=16" "-P" "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/cmp-test.cmake")
set_tests_properties(20-bit-perfect-65537-44100 PROPERTIES  _BACKTRACE_TRIPLES "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;32;add_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;54;add_cmp_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;0;")
add_test(28-bit-perfect-44100-192000 "/usr/bin/cmake" "-Dbits=28" "-DBIN=./" "-DEXAMPLES_BIN=../examples/" "-DlenToSkip=1" "-Dorate=192000" "-Dirate=44100" "-Dlen=16" "-P" "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/cmp-test.cmake")
set_tests_properties(28-bit-perfect-44100-192000 PROPERTIES  _BACKTRACE_TRIPLES "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;32;add_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;53;add_cmp_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;0;")
add_test(28-bit-perfect-192000-44100 "/usr/bin/cmake" "-Dbits=28" "-DBIN=./" "-DEXAMPLES_BIN=../examples/" "-DlenToSkip=1" "-Dorate=44100" "-Dirate=192000" "-Dlen=16" "-P" "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/cmp-test.cmake")
set_tests_properties(28-bit-perfect-192000-44100 PROPERTIES  _BACKTRACE_TRIPLES "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;32;add_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;54;add_cmp_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;0;")
add_test(28-bit-perfect-44100-65537 "/usr/bin/cmake" "-Dbits=28" "-DBIN=./" "-DEXAMPLES_BIN=../examples/" "-DlenToSkip=1" "-Dorate=65537" "-Dirate=44100" "-Dlen=16" "-P" "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/cmp-test.cmake")
set_tests_properties(28-bit-perfect-44100-65537 PROPERTIES  _BACKTRACE_TRIPLES "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;32;add_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;53;add_cmp_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;0;")
add_test(28-bit-perfect-65537-44100 "/usr/bin/cmake" "-Dbits=28" "-DBIN=./" "-DEXAMPLES_BIN=../examples/" "-DlenToSkip=1" "-Dorate=44100" "-Dirate=65537" "-Dlen=16" "-P" "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/cmp-test.cmake")
set_tests_properties(28-bit-perfect-65537-44100 PROPERTIES  _BACKTRACE_TRIPLES "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;32;add_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;54;add_cmp_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;0;")
add_test(1-delay-clear "./1-delay-clear")
set_tests_properties(1-delay-clear PROPERTIES  _BACKTRACE_TRIPLES "/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;62;add_test;/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeLists.txt;0;")
