# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003

# Include any dependencies generated for this target.
include tests/CMakeFiles/1-delay-clear.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/1-delay-clear.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/1-delay-clear.dir/flags.make

tests/CMakeFiles/1-delay-clear.dir/1-delay-clear.c.o: tests/CMakeFiles/1-delay-clear.dir/flags.make
tests/CMakeFiles/1-delay-clear.dir/1-delay-clear.c.o: tests/1-delay-clear.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object tests/CMakeFiles/1-delay-clear.dir/1-delay-clear.c.o"
	cd /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests && /usr/bin/gcc  -std=gnu99 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/1-delay-clear.dir/1-delay-clear.c.o   -c /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/1-delay-clear.c

tests/CMakeFiles/1-delay-clear.dir/1-delay-clear.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/1-delay-clear.dir/1-delay-clear.c.i"
	cd /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests && /usr/bin/gcc  -std=gnu99 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/1-delay-clear.c > CMakeFiles/1-delay-clear.dir/1-delay-clear.c.i

tests/CMakeFiles/1-delay-clear.dir/1-delay-clear.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/1-delay-clear.dir/1-delay-clear.c.s"
	cd /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests && /usr/bin/gcc  -std=gnu99 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/1-delay-clear.c -o CMakeFiles/1-delay-clear.dir/1-delay-clear.c.s

# Object files for target 1-delay-clear
1__delay__clear_OBJECTS = \
"CMakeFiles/1-delay-clear.dir/1-delay-clear.c.o"

# External object files for target 1-delay-clear
1__delay__clear_EXTERNAL_OBJECTS =

tests/1-delay-clear: tests/CMakeFiles/1-delay-clear.dir/1-delay-clear.c.o
tests/1-delay-clear: tests/CMakeFiles/1-delay-clear.dir/build.make
tests/1-delay-clear: src/libsoxr.a
tests/1-delay-clear: tests/CMakeFiles/1-delay-clear.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable 1-delay-clear"
	cd /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/1-delay-clear.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/1-delay-clear.dir/build: tests/1-delay-clear

.PHONY : tests/CMakeFiles/1-delay-clear.dir/build

tests/CMakeFiles/1-delay-clear.dir/clean:
	cd /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests && $(CMAKE_COMMAND) -P CMakeFiles/1-delay-clear.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/1-delay-clear.dir/clean

tests/CMakeFiles/1-delay-clear.dir/depend:
	cd /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003 /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003 /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests /home/akhila/fuzzer_packages/fuzzedpackages/bioacoustics/soxr-0.1.3.9003/tests/CMakeFiles/1-delay-clear.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/1-delay-clear.dir/depend

