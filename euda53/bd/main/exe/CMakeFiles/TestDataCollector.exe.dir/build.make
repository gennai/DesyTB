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
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/gennai/DESY/euda53

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/gennai/DESY/euda53/bd

# Include any dependencies generated for this target.
include main/exe/CMakeFiles/TestDataCollector.exe.dir/depend.make

# Include the progress variables for this target.
include main/exe/CMakeFiles/TestDataCollector.exe.dir/progress.make

# Include the compile flags for this target's objects.
include main/exe/CMakeFiles/TestDataCollector.exe.dir/flags.make

main/exe/CMakeFiles/TestDataCollector.exe.dir/src/TestDataCollector.cxx.o: main/exe/CMakeFiles/TestDataCollector.exe.dir/flags.make
main/exe/CMakeFiles/TestDataCollector.exe.dir/src/TestDataCollector.cxx.o: ../main/exe/src/TestDataCollector.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/gennai/DESY/euda53/bd/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object main/exe/CMakeFiles/TestDataCollector.exe.dir/src/TestDataCollector.cxx.o"
	cd /Users/gennai/DESY/euda53/bd/main/exe && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestDataCollector.exe.dir/src/TestDataCollector.cxx.o -c /Users/gennai/DESY/euda53/main/exe/src/TestDataCollector.cxx

main/exe/CMakeFiles/TestDataCollector.exe.dir/src/TestDataCollector.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestDataCollector.exe.dir/src/TestDataCollector.cxx.i"
	cd /Users/gennai/DESY/euda53/bd/main/exe && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/gennai/DESY/euda53/main/exe/src/TestDataCollector.cxx > CMakeFiles/TestDataCollector.exe.dir/src/TestDataCollector.cxx.i

main/exe/CMakeFiles/TestDataCollector.exe.dir/src/TestDataCollector.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestDataCollector.exe.dir/src/TestDataCollector.cxx.s"
	cd /Users/gennai/DESY/euda53/bd/main/exe && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/gennai/DESY/euda53/main/exe/src/TestDataCollector.cxx -o CMakeFiles/TestDataCollector.exe.dir/src/TestDataCollector.cxx.s

# Object files for target TestDataCollector.exe
TestDataCollector_exe_OBJECTS = \
"CMakeFiles/TestDataCollector.exe.dir/src/TestDataCollector.cxx.o"

# External object files for target TestDataCollector.exe
TestDataCollector_exe_EXTERNAL_OBJECTS =

main/exe/TestDataCollector.exe: main/exe/CMakeFiles/TestDataCollector.exe.dir/src/TestDataCollector.cxx.o
main/exe/TestDataCollector.exe: main/exe/CMakeFiles/TestDataCollector.exe.dir/build.make
main/exe/TestDataCollector.exe: main/lib/libEUDAQ.dylib
main/exe/TestDataCollector.exe: main/exe/CMakeFiles/TestDataCollector.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/gennai/DESY/euda53/bd/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable TestDataCollector.exe"
	cd /Users/gennai/DESY/euda53/bd/main/exe && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestDataCollector.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
main/exe/CMakeFiles/TestDataCollector.exe.dir/build: main/exe/TestDataCollector.exe

.PHONY : main/exe/CMakeFiles/TestDataCollector.exe.dir/build

main/exe/CMakeFiles/TestDataCollector.exe.dir/clean:
	cd /Users/gennai/DESY/euda53/bd/main/exe && $(CMAKE_COMMAND) -P CMakeFiles/TestDataCollector.exe.dir/cmake_clean.cmake
.PHONY : main/exe/CMakeFiles/TestDataCollector.exe.dir/clean

main/exe/CMakeFiles/TestDataCollector.exe.dir/depend:
	cd /Users/gennai/DESY/euda53/bd && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/gennai/DESY/euda53 /Users/gennai/DESY/euda53/main/exe /Users/gennai/DESY/euda53/bd /Users/gennai/DESY/euda53/bd/main/exe /Users/gennai/DESY/euda53/bd/main/exe/CMakeFiles/TestDataCollector.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : main/exe/CMakeFiles/TestDataCollector.exe.dir/depend

