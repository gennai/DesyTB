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
include main/exe/CMakeFiles/FileChecker.exe.dir/depend.make

# Include the progress variables for this target.
include main/exe/CMakeFiles/FileChecker.exe.dir/progress.make

# Include the compile flags for this target's objects.
include main/exe/CMakeFiles/FileChecker.exe.dir/flags.make

main/exe/CMakeFiles/FileChecker.exe.dir/src/FileChecker.cxx.o: main/exe/CMakeFiles/FileChecker.exe.dir/flags.make
main/exe/CMakeFiles/FileChecker.exe.dir/src/FileChecker.cxx.o: ../main/exe/src/FileChecker.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/gennai/DESY/euda53/bd/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object main/exe/CMakeFiles/FileChecker.exe.dir/src/FileChecker.cxx.o"
	cd /Users/gennai/DESY/euda53/bd/main/exe && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FileChecker.exe.dir/src/FileChecker.cxx.o -c /Users/gennai/DESY/euda53/main/exe/src/FileChecker.cxx

main/exe/CMakeFiles/FileChecker.exe.dir/src/FileChecker.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FileChecker.exe.dir/src/FileChecker.cxx.i"
	cd /Users/gennai/DESY/euda53/bd/main/exe && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/gennai/DESY/euda53/main/exe/src/FileChecker.cxx > CMakeFiles/FileChecker.exe.dir/src/FileChecker.cxx.i

main/exe/CMakeFiles/FileChecker.exe.dir/src/FileChecker.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FileChecker.exe.dir/src/FileChecker.cxx.s"
	cd /Users/gennai/DESY/euda53/bd/main/exe && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/gennai/DESY/euda53/main/exe/src/FileChecker.cxx -o CMakeFiles/FileChecker.exe.dir/src/FileChecker.cxx.s

# Object files for target FileChecker.exe
FileChecker_exe_OBJECTS = \
"CMakeFiles/FileChecker.exe.dir/src/FileChecker.cxx.o"

# External object files for target FileChecker.exe
FileChecker_exe_EXTERNAL_OBJECTS =

main/exe/FileChecker.exe: main/exe/CMakeFiles/FileChecker.exe.dir/src/FileChecker.cxx.o
main/exe/FileChecker.exe: main/exe/CMakeFiles/FileChecker.exe.dir/build.make
main/exe/FileChecker.exe: main/lib/libEUDAQ.dylib
main/exe/FileChecker.exe: main/exe/CMakeFiles/FileChecker.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/gennai/DESY/euda53/bd/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable FileChecker.exe"
	cd /Users/gennai/DESY/euda53/bd/main/exe && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/FileChecker.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
main/exe/CMakeFiles/FileChecker.exe.dir/build: main/exe/FileChecker.exe

.PHONY : main/exe/CMakeFiles/FileChecker.exe.dir/build

main/exe/CMakeFiles/FileChecker.exe.dir/clean:
	cd /Users/gennai/DESY/euda53/bd/main/exe && $(CMAKE_COMMAND) -P CMakeFiles/FileChecker.exe.dir/cmake_clean.cmake
.PHONY : main/exe/CMakeFiles/FileChecker.exe.dir/clean

main/exe/CMakeFiles/FileChecker.exe.dir/depend:
	cd /Users/gennai/DESY/euda53/bd && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/gennai/DESY/euda53 /Users/gennai/DESY/euda53/main/exe /Users/gennai/DESY/euda53/bd /Users/gennai/DESY/euda53/bd/main/exe /Users/gennai/DESY/euda53/bd/main/exe/CMakeFiles/FileChecker.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : main/exe/CMakeFiles/FileChecker.exe.dir/depend

