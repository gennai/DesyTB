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
include main/exe/CMakeFiles/OptionExample.exe.dir/depend.make

# Include the progress variables for this target.
include main/exe/CMakeFiles/OptionExample.exe.dir/progress.make

# Include the compile flags for this target's objects.
include main/exe/CMakeFiles/OptionExample.exe.dir/flags.make

main/exe/CMakeFiles/OptionExample.exe.dir/src/OptionExample.cxx.o: main/exe/CMakeFiles/OptionExample.exe.dir/flags.make
main/exe/CMakeFiles/OptionExample.exe.dir/src/OptionExample.cxx.o: ../main/exe/src/OptionExample.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/gennai/DESY/euda53/bd/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object main/exe/CMakeFiles/OptionExample.exe.dir/src/OptionExample.cxx.o"
	cd /Users/gennai/DESY/euda53/bd/main/exe && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/OptionExample.exe.dir/src/OptionExample.cxx.o -c /Users/gennai/DESY/euda53/main/exe/src/OptionExample.cxx

main/exe/CMakeFiles/OptionExample.exe.dir/src/OptionExample.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/OptionExample.exe.dir/src/OptionExample.cxx.i"
	cd /Users/gennai/DESY/euda53/bd/main/exe && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/gennai/DESY/euda53/main/exe/src/OptionExample.cxx > CMakeFiles/OptionExample.exe.dir/src/OptionExample.cxx.i

main/exe/CMakeFiles/OptionExample.exe.dir/src/OptionExample.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/OptionExample.exe.dir/src/OptionExample.cxx.s"
	cd /Users/gennai/DESY/euda53/bd/main/exe && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/gennai/DESY/euda53/main/exe/src/OptionExample.cxx -o CMakeFiles/OptionExample.exe.dir/src/OptionExample.cxx.s

# Object files for target OptionExample.exe
OptionExample_exe_OBJECTS = \
"CMakeFiles/OptionExample.exe.dir/src/OptionExample.cxx.o"

# External object files for target OptionExample.exe
OptionExample_exe_EXTERNAL_OBJECTS =

main/exe/OptionExample.exe: main/exe/CMakeFiles/OptionExample.exe.dir/src/OptionExample.cxx.o
main/exe/OptionExample.exe: main/exe/CMakeFiles/OptionExample.exe.dir/build.make
main/exe/OptionExample.exe: main/lib/libEUDAQ.dylib
main/exe/OptionExample.exe: main/exe/CMakeFiles/OptionExample.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/gennai/DESY/euda53/bd/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable OptionExample.exe"
	cd /Users/gennai/DESY/euda53/bd/main/exe && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/OptionExample.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
main/exe/CMakeFiles/OptionExample.exe.dir/build: main/exe/OptionExample.exe

.PHONY : main/exe/CMakeFiles/OptionExample.exe.dir/build

main/exe/CMakeFiles/OptionExample.exe.dir/clean:
	cd /Users/gennai/DESY/euda53/bd/main/exe && $(CMAKE_COMMAND) -P CMakeFiles/OptionExample.exe.dir/cmake_clean.cmake
.PHONY : main/exe/CMakeFiles/OptionExample.exe.dir/clean

main/exe/CMakeFiles/OptionExample.exe.dir/depend:
	cd /Users/gennai/DESY/euda53/bd && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/gennai/DESY/euda53 /Users/gennai/DESY/euda53/main/exe /Users/gennai/DESY/euda53/bd /Users/gennai/DESY/euda53/bd/main/exe /Users/gennai/DESY/euda53/bd/main/exe/CMakeFiles/OptionExample.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : main/exe/CMakeFiles/OptionExample.exe.dir/depend

