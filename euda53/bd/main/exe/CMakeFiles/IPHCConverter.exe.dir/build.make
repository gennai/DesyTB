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
include main/exe/CMakeFiles/IPHCConverter.exe.dir/depend.make

# Include the progress variables for this target.
include main/exe/CMakeFiles/IPHCConverter.exe.dir/progress.make

# Include the compile flags for this target's objects.
include main/exe/CMakeFiles/IPHCConverter.exe.dir/flags.make

main/exe/CMakeFiles/IPHCConverter.exe.dir/src/IPHCConverter.cxx.o: main/exe/CMakeFiles/IPHCConverter.exe.dir/flags.make
main/exe/CMakeFiles/IPHCConverter.exe.dir/src/IPHCConverter.cxx.o: ../main/exe/src/IPHCConverter.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/gennai/DESY/euda53/bd/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object main/exe/CMakeFiles/IPHCConverter.exe.dir/src/IPHCConverter.cxx.o"
	cd /Users/gennai/DESY/euda53/bd/main/exe && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/IPHCConverter.exe.dir/src/IPHCConverter.cxx.o -c /Users/gennai/DESY/euda53/main/exe/src/IPHCConverter.cxx

main/exe/CMakeFiles/IPHCConverter.exe.dir/src/IPHCConverter.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/IPHCConverter.exe.dir/src/IPHCConverter.cxx.i"
	cd /Users/gennai/DESY/euda53/bd/main/exe && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/gennai/DESY/euda53/main/exe/src/IPHCConverter.cxx > CMakeFiles/IPHCConverter.exe.dir/src/IPHCConverter.cxx.i

main/exe/CMakeFiles/IPHCConverter.exe.dir/src/IPHCConverter.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/IPHCConverter.exe.dir/src/IPHCConverter.cxx.s"
	cd /Users/gennai/DESY/euda53/bd/main/exe && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/gennai/DESY/euda53/main/exe/src/IPHCConverter.cxx -o CMakeFiles/IPHCConverter.exe.dir/src/IPHCConverter.cxx.s

# Object files for target IPHCConverter.exe
IPHCConverter_exe_OBJECTS = \
"CMakeFiles/IPHCConverter.exe.dir/src/IPHCConverter.cxx.o"

# External object files for target IPHCConverter.exe
IPHCConverter_exe_EXTERNAL_OBJECTS =

main/exe/IPHCConverter.exe: main/exe/CMakeFiles/IPHCConverter.exe.dir/src/IPHCConverter.cxx.o
main/exe/IPHCConverter.exe: main/exe/CMakeFiles/IPHCConverter.exe.dir/build.make
main/exe/IPHCConverter.exe: main/lib/libEUDAQ.dylib
main/exe/IPHCConverter.exe: main/exe/CMakeFiles/IPHCConverter.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/gennai/DESY/euda53/bd/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable IPHCConverter.exe"
	cd /Users/gennai/DESY/euda53/bd/main/exe && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/IPHCConverter.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
main/exe/CMakeFiles/IPHCConverter.exe.dir/build: main/exe/IPHCConverter.exe

.PHONY : main/exe/CMakeFiles/IPHCConverter.exe.dir/build

main/exe/CMakeFiles/IPHCConverter.exe.dir/clean:
	cd /Users/gennai/DESY/euda53/bd/main/exe && $(CMAKE_COMMAND) -P CMakeFiles/IPHCConverter.exe.dir/cmake_clean.cmake
.PHONY : main/exe/CMakeFiles/IPHCConverter.exe.dir/clean

main/exe/CMakeFiles/IPHCConverter.exe.dir/depend:
	cd /Users/gennai/DESY/euda53/bd && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/gennai/DESY/euda53 /Users/gennai/DESY/euda53/main/exe /Users/gennai/DESY/euda53/bd /Users/gennai/DESY/euda53/bd/main/exe /Users/gennai/DESY/euda53/bd/main/exe/CMakeFiles/IPHCConverter.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : main/exe/CMakeFiles/IPHCConverter.exe.dir/depend

