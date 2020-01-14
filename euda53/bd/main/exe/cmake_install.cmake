# Install script for directory: /Users/gennai/DESY/euda53/main/exe

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/gennai/DESY/euda53")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/ClusterExtractor.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ClusterExtractor.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ClusterExtractor.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ClusterExtractor.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ClusterExtractor.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ClusterExtractor.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/Converter.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Converter.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Converter.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Converter.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Converter.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Converter.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/ExampleProducer.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleProducer.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleProducer.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleProducer.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleProducer.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleProducer.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/ExampleSlowProducer.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleSlowProducer.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleSlowProducer.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleSlowProducer.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleSlowProducer.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleSlowProducer.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/ExampleReader.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleReader.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleReader.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleReader.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleReader.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExampleReader.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/FileChecker.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FileChecker.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FileChecker.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FileChecker.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FileChecker.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FileChecker.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/IPHCConverter.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/IPHCConverter.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/IPHCConverter.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/IPHCConverter.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/IPHCConverter.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/IPHCConverter.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/MagicLogBook.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MagicLogBook.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MagicLogBook.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MagicLogBook.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MagicLogBook.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MagicLogBook.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/OptionExample.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/OptionExample.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/OptionExample.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/OptionExample.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/OptionExample.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/OptionExample.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/RunListener.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/RunListener.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/RunListener.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/RunListener.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/RunListener.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/RunListener.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/TestDataCollector.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestDataCollector.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestDataCollector.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestDataCollector.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestDataCollector.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestDataCollector.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/TestLogCollector.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestLogCollector.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestLogCollector.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestLogCollector.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestLogCollector.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestLogCollector.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/TestMonitor.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestMonitor.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestMonitor.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestMonitor.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestMonitor.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestMonitor.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/TestProducer.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestProducer.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestProducer.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestProducer.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestProducer.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestProducer.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/TestReader.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestReader.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestReader.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestReader.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestReader.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestReader.exe")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/gennai/DESY/euda53/bd/main/exe/TestRunControl.exe")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestRunControl.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestRunControl.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/gennai/DESY/euda53/bd/main/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestRunControl.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/gennai/DESY/euda53/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestRunControl.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/TestRunControl.exe")
    endif()
  endif()
endif()

