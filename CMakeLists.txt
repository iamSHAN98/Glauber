cmake_minimum_required(VERSION 3.10)
project(Glauber VERSION 1.0 DESCRIPTION "Glauber Event Generator" LANGUAGES CXX)

# Set compile options
set(CMAKE_CXX_STANDARD_REQUIRED True)
if(NOT CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD STREQUAL "")
  set(CMAKE_CXX_STANDARD 14)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb3")
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# Cross-check compile options
message("Using C++ standard ${CMAKE_CXX_STANDARD}")
message("Build mode : ${CMAKE_BUILD_TYPE}")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message("Build options : ${CMAKE_CXX_FLAGS_DEBUG}")
else()
  message("Build options : ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

# External packages
# GSL : Integration & random number generation routines
find_package(GSL)

# YAML-CPP : Configuration parsing            
add_subdirectory(yaml-cpp)
set(YAML_CPP_INCLUDE yaml-cpp/include)

# DataStream : Data storage
add_subdirectory(DataStream)

# Project source
include_directories(include)
file(GLOB SOURCE "src/*.cc")

# Compilation
add_executable(Glauber main.cc ${SOURCE})
target_include_directories(Glauber PRIVATE ${GSL_INCLUDE_DIRS} 
                                           ${YAML_CPP_INCLUDE}
                                           ${DataStream_INCLUDE})
target_link_libraries(Glauber PRIVATE ${GSL_LIBRARIES}
                                      yaml-cpp
                                      DataStream)