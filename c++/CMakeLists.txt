cmake_minimum_required(VERSION 2.8)
project(LineDetect2)
message ("$ENV{CXXFLAGS}") 

#SET(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 $ENV{CXXFLAGS} -o0 -Wall  -ggdb  ")
#SET(CMAKE_CXX_FLAGS_RELEASE " -std=c++11 $ENV{CXXFLAGS} -O3 -Wall")
SET(CMAKE_CXX_FLAGS_DEBUG "-std=c++11  -DRAPIDJSON_HAS_STDSTRING   -O3 -Wall  -ggdb  ")
SET(CMAKE_CXX_FLAGS_RELEASE "  -DRAPIDJSON_HAS_STDSTRING   -std=c++0x  -O3 -Wall")

#find_package(OpenCV REQUIRED)
include_directories("${PROJECT_SOURCE_DIR}/")
add_definitions(" -D_LINUX   -std=c++0x ")
add_library(Line SHARED ConnTree.cpp DSCCLine.cpp DirLine.cpp intsigImage.cpp LineDetect.cpp)
#add_executable(LineDetect ConnTree.cpp DSCCLine.cpp DirLine.cpp main.cpp intsigImage.cpp LineDetect.cpp)
add_executable(LineDetect main.cpp)

# Find OpenCV, you may need to set OpenCV_DIR variable
# to the absolute path to the directory containing OpenCVConfig.cmake file
# via the command line or GUI
find_package(OpenCV REQUIRED)

# If the package has been found, several variables will
# be set, you can find the full list with descriptions
# in the OpenCVConfig.cmake file.
# Print some message showing some of them

message(STATUS "    OpenCV library status:")
message(STATUS "    version: ${OpenCV_VERSION}")
message(STATUS "    libraries: ${OpenCV_LIBS}")
message(STATUS "    include path: ${OpenCV_INCLUDE_DIRS}")

if(CMAKE_VERSION VERSION_LESS "2.8.11")
  # Add OpenCV headers location to your include paths
  include_directories(${OpenCV_INCLUDE_DIRS})
endif()

target_link_libraries(LineDetect ${OpenCV_LIBS} Line)

