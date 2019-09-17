# stop if cmake version is below 3.0
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# project name and enable C++ support
project(mzipsat CXX C)

set(CMAKE_BUILD_TYPE Release)

add_subdirectory(src)
add_subdirectory(example)

set_target_properties(dipoleamplitude PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
set_target_properties(mzipsat PROPERTIES
     LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
