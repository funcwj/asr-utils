cmake_minimum_required(VERSION 3.5)
project(WaveReader)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES test-reader.cc wave.cc wave.h)
add_executable(test-reader ${SOURCE_FILES})