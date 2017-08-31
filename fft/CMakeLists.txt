cmake_minimum_required(VERSION 3.5)
project(FFT)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES fft-computor-test.cpp fft-computer.cpp fft-computer.h)
add_executable(fft-computor-test ${SOURCE_FILES})