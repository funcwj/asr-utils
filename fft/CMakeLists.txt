cmake_minimum_required(VERSION 3.5)
project(FFT)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES fft-computer-test.cc fft-computer.cc fft-computer.h)
add_executable(fft-computor-test ${SOURCE_FILES})
