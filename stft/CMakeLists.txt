cmake_minimum_required(VERSION 3.5)
project(test)

set(CMAKE_CXX_STANDARD 11)
set(KALDI_DIR /gruntdata/jianwu.wj/kaldi)
add_definitions(-O3 -g -std=c++11 -DHAVE_CLAPACK)
include_directories(${KALDI_DIR}/tools/CLAPACK ${KALDI_DIR}/src ${KALDI_DIR}/tools/openfst/include)
link_directories(${KALDI_DIR}/src/lib)

add_executable(test_stft test_stft.cc stft.cc)

target_link_libraries(test_stft kaldi-base kaldi-util kaldi-matrix kaldi-feat)

