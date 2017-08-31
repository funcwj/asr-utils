//
// Created by wujian on 17-8-30.
//

#ifndef FFTCOMPUTER_H
#define FFTCOMPUTER_H

#include <complex>
#include <vector>
#include <algorithm>

const float PI = 3.14159265;

typedef std::vector< std::complex<float> > ComplexVector;

class FFTComputer {
public:
    FFTComputer() {};
    void ComplexFFT(ComplexVector &sig, bool invert);
    void ComplexFFT(float *real_samples, float *image_samples,
               size_t num_samples, bool invert);
    void RealFFT(float *real_samples, size_t num_samples);
};


#endif //FFTCOMPUTER_H
