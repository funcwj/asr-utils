#include <iostream>
#include "fft-computer.h"


int main() {
    float R[16] = {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float I[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    std::vector< std::complex<float> > s(R, R + 16);
    FFTComputer fftcomputer;

    fftcomputer.ComplexFFT(R, I, 16, false);
    fftcomputer.ComplexFFT(R, I, 16, true);
    for (size_t i = 0; i < 16; i++) {
        std::cout << "[" << static_cast<int>(R[i]) << ", " << static_cast<int>(I[i]) << "]" << std::endl;
    }

    fftcomputer.ComplexFFT(s, false);
    fftcomputer.ComplexFFT(s, true);
     for_each(s.begin(), s.end(), [](std::complex<float> s) {
         std::cout << "[" << static_cast<int>(s.real()) << ", "
                   << static_cast<int>(s.imag()) << "]" << std::endl;
     });

    return 0;
}