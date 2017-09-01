//
// Created by wujian on 17-8-30.
//

#include "fft-computer.h"

void FFTComputer::ComplexFFT(ComplexVector &sig, bool invert) {
    int N = sig.size();

    for (int i = 0, j = 0; i < N - 1; i++) {
        if(i < j) {
            swap(sig[i], sig[j]);
        }
        int k = N >> 1;
        while(j >= k) {
            j -= k;
            k = k >> 1;
        }
        j += k;
    }

    ComplexVector W(N / 2);

    for (int i = 0; i < W.size(); i++) {
        W[i] = std::polar(1.0f, 2 * PI * i / N);
        if (invert) W[i] = conj(W[i]);
    }

    for (int m = 1; m < N; m = m << 1) {
        for (int n = 0; n < N; n += (m << 1)) {
            int k = 0;
            for(int j = n; j < n + m; j++, k += N / (m << 1)) {
                std::complex<float> A = sig[j], B = sig[j + m];
                sig[j] = A + W[k] * B, sig[j + m] = A - W[k] * B;
            }
        }
    }

    if (invert) {
        for_each(sig.begin(), sig.end(), [=] (std::complex<float> &s) {s /= N;});
    }
}

void FFTComputer::ComplexFFT(float *real_samples, float *image_samples,
                             size_t num_samples, bool invert)
{
    int n, i, j, m, cnt, inc, k;
    float WR, WI, Ri, Ii, Rj, Ij;

    n = num_samples;
    float *R = real_samples, *I = image_samples;

    for(j = 0, i = 0; i < n - 1; i++) {
        if(i < j) {
            std::swap(R[i], R[j]);
            std::swap(I[i], I[j]);
        }
        m = n >> 1;
        while(j >= m) {
            j = j - m;
            m = m >> 1;
        }
        j = j + m;
    }

    m = 1;
    // 1, 2, 4 级
    while(m < n) {
        /*
            m = 1: [1, 2], [3, 4], [5, 6], [7, 8] 4
            m = 2: [1, 3], [2, 4], [5, 7], [6, 8] 2
            m = 4: [1, 5], [2, 6], [3, 7], [4, 8] 1
         */
        //printf("M = %d\n", m);
        cnt = 0, inc = n / (m << 1);
        // inc: 4 2 1
        // m  : 1 2 4
        // W递增inc
        while(cnt < inc) {
            // m = 1: 1 3 5 7
            // m = 2: 1 5
            // m = 4: 1
            i = cnt * m * 2;
            // W[0, n]: inc
            // 计算m次 迭代inc次
            for(int t = 0; t < m; t++, i++) {
                j = i + m;
                k = t * inc;
                // printf("[%3d, %3d] W[%3d, %3d]\n", i, j, k, nn);
                k == 0 ? WR = 1.0, WI = 0.0: WR = cos(PI * k * 2 / n), WI = sin(PI * k * 2 / n);
                if(invert) WI = - WI;
                //(R[i], I[i]) = (Ri, Ii) + W * (Rj, Ij)
                //(R[j], I[j]) = (Ri, Ii) - W * (Rj, Ij)
                Rj = R[j], Ij = I[j], Ri = R[i], Ii = I[i];
                R[i] = Ri + WR * Rj - WI * Ij, I[i] = Ii + WR * Ij + WI * Rj;
                R[j] = Ri - WR * Rj + WI * Ij, I[j] = Ii - WR * Ij - WI * Rj;
            }
            cnt++;
        }
        m = m << 1;
    }

    if (invert)
        for (i = 0; i < n; i++)
            R[i] = R[i] / n, I[i] = I[i] / n;
}

void FFTComputer::RealFFT(float *real_samples, size_t num_samples) {
    size_t n = num_samples >> 1;

    float *R = new float[n];
    float *I = new float[n];

    for (int i = 0; i < n; i++) {
        R[i] = real_samples[i * 2], I[i] = real_samples[i * 2 + 1];
    }

    ComplexFFT(R, I, n, false);
    float FR, FI, GR, GI, YR, YI, CYR, CYI, XR, XI, cosr, sinr;

    for (int r = 0; r < n; r++) {
        if(r == 0) {
            FR = R[r], FI = 0.0, GR = I[r], GI = 0.0;
        } else {
            YR  = R[r], YI = I[r];
            CYR = R[n - r], CYI = -I[n - r];
            FR  = (YR + CYR) / 2, FI = (YI + CYI) / 2;
            GR  = (YI - CYI) / 2, GI = (CYR - YR) / 2;
        }
        cosr = cos(r * PI / n);
        sinr = sin(r * PI / n);
        XR = FR + cosr * GR - sinr * GI;
        XI = FI + cosr * GI + sinr * GR;
        real_samples[r * 2] = XR, real_samples[r * 2 + 1] = XI;
    }
    delete[] R;
    delete[] I;
}
