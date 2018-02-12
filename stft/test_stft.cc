// test_stft.cc
// wujian@18.2.12

#include "matrix/matrix-lib.h"
#include "base/kaldi-common.h"
#include "feat/wave-reader.h"
#include "util/common-utils.h"
#include "stft.h"

using namespace kaldi;

int test_stft() {

    bool binary;
    Input wave_in("orig.wav", &binary);        
    WaveData wave_orig;
    wave_orig.Read(wave_in.Stream());
    
    // configs
    ShortTimeFTOptions opts;
    opts.frame_length = 1024;
    opts.frame_shift  = 256;
    opts.normalize    = false;
    opts.track_volume = true;


    ShortTimeFTComputer stft_computer(opts);

    Matrix<BaseFloat> specs;
    stft_computer.ShortTimeFT(wave_orig.Data(), &specs);

    Matrix<BaseFloat> recon;
    stft_computer.InverseShortTimeFT(specs, &recon);

    Output ko("copy.wav", binary, false);
    WaveData wave_copy(16000, recon);
    wave_copy.Write(ko.Stream());
    // std:: cout << vec << std::endl;
}

void test_realfft() {
    int32 dim = 16;

    Vector<BaseFloat> vec(dim);
    vec.SetRandn();
    std::cout << vec;
    RealFft(&vec, true);
    std::cout << vec;
    RealFft(&vec, false);
    vec.Scale(1.0 / dim);
    std::cout << vec;
}

void test_istft() {

    bool binary;
    Input wave_in("orig.wav", &binary);        
    WaveData wave_orig;
    wave_orig.Read(wave_in.Stream());
    
    // configs
    ShortTimeFTOptions opts;
    opts.frame_length = 1024;
    opts.frame_shift  = 256;
    opts.normalize    = false;
    opts.track_volume = true;
    opts.apply_log    = true;
    opts.power_spectrum = true;


    ShortTimeFTComputer stft_computer(opts);

    Matrix<BaseFloat> stft_orig;
    stft_computer.ShortTimeFT(wave_orig.Data(), &stft_orig);
    std::cout << stft_orig.Row(5);

    Matrix<BaseFloat> specs, angle, stft_recon;
    stft_computer.ComputeSpectrum(stft_orig, &specs);
    stft_computer.ComputeArg(stft_orig, &angle);
    stft_computer.RestoreShortTimeFT(specs, angle, &stft_recon);
    std::cout << stft_recon.Row(5);
    
    Matrix<BaseFloat> recon;
    stft_computer.InverseShortTimeFT(stft_recon, &recon);

    Output ko("copy.wav", binary, false);
    WaveData wave_copy(16000, recon);
    wave_copy.Write(ko.Stream());
    // std:: cout << vec << std::endl;
}


int main() {
    test_istft();
    return 0;
}
