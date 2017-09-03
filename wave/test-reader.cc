#include "wave.h"

int main() {
    wave::Wave<int16> wave;
    wave.Read("test.wav");
    wave.Write("copy.wav", wave.Data(), wave.NumSamples());
    return 0;
}
