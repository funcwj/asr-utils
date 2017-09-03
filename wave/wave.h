//
// Created by wujian on 17-9-3.
//

#ifndef WAVE_H
#define WAVE_H


#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstring>


typedef uint8_t  uint08;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef int8_t   int08;
typedef int16_t  int16;
typedef int32_t  int32;

#define CHAR_CAST(ptr) reinterpret_cast<char *>(ptr)

namespace wave {

    // sizeof(WaveInfo) == 36
    struct WaveInfo {
        int08  chunk_id[4];     // "RIFF"
        uint32 chunk_size;      // total bit size exclude chunk_id & chunk_size
        int08  format[4];       // "WAVE"
        int08  format_id[4];    // "fmt "
        uint32 format_size;     // format header total bit size exclude format_id & format_size
        int16  audio_format;    // 1
        int16  num_channels;    // 1/2
        uint32 sample_rate;     // 16000
        uint32 byte_rate;       // sample_rate * num_channel * sizeof(T)
        int16  block_align;     // num_channels * sizeof(T)
        int16  bits_per_sample; // sizeof(T) * 8
        // ...data part
    };

    template <typename T>
    class Wave {
    public:
        Wave(std::string filename): data_(NULL), num_bytes_(0) { Read(filename); };
        Wave(): data_(NULL), num_bytes_(0) {}

        ~Wave() { Destroy(); };

        void Read(std::string filename);
        void Write(std::string filename, T* data, uint32 num_samples,
                   int16 num_channels = 1, uint32 sample_rate = 16000);

        std::string Info() const;

        T* Data() const { return data_; }
        uint32 NumSamples() const { return num_bytes_ / sizeof(T); }
        uint32 SampleFreq() const { return info_.sample_rate; }
        int16 NumChannels() const { return info_.num_channels; }


    private:
        void ReadInfo(std::istream &is);
        void ReadData(std::istream &is);
        void ReadBinary(std::istream &is, char* ptr, uint32 num_bits);

        void Destroy() { if (data_ != NULL) delete[] data_; };

        WaveInfo info_;
        uint32 num_bytes_;
        T *data_;

    };

    template<typename T>
    std::ostream & operator << (std::ostream &out, const Wave<T> &wave) {
        out << wave.Info() << "Data           :\n\t[ ";
        T *data = wave.Data();
        for (int32 i = 0; i < wave.NumSamples(); i++)
            out << data[i] << " ";
        out << "]" << std::endl;
        return out;
    }
}





#endif // WAVE_H
