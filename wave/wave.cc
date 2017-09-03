//
// Created by wujian on 17-9-3.
//

#include "wave.h"

namespace wave {

    template<typename T>
    void Wave<T>::ReadData(std::istream &is) {
        // seek data part, cause there could be
        // other parameters between format part and data part
        is.seekg(info_.format_size - 16, std::ios::cur);

        char data_id[4];
        ReadBinary(is, data_id, sizeof(data_id));
        assert(strncmp(data_id, "data", sizeof(data_id)) == 0);

        ReadBinary(is, CHAR_CAST(&num_bytes_), sizeof(num_bytes_));
        data_ = new T[num_bytes_ / sizeof(T)];
        ReadBinary(is, CHAR_CAST(data_), num_bytes_);
        assert(info_.chunk_size == info_.format_size + 20 + num_bytes_);
    }

    template<typename T>
    void Wave<T>::ReadInfo(std::istream &is) {
        ReadBinary(is, CHAR_CAST(&info_), sizeof(info_));
        // check
        assert(info_.bits_per_sample / 8 == sizeof(T));
        assert(strncmp(CHAR_CAST(info_.chunk_id), "RIFF", 4) == 0);
        assert(strncmp(CHAR_CAST(info_.format), "WAVE", 4) == 0);
        assert(strncmp(CHAR_CAST(info_.format_id), "fmt ", 4) == 0);
        assert(info_.byte_rate == info_.sample_rate * info_.num_channels * sizeof(T));
        assert(info_.block_align == info_.num_channels * sizeof(T));
        assert(info_.format_size - 16 >= 0);
    }

    template<typename T>
    void Wave<T>::Read(std::string filename) {
        Destroy();
        std::ifstream input(filename);
        ReadInfo(input);
        ReadData(input);
        input.close();
    }

    template<typename T>
    void Wave<T>::ReadBinary(std::istream &is, char *ptr, uint32 num_bits) {
        is.read(ptr, num_bits);

        if (is.fail()) {
            std::cerr << "Read failure when loading " << num_bits << " bits" << std::endl;
            abort();
        }
    }

    template<typename T>
    void Wave<T>::Write(std::string filename, T* data,
                        uint32 num_samples, int16 num_channels,
                        uint32 sample_rate) {
        uint32 data_bitsize  = num_samples * sizeof(T);
        uint32 chunk_size    = data_bitsize + sizeof(WaveInfo);
        std::ofstream os(filename);
        WaveInfo info = {
            'R', 'I', 'F', 'F',
            chunk_size,
            'W', 'A', 'V', 'E',
            'f', 'm', 't', ' ',
            16, 1, num_channels,
            sample_rate,
            sample_rate * num_channels * sizeof(T),
            num_channels * sizeof(T),
            sizeof(T) * 8
        };
        // write chunk & format header
        os.write(CHAR_CAST(&info), sizeof(info));
        // write data header
        os.write("data", 4);
        os.write(CHAR_CAST(&data_bitsize), sizeof(data_bitsize));
        // write raw data
        os.write(CHAR_CAST(data), data_bitsize);
        os.close();
    }

    template<typename T>
    std::string Wave<T>::Info() const {
        std::stringstream ss;
        ss << "Num of channels: " << NumChannels() << "\n" <<
           "Sample rate    : " << SampleFreq() << "\n" <<
           "Bits per sample: " << sizeof(T) << "\n" <<
           "Num of samples : " << NumSamples() << "\n";
        std::string info = ss.str();
        return info;
    }


    template class Wave<int08>;
    template class Wave<int16>;
    template class Wave<int32>;

}