//
// Created by wujian on 17-9-5.
//

#ifndef MAT_H
#define MAT_H

#include <iostream>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <limits>
#include <cmath>

#include "cblas.h"

#define ALIGNMENTS_PARAM 16

#define MEM_ALLOC(size) aligned_alloc(ALIGNMENTS_PARAM, size)
#define MEM_FREED(addr) free(addr)


namespace mat {

    typedef float FloatType;
    typedef int32_t int32;

    class Mat {

    public:
        // default
        Mat(): row_(0), col_(0), step_(0), data_(NULL), own_(false) {}

        // matrix from ptr
        Mat(FloatType *addr, int32 row, int32 col, bool copy = false);

        // vector from ptr
        Mat(FloatType *addr, int32 dim, bool copy = false);

        // subvector or submatrix
        Mat(Mat &mat, int32 base_idx, int32 extend_dim = 1);

        // vector
        Mat(int32 dim) { New(1, dim); }

        // matrix
        Mat(int32 row, int32 col) { New(row, col); }

        // destroy
        ~Mat() { Destroy(); }

        // memset with 0
        void SetZero();

        // fill with value
        void Set(FloatType value);

        // vector dot
        FloatType VecDot(Mat &mat);

        // vector norm2
        FloatType VecNorm(FloatType p);

        // vector scale
        void VecScale(FloatType scale);

        // vector sum
        FloatType VecSum();

        // vector max
        FloatType VecMin();

        // vector min
        FloatType VecMax();

        // vector add vector
        void VecAddVec(Mat &mat, FloatType alpha);

        // vector add matrix x vector
        void VecAddMatVec(FloatType alpha, Mat &mat, Mat &vec,
                          bool trans, FloatType beta);

        // matrix add matrix x matrix
        void MatAddMatMat(FloatType alpha, Mat &matA, bool transA,
                          Mat &matB, bool transB, FloatType beta);

        FloatType  operator() (int32 index) const;
        // modify vec(i) = x
        FloatType& operator() (int32 index);
        FloatType  operator() (int32 rindex, int32 cindex) const;
        // modify mat(i, j) = x
        FloatType& operator() (int32 rindex, int32 cindex);
        // copy and realloc memory
        Mat&       operator=  (const Mat &other);

        int32 Row() const { return row_; }
        int32 Col() const { return col_; }
        int32 Step() const { return step_; }
        int32 Size() const { return col_ * row_; }

        FloatType* Data() const { return data_; }

        bool IsVec() const { return row_ == 1; }

        std::string Type() { return IsVec() ? std::string("Vec"): std::string("Mat"); }
        std::string Info();

    private:
        void Destroy();
        void New(int32 row, int32 col, int32 pad = 0);

        void SubVec(Mat &mat, int32 from, int32 dim);
        void SubMat(Mat &mat, int32 from_row, int32 num_row);

        // new memory and copy data
        void CopyFromMat(const Mat &mat);

        // new memory and copy data
        void CopyFromPtr(FloatType *addr, int32 row, int32 col);

        // reference addr
        void RefFromPtr(FloatType *addr, int32 row, int32 col, int32 step = 0);

        int32 row_, col_, step_;

        // true means take own of the memory and need to
        // free it when ~Mat()
        bool  own_;
        FloatType *data_;
    };

    std::ostream & operator << (std::ostream & out, Mat &mat);
}

#endif // MAT_H
