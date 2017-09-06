//
// Created by wujian on 17-9-5.
//

#include "mat.h"

namespace mat {

    Mat::Mat(Mat &mat, int32 base_idx, int32 extend_dim) {
        // vector => subvec
        // matrix => submat
        mat.IsVec() ? SubVec(mat, base_idx, extend_dim): SubMat(mat, base_idx, extend_dim);
    }

    Mat::Mat(FloatType *addr, int32 dim, bool copy) {
        copy ? CopyFromPtr(addr, 1, dim): RefFromPtr(addr, 1, dim);
    }

    Mat::Mat(FloatType *addr, int32 row, int32 col, bool copy) {
        copy ?  CopyFromPtr(addr, row, col): RefFromPtr(addr, row, col);
    }


    void Mat::SubVec(Mat &mat, int32 from, int32 dim) {
        assert(from + dim <= mat.Col());
        FloatType *base_addr = mat.Data();
        RefFromPtr(base_addr + from, 1, dim);
    }

    void Mat::SubMat(Mat &mat, int32 from_row, int32 num_row) {
        assert(from_row + num_row <= mat.Row());
        FloatType *base_addr = mat.Data();
        RefFromPtr(base_addr + from_row * mat.Step(), num_row, mat.Col(), mat.Step());
    }

    // simple: first destroy, then new
    Mat& Mat::operator = (const Mat &other) {
        Destroy();
        CopyFromMat(other);
        return *this;
    }

    void Mat::CopyFromMat(const Mat &mat) {
        New(mat.Row(), mat.Col(), mat.Step());
        FloatType *src_addr = mat.Data();
        FloatType *dst_addr = data_;
        std::clog << "Copy " << mat.Row() * mat.Step() << " floats" << std::endl;
        memcpy(dst_addr, src_addr, sizeof(FloatType) * mat.Row() * mat.Step());
    }

    void Mat::CopyFromPtr(FloatType *addr, int32 row, int32 col) {
        Destroy();
        New(row, col, 0);
        if (row == 1) {
            memcpy(data_, addr, sizeof(FloatType) * col);
            return;
        }
        for (int32 i = 0; i < row_; i++) {
            for (int32 j = 0; j < col_; j++)
                memcpy(data_ + i * step_, addr + i * col, sizeof(FloatType) * col);
        }
    }

    void Mat::Destroy() {
        if (data_ != NULL && own_) {
            std::clog << "Free memory at " << data_ << std::endl;
            MEM_FREED(reinterpret_cast<void*>(data_));
        }
        col_ = row_ = step_ = 0;
    }

    void Mat::SetZero() {
        if (col_ == step_) {
            memset(data_, 0, sizeof(FloatType) * col_ * row_);
            return;
        }
        for (int32 i = 0; i < row_; i++)
            memset(data_ + i * step_, 0, sizeof(FloatType) * col_);
    }

    void Mat::Set(FloatType value) {
        for (int32 i = 0; i < row_; i++)
            for(int j = 0; j < col_; j++)
                (*this)(i, j) = value;
    }

    // vec dot vec
    FloatType Mat::VecDot(Mat &mat) {
        assert(mat.IsVec() && IsVec());
        assert(mat.Size() == Size());
        int32 dim = static_cast<int32>(Size());
        assert(dim);
        FloatType dot = cblas_sdot(dim, data_, 1, mat.Data(), 1);
        return dot;
    }

    // vec p-norm
    FloatType Mat::VecNorm(FloatType p) {
        assert(p >= 1.0 && IsVec());
        int32 dim = static_cast<int32>(Size());
        assert(dim);
        FloatType norm;
        // sum_i abs(y[i])
        if (p == 1.0)
            norm = cblas_sasum(dim, data_, 1);
        // sqrt(sum_i y[i]^2)
        else if (p == 2.0)
            norm = cblas_snrm2(dim, data_, 1);
        // max_i abs(y[i])
        else if (p == std::numeric_limits<FloatType>::infinity()) {
            norm = -p;
            for (int32 i = 0; i < dim; i++)
                norm = std::max(norm, std::abs((*this)(i)));
        } else {
            norm = 0.0;
            for (int32 i = 0; i < dim; i++)
                norm += pow(std::abs((*this)(i)), p);
            norm = pow(norm, 1.0 / p);
        }
        return norm;
    }

    FloatType Mat::VecMax() {
        assert(IsVec());
        int32 dim = static_cast<int32>(Size());
        assert(dim);
        FloatType maxn = -std::numeric_limits<FloatType>::infinity();
        for (int32 i = 0; i < dim; i++)
            maxn = std::max(maxn, (*this)(i));
        return maxn;
    }

    FloatType Mat::VecMin() {
        assert(IsVec());
        int32 dim = static_cast<int32>(Size());
        assert(dim);
        FloatType minn = std::numeric_limits<FloatType>::infinity();
        for (int32 i = 0; i < dim; i++)
            minn = std::min(minn, (*this)(i));
        return minn;
    }

    // y = y * scale
    void Mat::VecScale(FloatType scale) {
        assert(IsVec());
        int32 dim = static_cast<int32>(Size());
        assert(dim);
        cblas_sscal(dim, scale, data_, 1);
    }

    // sum(v[i])
    FloatType Mat::VecSum() {
        assert(IsVec());
        int32 dim = static_cast<int32>(Size());
        assert(dim);
        FloatType sum = 0.0;
        for (int i = 0; i < dim; i++)
            sum += (*this)(i);
        return sum;
    }

    // y += alpha * v
    void Mat::VecAddVec(Mat &vec, FloatType alpha) {
        assert(vec.IsVec() && IsVec());
        assert(vec.Size() == Size());
        int32 dim = static_cast<int32>(Size());
        assert(dim);
        cblas_saxpy(dim, alpha, vec.Data(), 1, data_, 1);
    }

    // y += bata * y + alpha * Mv
    void Mat::VecAddMatVec(FloatType alpha, Mat &mat, Mat &vec,
                           bool trans, FloatType beta) {
        assert(IsVec() && vec.IsVec());
        std::clog << Size() << " = " << "[" << mat.Row() << " x " << mat.Col() << "] * "
                << vec.Size() << std::endl;
        assert((!trans && mat.Col() == vec.Size() && mat.Row() == Size()) ||
               (trans  && mat.Row() == vec.Size() && mat.Col() == Size()));
        assert(&vec != this);
        CBLAS_TRANSPOSE trans_type = trans ? CblasTrans: CblasNoTrans;
        cblas_sgemv(CblasRowMajor, trans_type, mat.Row(), mat.Col(), alpha,
                    mat.Data(), mat.Step(), vec.Data(), 1, beta, data_, 1);
    }

    // C = beta * AB + beta * C
    void Mat::MatAddMatMat(FloatType alpha, Mat &matA, bool transA,
                           Mat &matB, bool transB, FloatType beta) {
        assert(&matA != this && &matB != this);
        assert(Size());
        // AC-AR BC-BR = AC-BR == CR-CC
        if (transA && transB)
            assert(matA.Col() == row_ && matB.Row() == col_ && matA.Row() == matB.Col());
        // AC-AR BR-BC = AC-BC == CR-CC
        if (transA && !transB)
            assert(matA.Col() == row_ && matB.Col() == col_ && matA.Row() == matB.Row());
        // AR-AC BC-BR = AR-BR == CR-CC
        if (!transA && transB)
            assert(matA.Row() == row_ && matB.Row() == col_ && matA.Col() == matB.Col());
        // AR-AC BR-BC = AR-BC == CR-CC
        if (!transA && !transB)
            assert(matA.Row() == row_ && matB.Col() == col_ && matA.Col() == matB.Row());
        CBLAS_TRANSPOSE trans_type_a = transA ? CblasTrans: CblasNoTrans;
        CBLAS_TRANSPOSE trans_type_b = transB ? CblasTrans: CblasNoTrans;

        // std::clog << "cblas_sgemm: a_step = " << matA.Step() << "; " <<
                  // "b_step = " << matB.Step() << "; c_step = " << step_ << std::endl;
        cblas_sgemm(CblasRowMajor, trans_type_a, trans_type_b, row_, col_,
                    transA ? matA.Row(): matA.Col(), alpha, matA.Data(),
                    matA.Step(), matB.Data(), matB.Step(), beta, data_, step_);
    }

    // set private data
    void Mat::New(const int32 row, const int32 col, const int32 step) {
        col_ = col, row_ = row;
        own_ = true;
        assert(col_ + row_);
        int32 float_size = sizeof(FloatType);
        // step = 0 means we need to calculate it
        step_ = (step == 0) ? col_ + (float_size - col_ % float_size): step;
        void* addr = MEM_ALLOC(row_ * step_ * float_size);
        assert(addr != NULL);
        std::clog << "New memory [" << row_ << " X " << step_ << "] at " << addr << std::endl;
        data_ = reinterpret_cast<FloatType*>(addr);
    }

    // set private data
    void Mat::RefFromPtr(FloatType *addr, int32 row, int32 col, int32 step) {
        own_ = false;
        data_ = addr;
        row_ = row;
        col_ = col;
        step_ = (step == 0 ? col: step);
    }

    FloatType Mat::operator()(int32 index) const {
        assert(IsVec() && index < col_);
        return *(data_ + index);
    }

    FloatType& Mat::operator()(int32 index) {
        assert(IsVec() && index < col_);
        return *(data_ + index);
    }

    FloatType Mat::operator()(int32 rindex, int32 cindex) const {
        assert(!IsVec());
        assert(rindex < row_ && cindex < col_);
        return *(data_ + rindex * step_ + cindex);
    }

    FloatType& Mat::operator()(int32 rindex, int32 cindex) {
        assert(!IsVec());
        assert(rindex < row_ && cindex < col_);
        return *(data_ + rindex * step_ + cindex);
    }

    std::string Mat::Info() {
        std::ostringstream ss;
        ss << Type() << ": ";
        IsVec() ? ss << Size(): ss << Row() << " X " << Col();
        ss << " (step = " << Step() << ")";
        return ss.str();
    }

    std::ostream & operator << (std::ostream & out, Mat &mat) {
        if (mat.IsVec()) {
            std::cout << "[";
            for (int32 i = 0; i < mat.Size(); i++)
                std::cout << " " << mat(i);
            std::cout << " ]";
        } else {
            std::cout << "[";
            for (int32 i = 0; i < mat.Row(); i++) {
                Mat vec(mat, i);
                std::cout << (i == 0 ? "": " ") << vec;
                std::cout << (i == mat.Row() - 1 ? "": "\n");
            }
            std::cout << "]";
        }
    }

}