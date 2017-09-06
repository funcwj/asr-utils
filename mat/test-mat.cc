#include <iostream>
#include "mat.h"

void BasicTest() {
    float data[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    mat::Mat mat1(data, 3, 3);
    std::cout << "mat1:\n";
    std::cout << mat1 << std::endl;
    std::cout << "mat1 info:\n";
    std::cout << mat1.Info() << std::endl;
    mat::Mat vec0(data, 9);
    std::cout << vec0 << std::endl;
    std::cout << "vec0 info:\n";
    std::cout << vec0.Info() << std::endl;
    mat::Mat tmp;
    tmp = mat1;
    std::cout << "tmp info:\n";
    std::cout << tmp.Info() << std::endl;
    std::cout << "tmp = mat1:\n";
    std::cout << tmp << std::endl;
    mat::Mat vec1(mat1, 2);
    std::cout << "vec1(mat1, 2):\n";
    std::cout << vec1 << std::endl;
    std::cout << "after set vec1(2) = 10:\n";
    vec1(2) = 10;
    std::cout << vec1 << std::endl;
    mat1(1, 1) = 23;
    std::cout << "after set mat1(1, 1) = 23:\n";
    std::cout << mat1 << std::endl;
    std::cout << "tmp:\n";
    std::cout << tmp << std::endl;
}

void VVMathTest() {
    float data1[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    float data2[9] = {-2, -3, 4, 5, 6, 7, 8, 9, 10};

    mat::Mat vec1(data1, 9);
    mat::Mat vec2(data2, 9);
    std::cout << "vec1: " << vec1 << std::endl;
    std::cout << "vec2: " << vec2 << std::endl;
    std::cout << "vec1 * vec2 = " << vec1.VecDot(vec2) << std::endl;
    std::cout << "vec1 norm1  = " << vec1.VecNorm(1) << std::endl;
    std::cout << "vec1 norm2  = " << vec1.VecNorm(2) << std::endl;
    std::cout << "vec1 normi  = " << vec1.VecNorm(std::numeric_limits<float>::infinity()) << std::endl;

    vec1.VecScale(0.5);
    std::cout << "vec1 * 0.5  = " << vec1 << std::endl;
    vec1.VecAddVec(vec2, 1);
    std::cout << "vec1 + vec2 = " << vec1 << std::endl;
    std::cout << "sum(vec2)   = " << vec2.VecSum() << std::endl;
    std::cout << "vec2: " << vec2 << std::endl;
    std::cout << "vec2 norm1  = " << vec2.VecNorm(1) << std::endl;
    std::cout << "vec2 norm2  = " << vec2.VecNorm(2) << std::endl;
    std::cout << "vec2 normi  = " << vec2.VecNorm(std::numeric_limits<float>::infinity()) << std::endl;
    std::cout << "vec2 max  = " << vec2.VecMax() << std::endl;
    std::cout << "vec2 min  = " << vec2.VecMin() << std::endl;

}


void MVMathTest() {
    float mat_data[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    // 2 x 4
    mat::Mat mat(mat_data, 2, 4);
    std::cout << "mat: " << std::endl << mat << std::endl;
    float vec1_data[4] = {1, 2, 3, 4};
    mat::Mat vec1(vec1_data, 4);
    std::cout << "vec1: " << std::endl << vec1 << std::endl;
    mat::Mat ans1(2);
    ans1.SetZero();
    ans1.VecAddMatVec(0.5, mat, vec1, false, 0);
    std::cout << "0.5 * mat * vec1 = " << ans1 << std::endl;
    ans1.VecAddMatVec(0.5, mat, vec1, false, 0.5);
    std::cout << "0.5 * mat * vec1 + 0.5 * ans = " << ans1 << std::endl;
}


void MMMathTest() {
    // 2 x 3
    float data_mat1[6] = {2, 4, 6, 1, 3, 5};
    // 3 x 2
    float data_mat2[6] = {1, 2, 3, 4, 5, 6};
    float data_ans[4]  = {2, 3, 4, 5};

    mat::Mat mat1(data_mat1, 2, 3);
    mat::Mat mat2(data_mat2, 3, 2);
    // mat::Mat ans1(data_ans, 2, 2);
    mat::Mat ans1(2, 2);
    ans1.Set(2.5);
    std::cout << "mat1: " << std::endl << mat1 << std::endl;
    std::cout << "mat2: " << std::endl << mat2 << std::endl;
    std::cout << "ans1: " << std::endl << ans1 << std::endl;
    ans1.MatAddMatMat(0.5, mat1, false, mat2, false, 0.5);
    std::cout << "0.5 * mat1 * mat2 + ans1 * 0.5: " << std::endl << ans1 << std::endl;

    mat::Mat ans2(3, 3);
    ans2.Set(1);
    std::cout << "ans2: " << std::endl << ans2 << std::endl;
    ans2.MatAddMatMat(0.5, mat1, true, mat2, true, 0.5);
    std::cout << "0.5 * mat2^T * mat1^T + ans2 * 0.5: " << std::endl << ans2 << std::endl;

}


int main() {
    BasicTest();
    VVMathTest();
    MVMathTest();
    MMMathTest();
    return 0;
}