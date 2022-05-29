#include <iostream>
#include "matrix.h"
#include "iterative.h"
#include "utils.h"

template<typename T>
void checkResult(std::vector<T> solutionVector, int stepCount){
    if (stepCount > 0) {
        std::cout << "迭代次数: " << stepCount << '\n';
        std::cout << "解向量如下: \n";
        printXi(solutionVector);
    } else {
        std::cout << "迭代过程中解的绝对值或迭代次数超过所指定的上限；可能原因：迭代不收敛，或所指定的上限值过小。" << '\n';
    }
}

template<typename T>
void test_Hilbert(int n) {
    auto solveByJacobi = [&](matrix<T> A, matrix<T> b) {
        std::cout << "Jacobi迭代求解结果如下：\n";
        auto[solutionVector, stepCount] = Jacobi<T>(A, b, 1e-8, 1e20, 1e6);
        checkResult(solutionVector, stepCount);
    };
    auto solveBySOR = [&](matrix<T> A, matrix<T> b, T mu) {
        std::cout << "参数为 " << mu << " 的SOR迭代求解结果如下：\n";
        auto[solutionVector, stepCount] = SOR<T>(A, b, mu, 1e-8, 1e20, 1e6);
        checkResult(solutionVector, stepCount);
    };
    auto hilbertMatrix = matrix<T>(n, n);
    for (size_t i = 1; i <= n; i++) {
        for (size_t j = 1; j <= n; j++) {
            hilbertMatrix[i - 1][j - 1] = 1.0 / (i + j - 1);
        }
    }
    matrix<T> hilbertBi = matrix<T>(n, 1);
    hilbertBi.fillMatrix(1);
    hilbertBi = hilbertMatrix * hilbertBi;
    std::cout << "【此时处理的是 " << n << " 阶Hilbert矩阵】\n";
    solveByJacobi(hilbertMatrix, hilbertBi);
    solveBySOR(hilbertMatrix, hilbertBi, 1);
    solveBySOR(hilbertMatrix, hilbertBi, 1.25);
    solveBySOR(hilbertMatrix, hilbertBi, 1.5);
}

int main() {
    test_Hilbert<long double>(6);
    test_Hilbert<long double>(8);
    test_Hilbert<long double>(10);
}
