#include <iostream>
#include "matrix.h"
#include "iterative.h"
#include "utils.h"

template<typename T, typename F>
void testIteration(T start, const F &fun, T eps, T allowMaxValue, std::initializer_list<T> &xList) {
    try {
        auto[result, count] = functionIteration<T>(start, fun, eps, allowMaxValue, xList);
        std::cout << "迭代结果为: " << result << '\n';
        std::cout << "迭代次数为: " << count << '\n';
    } catch (std::domain_error &e) {
        std::cout << "Exception caught: " << e.what() << '\n';
    }

}

template<typename T>
void testEx1() {
    static const T zero_eps = 1e-10;
    std::cout << "【不动点迭代法 - 1】\n";
    std::initializer_list<T> fixedPointDomain_1 = {-1e20, -zero_eps, zero_eps, 1e20};
    testIteration((T) 1.5, [&](T x) { return 1 + 1.0 / (x * x); }, (T) 1e-15, (T) 1e20, fixedPointDomain_1);
    std::initializer_list<T> fixedPointDomain_2 = {1, 1e20};
    std::cout << "【不动点迭代法 - 2】\n";
    testIteration((T) 1.5, [&](T x) { return 1.0 / sqrt(x - 1); }, (T) 1e-15, (T) 1e20, fixedPointDomain_2);
    std::cout << "【牛顿迭代法】\n";
    std::initializer_list<T> newtonDomain = {-1e20, -zero_eps, zero_eps, 2.0 / 3, 2.0 / 3 + zero_eps, 1e20};
    testIteration((T) 1.5, [&](T x) { return x - (x * x * x - x * x - 1) / (3 * x * x - 2 * x); }, (T) 1e-15, (T) 1e20, newtonDomain);
}

template<typename T>
void solveEx2(matrix<T> A, matrix<T> b) {
    auto checkResult = [&](std::vector<T> solutionVector, int stepCount) {
        if (stepCount > 0) {
            std::cout << "迭代次数: " << stepCount << '\n';
            std::cout << "解向量如下: \n";
            printXi(solutionVector, 10);
        } else {
            std::cout << "迭代过程中解的绝对值或迭代次数超过所指定的上限；可能原因：迭代不收敛，或所指定的上限值过小。" << '\n';
        }
    };
    auto solveByJacobi = [&](matrix<T> A, matrix<T> b) {
        std::cout << "Jacobi迭代求解结果如下：\n";
        auto[solutionVector, stepCount] = Jacobi<T>(A, b, 1e-4, 1e20, 1e9);
        checkResult(solutionVector, stepCount);
    };
    auto solveBySOR = [&](matrix<T> A, matrix<T> b, T mu) {
        std::cout << "Gauss-Seidel迭代求解结果如下：\n";
        auto[solutionVector, stepCount] = SOR<T>(A, b, mu, 1e-4, 1e20, 1e9);
        checkResult(solutionVector, stepCount);
    };
    solveByJacobi(A, b);
    solveBySOR(A, b, 1);
}

template<typename T>
void testEx2() {
    matrix<T> A, b;
    A.readMatrix();
    b.readMatrix();
    solveEx2(A, b);
}

int main() {
    testEx1<long double>();
    testEx2<long double>();
    return 0;
}
/*
4 4
10 -7 0 1
-3 2.099999 6 2
5 -1 5 -1
2 1 0 2
4 1
8 5.900001 5 1
*/
