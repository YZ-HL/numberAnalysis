#include <iostream>
#include "matrix.h"
#include "iterative.h"
#include "utils.h"

template<typename T, typename F>
void testIteration(T start, const F &fun, T eps, T allowMaxValue, std::initializer_list<T> &xList) {
    try {
        auto[result, count] = functionIteration<T>(start, fun, eps, allowMaxValue, xList);
        std::cout << "�������Ϊ: " << result << '\n';
        std::cout << "��������Ϊ: " << count << '\n';
    } catch (std::domain_error &e) {
        std::cout << "Exception caught: " << e.what() << '\n';
    }

}

template<typename T>
void testEx1() {
    static const T zero_eps = 1e-10;
    std::cout << "������������� - 1��\n";
    std::initializer_list<T> fixedPointDomain_1 = {-1e20, -zero_eps, zero_eps, 1e20};
    testIteration((T) 1.5, [&](T x) { return 1 + 1.0 / (x * x); }, (T) 1e-15, (T) 1e20, fixedPointDomain_1);
    std::initializer_list<T> fixedPointDomain_2 = {1, 1e20};
    std::cout << "������������� - 2��\n";
    testIteration((T) 1.5, [&](T x) { return 1.0 / sqrt(x - 1); }, (T) 1e-15, (T) 1e20, fixedPointDomain_2);
    std::cout << "��ţ�ٵ�������\n";
    std::initializer_list<T> newtonDomain = {-1e20, -zero_eps, zero_eps, 2.0 / 3, 2.0 / 3 + zero_eps, 1e20};
    testIteration((T) 1.5, [&](T x) { return x - (x * x * x - x * x - 1) / (3 * x * x - 2 * x); }, (T) 1e-15, (T) 1e20, newtonDomain);
}

template<typename T>
void solveEx2(matrix<T> A, matrix<T> b) {
    auto checkResult = [&](std::vector<T> solutionVector, int stepCount) {
        if (stepCount > 0) {
            std::cout << "��������: " << stepCount << '\n';
            std::cout << "����������: \n";
            printXi(solutionVector, 10);
        } else {
            std::cout << "���������н�ľ���ֵ���������������ָ�������ޣ�����ԭ�򣺵���������������ָ��������ֵ��С��" << '\n';
        }
    };
    auto solveByJacobi = [&](matrix<T> A, matrix<T> b) {
        std::cout << "Jacobi������������£�\n";
        auto[solutionVector, stepCount] = Jacobi<T>(A, b, 1e-4, 1e20, 1e9);
        checkResult(solutionVector, stepCount);
    };
    auto solveBySOR = [&](matrix<T> A, matrix<T> b, T mu) {
        std::cout << "Gauss-Seidel������������£�\n";
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
