#ifndef NUMBERANALYSISTHIRD_UTILS_H
#define NUMBERANALYSISTHIRD_UTILS_H

#include <initializer_list>
#include "matrix.h"

template<typename T>
//����������Сֵ
T extraMin(const std::initializer_list<T> &args) {
    T mx = *args.begin();
    for (auto x: args) {
        if (mx < x) {
            mx = x;
        }
    }
    return mx;
}

//����A,b�������������
template<typename T>
std::pair<matrix<T>, bool> constructAugmentedMatrix(matrix<T> &A, matrix<T> &b) {
    if (A.n != b.n || extraMin<int>({A.n, A.m, b.n, b.m}) <= 0 || b.m > 1) {
        return {matrix<T>(), false};
    }
    int level = 0, n = A.n, m = A.m + 1;
    matrix<T> result;
    result.n = n, result.m = m;
    for (auto xi: A.G) {
        result.G.push_back(std::vector<T>());
        for (auto yi: xi) {
            (*--result.G.end()).push_back(yi);
        }
        (*--result.G.end()).push_back(b[level++][0]);
    }
    return {result, true};
}

template<typename T>
void printXi(std::vector<T> &result) {
    int cnt = 0;
    for (auto x: result) {
        std::cout << "X" << (++cnt) << "=" << x << '\n';
    }
}

template<typename T>
void printXi(std::vector<T> &result, int PC) {
    int cnt = 0;
    for (auto x: result) {
        std::cout << std::fixed << std::setprecision(PC) << "X" << (++cnt) << " = " << x << '\n';
    }
}

template<class R, class C, class... Args>
constexpr size_t getFunctionArgsNum(R(C::*)(Args...)) {
    return sizeof...(Args);
}

#endif