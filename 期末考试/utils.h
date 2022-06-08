#ifndef UTILS_H
#define UTILS_H

#include <initializer_list>
#include "matrix.h"

//多变量求解最小值
template<typename T>
T extraMin(const std::initializer_list<T> &args) {
    T mx = *args.begin();
    for (auto x: args) {
        if (mx < x) {
            mx = x;
        }
    }
    return mx;
}

//多变量求解最大值
template<typename T>
T extraMax(const std::initializer_list<T> &args) {
    T mx = *args.begin();
    for (auto x: args) {
        if (mx > x) {
            mx = x;
        }
    }
    return mx;
}

//给定A,b，构造增广矩阵
template<typename T>
std::pair<matrix<T>, bool> constructAugmentedMatrix(matrix<T> &A, matrix<T> &b) {
    if (A.n != b.n || extraMin<int>({A.n, A.m, b.n, b.m}) <= 0 || b.m > 1) {
        return {matrix<T>(), false};
    }
    int level = 0;
    size_t n = A.n, m = A.m + 1;
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
    size_t cnt = 0;
    for (auto x: result) {
        std::cout << "X" << (++cnt) << "=" << x << '\n';
    }
}

template<typename T>
T Tabs(T value) {
    return value >= 0 ? value : -value;
}

template<typename T>
void printXi(std::vector<T> &result, int PC) {
    size_t cnt = 0;
    for (auto x: result) {
        std::cout << std::fixed << std::setprecision(PC) << "X" << (++cnt) << " = " << x << '\n';
    }
}

template<class R, class C, class... Args>
constexpr size_t getFunctionArgsNum(R(C::*)(Args...)) {
    return sizeof...(Args);
}

#endif
