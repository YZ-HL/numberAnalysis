#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <exception>

template<typename T>
struct matrix {
    size_t n, m;
    std::vector<std::vector<T>> G;

    matrix() {}

    matrix(size_t _n, size_t _m) {
        n = _n;
        m = _m;
        G = std::vector<std::vector<T>>(n, std::vector<T>(m, 0));
    }

    void fillMatrix(T val) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < m; j++) {
                G[i][j] = val;
            }
        }
    }

    std::vector<T> to_vector() {
        if (m != 1) {
            throw std::runtime_error("The dimensions of the matrix do not correspond. Can't convert to vector.");
        }
        std::vector<T> result;
        for (auto x: G) {
            for (auto y: x) {
                result.push_back(y);
            }
        }
        return result;
    }

    T getInfinityNorm() {
        T result = 0;
        for (size_t i = 0; i < m; i++) result += G[0][i];
        for (size_t i = 0; i < n; i++) {
            T nowValue = 0;
            for (size_t j = 0; j < m; j++) {
                nowValue += G[i][j];
            }
            result = max(result, nowValue);
        }
        return result;
    }

    matrix<T> operator*(const matrix &a) {
        if (m != a.n) {
            throw std::runtime_error("The dimensions of the matrix do not correspond. Can't convert to vector.");
        }
        matrix<T> result(n, a.m);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < a.m; j++) {
                for (size_t k = 0; k < m; k++) {
                    result[i][j] += G[i][k] * a.G[k][j];
                }
            }
        }
        return result;
    }

    std::vector<T> &operator[](size_t k) { return G[k]; }

    //读入一个n*m的矩阵
    void readMatrix() {
        std::cin >> n >> m;
        *this = matrix<T>(n, m);
        for (size_t r = 0; r < n; r++) {
            for (size_t c = 0; c < m; c++) {
                std::cin >> G[r][c];
            }
        }
    }

    //判断当前矩阵是否为合法方阵
    bool isSquareMatrix() {
        if (n != m || std::min(n, m) <= 0) {
            return false;
        }
        return true;
    }

    //默认精度为3
    void printData() {
        std::cout << n << " " << m << '\n';
        for (auto r: G) {
            for (auto c: r) {
                std::cout << std::fixed << std::setprecision(3) << c << ' ';
            }
            std::cout << '\n';
        }
    }

    //自定义精度
    void printData(int PC) {
        std::cout << n << " " << m << '\n';
        for (auto r: G) {
            for (auto c: r) {
                std::cout << std::fixed << std::setprecision(PC) << c << ' ';
            }
            std::cout << '\n';
        }
    }
};

#endif