#ifndef NUMBERANALYSISTHIRD_MATRIX_H
#define NUMBERANALYSISTHIRD_MATRIX_H

#include <iostream>
#include <vector>
#include <iomanip>

template<typename T>
struct matrix {
    int n, m;
    std::vector<std::vector<T>> G;

    matrix() {}

    matrix(int _n, int _m) {
        n = _n;
        m = _m;
        G = std::vector<std::vector<T>>(n, std::vector<T>(m, 0));
    }

    std::vector<T> &operator[](int k) { return G[k]; }

    //读入一个n*m的矩阵
    void readMatrix() {
        std::cin >> n >> m;
        *this = matrix<T>(n, m);
        for (int r = 0; r < n; r++) {
            for (int c = 0; c < m; c++) {
                std::cin >> G[r][c];
            }
        }
    }

    //判断当前矩阵是否为合法方阵
    bool isSquareMatrix() {
        if (n != m || std::min(n, m) <= 0) return false;
        return true;
    }

    //默认精度为3
    void printData() {
        for (auto r: G) {
            for (auto c: r) {
                std::cout << std::fixed << std::setprecision(3) << c << ' ';
            }
            std::cout << '\n';
        }
    }

    void printData(int PC) {
        for (auto r: G) {
            for (auto c: r) {
                std::cout << std::fixed << std::setprecision(PC) << c << ' ';
            }
            std::cout << '\n';
        }
    }

    //递归计算行列式G的值
    std::pair<T, bool> getDeterminant(matrix<T> &G) {
        if (!G.isSquareMatrix()) return {0, false};
        int Kth = G.n;
        if (Kth == 1) return {G[0][0], true};
        if (Kth == 2) return {G[0][0] * G[1][1] - G[0][1] * G[1][0], true};
        T result = 0;
        int sgn = (Kth % 2 == 0) ? 1 : -1;
        for (int i = 0; i < Kth; i++) {
            matrix<T> recursionMatrix = matrix<T>();
            recursionMatrix.G.resize(Kth - 1);
            recursionMatrix.n = recursionMatrix.m = Kth - 1;
            for (int j = 0; j < Kth - 1; j++) {
                for (int k = 0; k < Kth; k++) {
                    if (k != i) recursionMatrix[j].push_back(G[j][k]);
                }
            }
            auto Det = getDeterminant(recursionMatrix);
            result += sgn * Det.first * G[Kth - 1][i];
            sgn *= -1;
        }
        return {result, true};
    }

    //求取顺序主子式并返回计算结果
    std::pair<std::vector<T>, bool> getLeadingPrincipleMinor() {
        if (!(*this).isSquareMatrix()) return {std::vector<T>(), false};
        std::vector<T> result;
        for (int i = 0; i < n; i++) {
            matrix<T> recursionMatrix = matrix<T>(i + 1, i + 1);
            for (int r = 0; r <= i; r++) {
                for (int c = 0; c <= i; c++) {
                    recursionMatrix[r][c] = (*this)[r][c];
                }
            }
            auto result_i = getDeterminant(recursionMatrix);
            if (result_i.second == false) return {std::vector<T>(), false};
            result.push_back(result_i.first);
        }
        return {result, true};
    }
};

#endif