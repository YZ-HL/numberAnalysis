#ifndef NUMBERANALYSISTHIRD_GAUSS_H
#define NUMBERANALYSISTHIRD_GAUSS_H

#include <vector>
#include <cmath>
#include "utils.h"

template<typename T>
struct gauss {
    matrix<T> G;
    const long double eps = 1e-8;

    //传入Ax=b，求解线性方程组,返回自由变量的个数，-1代表无解
    std::pair<std::vector<T>, int> eliminationWithColumnPivoting(matrix<T> &A, matrix<T> &b) {
        auto resultG = constructAugmentedMatrix(A, b);
        if (!resultG.second) return {std::vector<T>(), -1};
        matrix<T> G = resultG.first;
        int r = 0, c = 0, n = G.n, m = G.m - 1;
        std::vector<T> xi(n);
        std::vector<bool> free_x(n);
        while (r < n && c < m) {
            int mr = r;
            for (int i = r + 1; i < n; i++) {
                if (fabs(G[i][c]) > fabs(G[mr][c])) {
                    mr = i;
                }
            }
            if (mr != r) {
                for (int j = c; j < m + 1; j++) {
                    std::swap(G[r][j], G[mr][j]);
                }
            }
            if (fabs(G[r][c]) < eps) {
                G[r][c] = 0;
                ++c;
                continue;
            }
            for (int i = r + 1; i < n; i++) {
                if (G[r][c]) {
                    T tmp = G[i][c] / G[r][c];
                    for (int j = c; j < m + 1; j++) {
                        G[i][j] -= G[r][j] * tmp;
                    }
                }
            }
            ++r;
            ++c;
        }
        for (int i = r; i < n; i++) {
            if (std::fabs(G[i][m]) > eps) {
                return {std::vector<T>(), -1};
            }
        }
        if (r < m) {
            for (int i = r - 1; i >= 0; i--) {
                int free_count = 0, k = -1;
                for (int j = 0; j < m; j++) {
                    if (fabs(G[i][j]) > eps && free_x[j]) {
                        ++free_count;
                        k = j;
                    }
                }
                if (free_count > 0) continue;
                T sum = G[i][m];
                for (int j = 0; j < m; j++) {
                    if (j != k) {
                        sum -= G[i][j] * xi[j];
                    }
                    xi[k] = sum / G[i][k];
                    free_x[k] = 0;
                }
            }
            return {xi, m - r};
        }
        for (int i = m - 1; i >= 0; i--) {
            T sum = G[i][m];
            for (int j = i + 1; j < m; j++) {
                sum -= G[i][j] * xi[j];
            }
            xi[i] = sum / G[i][i];
        }
        return {xi, 0};
    }
};

#endif