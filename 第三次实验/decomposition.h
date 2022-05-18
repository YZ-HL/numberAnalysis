#ifndef NUMBERANALYSISTHIRD_DECOMPOSiTION_H
#define NUMBERANALYSISTHIRD_DECOMPOSiTION_H

#include <vector>
#include <iostream>
#include "matrix.h"

template<typename T>
struct LUD {
    matrix<T> G, L, U;
    bool LU_Decomposition_Done = false;

    LUD() {}

    LUD(matrix<T> &R) : G(R) {}

    //检查方阵是否能被LU分解
    std::pair<std::vector<T>, bool> checkIsCanLU(matrix<T> &G) {
        if (!G.isSquareMatrix()) return {std::vector<T>(), false};
        auto leadingPrincipleMinor = G.getLeadingPrincipleMinor();
        if (leadingPrincipleMinor.second == false) return {std::vector<T>(), false};
        for (auto v: leadingPrincipleMinor.first) {
            if (v == 0) return {std::vector<T>(), false};
        }
        return {leadingPrincipleMinor.first, true};
    }

    //尝试LU分解，返回是否成功
    std::pair<std::vector<T>, bool> tryIsCanLU() {
        auto result = checkIsCanLU(G);
        if (result.second) {
            decompositionByLU(G, L, U);
            return result;
        }
        return {std::vector<T>(), false};
    }

    //LU分解，R -> 原矩阵
    void decompositionByLU(matrix<T> &R, matrix<T> &L, matrix<T> &U) {
        int N = R.n;
        L = matrix<T>(N, N);
        U = matrix<T>(N, N);
        for (int i = 0; i < N; i++) U[0][i] = R[0][i], L[i][i] = 1;
        for (int i = 1; i < N; i++) L[i][0] = R[i][0] / U[0][0];
        for (int r = 1; r < N; r++) {
            for (int i = r; i < N; i++) U[r][i] = R[r][i];
            for (int k = 0; k < r; k++) {
                for (int i = r; i < N; i++) {
                    U[r][i] -= L[r][k] * U[k][i];
                }
            }
            if (r + 1 < N) {
                for (int i = r + 1; i < N; i++) L[i][r] = R[i][r];
                for (int k = 0; k < r; k++) {
                    for (int i = r + 1; i < N; i++) {
                        L[i][r] -= L[i][k] * U[k][r];
                    }
                }
                for (int i = r + 1; i < N; i++) {
                    L[i][r] /= U[r][r];
                }
            }
        }
        LU_Decomposition_Done = true;
    }

    //LU分解完成后，求解Ly=b,Ux=y
    std::pair<std::vector<T>, bool> solveLinearEquationsByLU(matrix<T> &b) {
        if (!LU_Decomposition_Done || b.n != G.n || b.m != 1) {
            return {std::vector<T>(), false};
        }
        int Kth = G.n;
        std::vector<T> xi(Kth), yi(Kth);
        for (int i = 0; i < Kth; i++) {
            yi[i] = b[i][0];
            for (int k = 0; k < i; k++) {
                yi[i] -= L[i][k] * yi[k];
            }
        }
        for (int i = Kth - 1; i >= 0; i--) {
            xi[i] = yi[i];
            for (int k = i + 1; k < Kth; k++) {
                xi[i] -= U[i][k] * xi[k];
            }
            xi[i] /= U[i][i];
        }
        return {xi, true};
    }
};

#endif