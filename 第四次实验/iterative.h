#ifndef ITERATIVE_H
#define ITERATIVE_H

#include <vector>
#include <exception>
#include <cfloat>
#include <cmath>
#include "matrix.h"
#include "utils.h"

template<typename T>
void iterativeBasicCheck(size_t sizA, size_t sizB, T eps, T allowMaxValue) {
    if (sizA != sizB) {
        throw std::runtime_error("The dimensions of the matrix A and the matrix B are not correspond.");
    }
    if (eps < 0) {
        throw std::runtime_error("eps must greater than zero!");
    }
    if (allowMaxValue < 0) {
        throw std::runtime_error("allowMaxValue must greater than zero!");
    }
}

template<typename T>
std::pair<std::vector<T>, size_t> Jacobi(matrix<T> A, matrix<T> b, T eps, T allowMaxValue, size_t allowMaxIterativeCount) {
    size_t sizA = A.n, sizB = b.n;
    iterativeBasicCheck(sizA, sizB, eps, allowMaxValue);

    auto xi_old = b.to_vector();
    auto xi_new = std::vector<T>(sizB);

    for (size_t i = 0; i < sizB; i++) {
        xi_old[i] /= A[i][i];
    }

    T eps_i = DBL_MAX;
    size_t iterative_count = 1;
    while (eps_i > eps) {
        iterative_count++;
        if (iterative_count > allowMaxIterativeCount) {
            return {std::vector<T>(), 0};
        }
        for (size_t i = 0; i < sizB; i++) {
            xi_new[i] = b[i][0];
            for (size_t j = 0; j < sizA; j++) {
                if (i == j) continue;
                xi_new[i] -= A[i][j] * xi_old[j];
            }
            xi_new[i] /= A[i][i];
        }
        eps_i = 0;
        for (size_t i = 0; i < sizB; i++) {
            eps_i = std::max(eps_i, (T) fabs(xi_new[i] - xi_old[i]));
            if ((T)fabs(xi_new[i]) > allowMaxValue) {
                return {std::vector<T>(), 0};
            }
        }
        xi_old = xi_new;
    }
    return {xi_new, iterative_count};
}

template<typename T>
std::pair<std::vector<T>, size_t> SOR(matrix<T> A, matrix<T> b, T mu, T eps, T allowMaxValue, size_t allowMaxIterativeCount) {
    size_t sizA = A.n, sizB = b.n;
    iterativeBasicCheck(sizA, sizB, eps, allowMaxValue);

    auto xi_old = std::vector<T>(sizB, 0);
    auto xi_new = std::vector<T>(sizB, 0);

    T eps_i = DBL_MAX;
    size_t iterative_count = 1;
    while (eps_i > eps) {
        iterative_count++;
        if (iterative_count > allowMaxIterativeCount) {
            return {std::vector<T>(), 0};
        }
        for (size_t i = 0; i < sizB; i++) {
            xi_new[i] = xi_old[i];
            T tmpValue = b[i][0];
            for (size_t j = 0; j < i; j++) {
                tmpValue -= A[i][j] * xi_new[j];
            }
            for (size_t j = i; j < sizB; j++) {
                tmpValue -= A[i][j] * xi_old[j];
            }
            xi_new[i] += mu * tmpValue / A[i][i];
        }
        eps_i = 0;
        for (size_t i = 0; i < sizB; i++) {
            eps_i = std::max(eps_i, (T) fabs(xi_new[i] - xi_old[i]));
            if ((T)fabs(xi_new[i]) > allowMaxValue) {
                return {std::vector<T>(), 0};
            }
        }
        xi_old = xi_new;
    }
    return {xi_new, iterative_count};
}

#endif
