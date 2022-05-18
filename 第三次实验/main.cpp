#include "matrix.h"
#include "gauss.h"
#include "decomposition.h"


matrix<long double> A, b;

template<typename T>
void test_LUD() {
    //传入Ax=b中的向量b
    auto LU_SolveEquations = [&](LUD<T> &GA, matrix<T> &b) {
        std::cout << "求解方程组，b矩阵为：\n";
        b.printData();
        auto result = GA.solveLinearEquationsByLU(b);
        if (result.second) {
            std::cout << "求解结果Xi为：\n";
            printXi(result.first, 15);
            return true;
        } else {
            std::cout << "未进行LU分解，或给定b方程组维数不匹配。\n";
            return false;
        }
    };
    LUD<T> GA = LUD<T>(A);
    auto resultLU = GA.tryIsCanLU();
    if (resultLU.second) {
        std::cout << "检查通过，顺序主子式为：\n";
        int count = 0;
        for (auto x: resultLU.first){
            std::cout << std::fixed << std::setprecision(10) << "Det " << (++count) << ": " << x << '\n';
        }
        std::cout << "即该矩阵可以被LU分解，分解结果为：\n";
        std::cout << "矩阵L：\n";
        GA.L.printData();
        std::cout << "矩阵U：\n";
        GA.U.printData();
        if (LU_SolveEquations(GA, b)) {
            std::cout << "求解结束，结果如上所示。\n";
        } else {
            std::cout << "出现异常，请检查向量b的数据后再运行程序\n";
        }
    } else {
        std::cout << "数据非法；或存在顺序主子式为0，不能进行LU分解\n";
    }
}

template<typename T>
void test_ColPivot() {
    gauss<T> gaussTest;
    auto res = gaussTest.eliminationWithColumnPivoting(A, b);
    std::cout << "自由变量个数： " << res.second << '\n';
    if (!res.second) {
        printXi(res.first, 15);
    } else if (res.second > 0) {
        std::cout << "方程组有无数个解，其自由变量个数为：" << res.second << '\n';
    } else {
        std::cout << "输入数据有误，或方程组无解\n";
    }
}

int main() {
    A.readMatrix();
    b.readMatrix();
    std::cout << "【三角分解求解】\n";
    test_LUD<long double>();
    std::cout << "【列主元消去法求解】\n";
    test_ColPivot<long double>();
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

3 3
3.01 6.03 1.99
1.27 4.16 -1.23
0.987 -4.81 9.34
3 1
1 1 1

3 3
3.00 6.03 1.99
1.27 4.16 -1.23
0.990 -4.81 9.34
3 1
1 1 1
*/