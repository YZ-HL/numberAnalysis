#include "matrix.h"
#include "gauss.h"
#include "decomposition.h"

matrix<long double> A, b;

template<typename T>
void test_LUD() {
    //����Ax=b�е�����b
    auto LU_SolveEquations = [&](LUD<T> &GA, matrix<T> &b) {
        std::cout << "��ⷽ���飬b����Ϊ��\n";
        b.printData();
        auto result = GA.solveLinearEquationsByLU(b);
        if (result.second) {
            std::cout << "�����XiΪ��\n";
            printXi(result.first, 15);
            return true;
        } else {
            std::cout << "δ����LU�ֽ⣬�����b������ά����ƥ�䡣\n";
            return false;
        }
    };
    LUD<T> GA = LUD<T>(A);
    auto resultLU = GA.tryIsCanLU();
    if (resultLU.second) {
        std::cout << "���ͨ����˳������ʽΪ��\n";
        int count = 0;
        for (auto x: resultLU.first) {
            std::cout << std::fixed << std::setprecision(10) << "Det " << (++count) << ": " << x << '\n';
        }
        std::cout << "���þ�����Ա�LU�ֽ⣬�ֽ���Ϊ��\n";
        std::cout << "����L��\n";
        GA.L.printData();
        std::cout << "����U��\n";
        GA.U.printData();
        if (LU_SolveEquations(GA, b)) {
            std::cout << "�����������������ʾ��\n";
        } else {
            std::cout << "�����쳣����������b�����ݺ������г���\n";
        }
    } else {
        std::cout << "���ݷǷ��������˳������ʽΪ0������������Է����������LU�ֽ�\n";
    }
}

template<typename T>
void test_ColPivot() {
    gauss<T> gaussTest;
    std::vector<std::pair<int, int>> swapStep;
    auto res = gaussTest.eliminationWithColumnPivoting(A, b, swapStep);
    std::cout << "���ɱ��������� " << res.second << '\n';
    if (!res.second) {
        printXi(res.first, 15);
        std::cout << "�н��������У�\n";
        for (auto row: swapStep) {
            std::cout << "Swap: " << row.first << " " << row.second << '\n';
        }
    } else if (res.second > 0) {
        std::cout << "���������������⣬�����ɱ�������Ϊ��" << res.second << '\n';
    } else {
        std::cout << "�����������󣬻򷽳����޽�\n";
    }
}

int main() {
    A.readMatrix();
    b.readMatrix();
    std::cout << "�����Ƿֽ���⡿\n";
    test_LUD<long double>();
    std::cout << "������Ԫ��ȥ����⡿\n";
    test_ColPivot<long double>();
    system("pause");
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

11 11
68 81 94 107 120 1 14 27 40 53 66
80 93 106 119 11 13 26 39 52 65 67
92 105 118 10 12 25 38 51 64 77 79
104 117 9 22 24 37 50 63 76 78 91
116 8 21 23 36 49 62 75 88 90 103
7 20 33 35 48 61 74 87 89 102 115
19 32 34 47 60 73 86 99 101 114 6
31 44 46 59 72 85 98 100 113 5 18
43 45 58 71 84 97 110 112 4 17 30
55 57 70 83 96 109 111 3 16 29 42
56 69 82 95 108 121 2 15 28 41 54
11 1
671 671 671 671 671 671 671 671 671 671 671

9 9
47 58 69 80 1 12 23 34 45
57 68 79 9 11 22 33 44 46
67 78 8 10 21 32 43 54 56
77 7 18 20 31 42 53 55 66
6 17 19 30 41 52 63 65 76
16 27 29 40 51 62 64 75 5
26 28 39 50 61 72 74 4 15
36 38 49 60 71 73 3 14 25
37 48 59 70 81 2 13 24 35
9 1
369 369 369 369 369 369 369 369 369
*/