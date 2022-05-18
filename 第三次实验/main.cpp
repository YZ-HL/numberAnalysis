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
        for (auto x: resultLU.first){
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
        std::cout << "���ݷǷ��������˳������ʽΪ0�����ܽ���LU�ֽ�\n";
    }
}

template<typename T>
void test_ColPivot() {
    gauss<T> gaussTest;
    auto res = gaussTest.eliminationWithColumnPivoting(A, b);
    std::cout << "���ɱ��������� " << res.second << '\n';
    if (!res.second) {
        printXi(res.first, 15);
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