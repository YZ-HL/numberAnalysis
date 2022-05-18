#include <iostream>
#include <vector>
template<typename T> struct Interpolation{
    std::vector<T> generateEquidistantPoint(T start, T D, int number){
        std::vector<T> result;
        T initVal = start;
        for(int i = 0; i < number; i++, initVal += D){
            result.push_back(initVal);
        }
        return result;
    }
    //手动指定计算用函数。
    std::vector<T> getFunctionResult(std::vector<T> valX){
        std::vector<T> result;
        for(auto x : valX){
            //Runge Function
            result.push_back((T)1.0 / (1.0 + 25 * x * x));
        }
        return result;
    }
    std::string returnNewtonPolynomial(std::vector<T> valX, std::vector<T> helpFunction){
        int siz = valX.size();
        std::string result = "f(x)=";
        for(int i = 0; i < siz; i++){
            result += "(" + std::to_string(helpFunction[i]) + ")";
            for(int j = 0; j < i; j++){
                result = result + "*";
                result = result + "(x - (" + std::to_string(valX[j]) + "))";
            }
            if(i != siz - 1)    result += "+";
        }
        return result;
    }
    std::string returnLagrangePolynomial(std::vector<T> &valX, std::vector<T> &valY, std::vector<T> &omegaS){
        int siz = valX.size();
        std::string result = "f(x)=";
        for(int i = 0; i < siz; i++){
            std::string numerator = "(" + std::to_string(valY[i]) + ")";
            std::string denominator = "(x - (" + std::to_string(valX[i]) + "))";
            for(int j = 0; j < siz; j++){
                numerator = numerator + "*";
                numerator = numerator + "(x - (" + std::to_string(valX[j]) + "))";
            }
            for(int j = 0; j < siz; j++){
                if(j == i)    continue;
                denominator = denominator + "*";
                denominator = denominator + "(" + std::to_string(valX[i] - valX[j]) + ")";
            }
            if(i == 0)
                result = result + "(" + numerator + ")" + "/" + "(" + denominator + ")";
            else
                result = result + "+" + "(" + numerator + ")" + "/" + "(" + denominator + ")";
        }
        return result;
    }
    std::string newtonInterpolation(std::vector<T> &x, std::vector<T> &fx){
        int sizX = x.size(), sizF = fx.size();
        if(sizX != sizF || !sizX)    return "Error!";
        std::vector<std::vector<T> > diff(sizX, std::vector<T>(sizX, 0));
        for(int i = 0; i < sizX; i++) {
            diff[0][i] = fx[i];
        }
        for(int i = 1; i < sizX; i++){
            for(int j = 0; j < sizX - i; j++){
                diff[i][j] = (diff[i - 1][j + 1] - diff[i - 1][j]) / (x[j + i] - x[j]);
            }
        }
        std::vector<T> useDiff;
        for(int i = 0; i < sizX; i++)    useDiff.push_back(diff[i][0]);
        return returnNewtonPolynomial(x, useDiff);
    }
    std::string lagrangeInterpolation(std::vector<T> &x, std::vector<T> &fx){
        int sizX = x.size(), sizF = fx.size();
        if(sizX != sizF || !sizX)    return "Error!";
        std::vector<T> omegaS;
        for(int i = 0; i < sizX; i++){
            T sum = 1;
            for(int j = 0; j < sizX; j++){
                if(i == j)    continue;
                sum = sum * (x[i] - x[j]);
            }
            omegaS.push_back(sum);
        }
        return returnLagrangePolynomial(x, fx, omegaS);
    }
};
void problem1(){
    Interpolation<long double> solve;
    std::vector<long double> xi = {0.2, 0.4, 0.6, 0.8, 1.0};
    std::vector<long double> fx = {0.98, 0.92, 0.81, 0.64, 0.38};
    auto result = solve.lagrangeInterpolation(xi, fx);
    std::cout << result << '\n';
};
void problem2_n(long double start, long double D, int number){
    Interpolation<long double> solve;
    auto xi = solve.generateEquidistantPoint(start, D, number);
    auto fx = solve.getFunctionResult(xi);
    for(auto x : xi)    std::cout << x << '\n';
    for(auto y : fx)    std::cout << y << '\n';
    auto result = solve.lagrangeInterpolation(xi, fx);
    std::cout << result << '\n';
}
void problem3(){
    Interpolation<long double> solve;
    std::vector<long double> xi = {0, 1, 4, 9, 16, 25, 36, 49, 64};
    std::vector<long double> fx = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    auto result = solve.lagrangeInterpolation(xi, fx);
    std::cout << result << '\n';
}
int main() {
    //problem1();
    //problem2_n(-0.5,0.1,10);
    //problem2_n(-1,0.2,10);
    //problem2_n(-0.5, 0.05, 20);
    //problem3();
    return 0;
}