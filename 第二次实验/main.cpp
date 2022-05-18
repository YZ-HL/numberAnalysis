#include <bits/stdc++.h>
template<typename T> struct Integration{
    std::vector<T> generateEquidistantPoint(T start, T D, int siz){
        std::vector<T> result;
        T initVal = start;
        for(int i = 0; i < siz; i++, initVal += D){
            result.push_back(initVal);
        }
        return result;
    }
    std::vector<T> generateRandomPoint(T leftPoint, T rightPoint, int siz){
        std::default_random_engine e;
        std::uniform_real_distribution<T> u(leftPoint, rightPoint);
        std::vector<T> result;
        for(int i = 0; i < siz; i++){
            result.push_back(u(e));
        }
        return result;
    }
    template<typename F> std::pair<T, bool> compositeTrapezoidal(const F _function, std::vector<T> &setX){
        if(!std::is_sorted(setX.begin(), setX.end()))    std::sort(setX.begin(), setX.end());
        setX.erase(std::unique(setX.begin(), setX.end()), setX.end());
        int sizX = setX.size();
        if(sizX <= 2)    return {0, false};
        T result = 0, h = (setX.back() - setX.front()) / (sizX - 1);
        for(int i = 1; i < sizX - 1; i++){
            result += 2 * _function(setX[i]);
        }
        result = h / 2 * (_function(setX.front()) + _function(setX.back()) + result);
        return {result, true};
    }
    template<typename F> std::pair<T, bool> compositeSimpson(const F _function, std::vector<T> &setX){
        if(!std::is_sorted(setX.begin(), setX.end()))    std::sort(setX.begin(), setX.end());
        setX.erase(std::unique(setX.begin(), setX.end()), setX.end());
        int sizX = setX.size();
        if(sizX <= 2)    return {0, false};
        T result = 0, h = (setX.back() - setX.front()) / (sizX - 1);
        for(int i = 0; i < sizX - 1; i++){
            result += 4 * _function((setX[i] + setX[i + 1]) / 2);
        }
        for(int i = 1; i < sizX - 1; i++){
            result += 2 * _function(setX[i]);
        }
        result = h / 6 * (_function(setX.front()) + _function(setX.back()) + result);
        return {result, true};
    }
}; using LDI = Integration<long double>;
//#define DEBUG
template<typename T> void printIntegrationResult(std::vector<T> setX, std::pair<T, bool> result){
    if(!result.second){
        std::cout << "ERROR\n";
        return;
    }
#ifdef DEBUG
    int cnt = 0;
    for(auto x : setX){
        std::cout << (++cnt) << " : " << x << '\n';
    }
#endif
    std::cout << "val: " << result.first << '\n';
}
int main(void){
    auto Function1 = [&](long double x){ return sqrt(x) * log(x); };
    auto Function2 = [&](long double x){ return 4.0 / (1 + x * x); };
    LDI testIntegration;
    auto pointSet1 = testIntegration.generateEquidistantPoint(0.001, 0.001, 1000);
    auto pointSet2 = testIntegration.generateRandomPoint(0.001, 1, 1000);
    std::cout << "【问题1】";
    std::cout << "在(0,1]上等距选点得到的结果：" << '\n';
    auto ans11 = testIntegration.compositeTrapezoidal(Function1, pointSet1);
    auto ans12 = testIntegration.compositeSimpson(Function1, pointSet1);
    printIntegrationResult(pointSet1, ans11);
    printIntegrationResult(pointSet2, ans12);

    std::cout << "在(0,1]上随机选点得到的结果：" << '\n';
    auto ans13 = testIntegration.compositeTrapezoidal(Function1, pointSet2);
    auto ans14 = testIntegration.compositeSimpson(Function1, pointSet2);
    printIntegrationResult(pointSet1, ans13);
    printIntegrationResult(pointSet2, ans14);

    std::cout << "【问题2】";
    std::cout << "在(0,1]上等距选点得到的结果：" << '\n';
    auto ans21 = testIntegration.compositeTrapezoidal(Function2, pointSet1);
    auto ans22 = testIntegration.compositeSimpson(Function2, pointSet1);
    printIntegrationResult(pointSet1, ans21);
    printIntegrationResult(pointSet2, ans22);

    std::cout << "在(0,1]上随机选点得到的结果：" << '\n';
    auto ans23 = testIntegration.compositeTrapezoidal(Function2, pointSet2);
    auto ans24 = testIntegration.compositeSimpson(Function2, pointSet2);
    printIntegrationResult(pointSet1, ans23);
    printIntegrationResult(pointSet2, ans24);

    return 0;
}