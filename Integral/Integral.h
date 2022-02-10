#ifndef INTEGRAL_H_
#define INTEGRAL_H_
#include <iostream>
#include <cmath>
#include <stdexcept>

//梯形法则
template <class T>
double trapezoid(T (*f)(T), T left, T right, int nodenum) //nodenum 包含左右端点
{
    if (nodenum < 2)
        throw std::range_error("The number of node is not enough!");
    T step = (right - left) / (nodenum - 1);
    T sum = 0;
    for (int i = 1; i < nodenum - 1; i++)
        sum += f(left + i * step);
    sum += (f(left) + f(right)) / 2;
    return sum * step;
}

// Simpson 法则
template <class T>
double Simpson(T (*f)(T), T left, T right, int nodenum) //nodenum 包含左右端点
{
    if (nodenum < 2)
        throw std::range_error("The number of node is not enough!");
    T step = (right - left) / (nodenum - 1);
    T sum1 = 0, sum2 = 0;
    for (int i = 1; i < nodenum - 1; i++)
        sum1 += f(left + i * step);
    for (int i = 1; i < nodenum; i++)
        sum2 += f(left + i * step - step / 2);
    return (sum1 / 3 + 2 * sum2 / 3 + (f(left) + f(right)) / 6) * step;
}

#endif //end INTEGRAL_H_