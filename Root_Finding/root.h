#ifndef ROOT_H_
#define ROOT_H_
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#define ZEROLIM 0.000000000000001

//二分法求零点
template <class T>
T Bisection(T (*f)(T), T left, T right, T zerolim = ZEROLIM)
{
    T f_l = f(left);
    T f_r = f(right);
    if (f_l * f_r > 0)
        throw std::runtime_error("The two function values have the same sign!");
    T m = (left + right) / 2;
    T f_m = f(m);
    while (right - left > zerolim)
    {
        //std::cout<<std::setprecision(8)<<left<<' '<<right<<std::endl;
        if (f_l * f_m <= 0)
        {
            f_r = f_m;
            right = m;
        }
        else
        {
            f_l = f_m;
            left = m;
        }
        m = (left + right) / 2;
        f_m = f(m);
    }
    return m;
}

//割线法求零点
template <class T>
T Secant(T (*f)(T), T x_0, T x_1, int steplim, T zerolim = ZEROLIM)
{
    T x_next, x_now = x_1, x_last = x_0;
    T f_now = f(x_1), f_last = f(x_0);
    int stepnum = 0;
    for (stepnum = 1; stepnum <= steplim; stepnum++)
    {
        std::cout << std::setprecision(8) << x_now << std::endl;
        if (std::abs(f_now) < zerolim)
            break;
        x_next = x_now - f_now / (f_now - f_last) * (x_now - x_last);
        x_last = x_now;
        x_now = x_next;
        f_last = f_now;
        f_now = f(x_now);
    }
    return x_now;
}

//牛顿法求根
template <class T>
T Newton(T (*f)(T), T (*df)(T), T x_0, int stepnum, T zerolim = ZEROLIM)
{
    for (int i = 0; i < stepnum; i++)
    {
        T f_0 = f(x_0);
        std::cout << std::setprecision(8) << x_0 << std::endl;
        if (std::abs(f_0) < zerolim)
            break;
        x_0 -= f_0 / df(x_0);
    }
    return x_0;
}

template <class T>
T Dekker(T (*f)(T), T a, T b, int stepnum, T zerolim = ZEROLIM)
{
    T f_a = f(a); //当前步的 f(a)
    T f_b = f(b); //当前步的 f(b)
    if (f_a * f_b > 0)
    {
        std::cout << "f(" << a << ") 和 f(" << b << ") 同号" << std::endl;
        return 0;
    }
    T b_last = a;                     //上一步的 b
    T f_b_last = f_a;                 //上一步的 f(b)
    for (int k = 0; k < stepnum; k++) //每一次循环只需要计算一次函数值
    {
        if (std::abs(f_b) > std::abs(f_a)) //如果 |f(b)|>|f(a)|, 则交换 a,b, 使 |f(b)| 更小
        {
            T temp = a;
            a = b;
            b = temp;
            temp = f_a;
            f_a = f_b;
            f_b = temp;
        }
        T s = b - f_b * (b - b_last) / (f_b - f_b_last);
        T m = (a + b) / 2;
        b_last = b;
        f_b_last = f_b;
        if ((b <= s && s <= m) || (m <= s && s <= b))
            b = s;
        else
            b = m;
        f_b = f(b); //更新 f(b)
        if (f_b_last * f_b < 0)
        {
            a = b_last;
            f_a = f_b_last;
        }
        std::cout << k << ":" << b << std::endl;
    }
    return b;
}

//逆二次插值法
template <class T>
T Inverse_Quadratic_Interpolation(const T &a, const T &b, const T &c, const T &f_a, const T &f_b, const T &f_c)
{
    T s = a * f_b / (f_a - f_b) * f_c / (f_a - f_c);
    s += b * f_a / (f_b - f_a) * f_c / (f_b - f_c);
    s += c * f_a / (f_c - f_a) * f_b / (f_c - f_b);
    return s;
}

//s 在不在 a,b 二者之间
template <class T>
bool IsBetween(const T &s, const T &a, const T &b)
{
    if (s < a && s < b)
        return false;
    if (s > a && s > b)
        return false;
    return true;
}

template <class T>
T Brent(T (*f)(T), T a, T b, int stepnum, T delta, T zerolim = ZEROLIM)
{
    T f_a = f(a); //当前步的 f(a)
    T f_b = f(b); //当前步的 f(b)
    if (f_a * f_b > 0)
    {
        std::cout << "f(" << a << ") 和 f(" << b << ") 同号" << std::endl;
        return 0;
    }
    bool IsBisect = true;
    T b_last = a;
    T f_b_last = f_a; //上一步的 f(b)
    T b_last_last = b_last;
    T s, f_s;
    for (int k = 0; k < stepnum; k++)
    {
        if (abs(a - b) < zerolim)
            break;
        if (abs(f_a) < abs(f_b))
        {
            swap(a, b);
            swap(f_a, f_b);
        }
        if (std::abs(f_a - f_b_last) >= zerolim && std::abs(f_b - f_b_last) >= zerolim)
            s = Inverse_Quadratic_Interpolation(a, b, b_last, f_a, f_b, f_b_last);
        else
            s = b - f_b / (f_b - f_a) * (b - a);
        if (!IsBetween(s, b, 0.75 * a + 0.25 * b) || (IsBisect && (std::abs(b - b_last) < delta || std::abs(s - b) >= std::abs(b - b_last) / 2)) || (!IsBisect && (std::abs(b_last - b_last_last) < delta || std::abs(s - b) >= std::abs(b_last - b_last_last) / 2)))
        {
            s = a / 2 + b / 2;
            IsBisect = true;
        }
        else
            IsBisect = false;
        b_last_last = b_last;
        b_last = b;
        f_b_last = f_b;
        f_s = f(s);
        if (f_a * f_s < 0)
        {
            b = s;
            f_b = f_s;
        }
        else
        {
            a = s;
            f_a = f_s;
        }
        std::cout << k + 1 << ":" << b << std::endl;
    }
    return b;
}

template <class T>
bool Iterative_Root(T (*phi)(T), T &x, int stepnum, T (*f)(T) = NULL, T zerolim = ZEROLIM)
{
    T phi_x = phi(x);
    for (int k = 0; k < stepnum; k++)
    {
        if ((!f || std::abs(f(x)) < zerolim) && std::abs(phi_x - x) < zerolim)
            return true;
        x = phi_x;
        phi_x = phi(x);
        std::cout << k + 1 << ":" << x << std::endl;
    }
    return false;
}

template <class T>
bool Iterative_Aitken(T (*phi)(T), T &x, int stepnum, T (*f)(T) = NULL, T zerolim = ZEROLIM)
{
    T x_last = x;
    T phi_x = phi(x);
    T xi_last = x;
    T xi = x;
    for (int k = 0; k < stepnum; k++)
    {
        if ((!f || std::abs(f(xi)) < zerolim) && std::abs(xi - xi_last) < zerolim)
        {
            x = xi;
            return true;
        }
        x_last = x;
        x = phi_x;
        phi_x = phi(x);
        if (k > 1)
        {
            xi_last = xi;
            xi = x_last - (x - x_last) / (phi_x - x - x + x_last) * (x - x_last);
        }
        std::cout << k + 1 << ": x=" << std::setprecision(5) << x << " xi=" << xi << std::endl;
    }
    return false;
}

#endif //ROOT_H_