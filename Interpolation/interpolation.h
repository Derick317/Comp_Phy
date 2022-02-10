#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_
#include "../Linear_Funtions/matrix.h" //<iostream>, <cmath>, <complex>, <stdexcept>
#include <vector>

//Cubic Spline, expressed by second derivative
template <class T>
class cubic_spline
{
private:
    int nodenum;      // number of nodes
    std::vector<T> m; // the second derivative of the interpolation function at the nodes
    std::vector<T> x; // the abscissa of the interpolation node
    std::vector<T> y; // the ordinate of the interpolation node

    int bisearch(T somex) const; //Binary search, returns the stored index of the left endpoint of the interval where somex is located (the index starts with 0)

public:
    cubic_spline() : nodenum(0), m(0), x(0), y(0){};                                                          //默认构造函数
    cubic_spline(const T *xp, const T *yp, int n, char c, T cond1 = 0.0, T cond2 = 0.0, T zerolim = ZREOLIM); //构造函数

    T f(T somex) const;               // get the value of the spline function at some x
    T getcoeff(int n, int deg) const; //
};

template <class T>
cubic_spline<T>::cubic_spline(const T *xp, const T *yp, int n, char c, T cond1, T cond2, T zerolim) : nodenum(n), x(xp, xp + n), y(yp, yp + n), m(n) //构造函数
{
    if (n < 2)
        throw std::runtime_error("The number of point is not enough!");
    for (int i = 1; i < n; i++) //check if x is incrementing
        if (xp[i] <= xp[i - 1] + zerolim)
            throw std::runtime_error("Array x is not monotonic increasing!");
    if (c == '1') //1st Boundary Condition
    {
        tridiagmatrix<T> A(n);          //an n-dimensional tridiagonal matrix
        vec<T> d(n);                    //an n-dimensional vector
        for (int i = 1; i < n - 1; i++) //change the element of the matrix and the vector
        {
            T l = (xp[i + 1] - xp[i]) / (xp[i + 1] - xp[i - 1]);
            A.changeelem(i, i, 2);
            A.changeelem(i + 1, i + 2, l);
            A.changeelem(i + 1, i, 1 - l);
            d.changeelem(i + 1, 6 * ((yp[i + 1] - yp[i]) / (xp[i + 1] - xp[i]) - (yp[i] - yp[i - 1]) / (xp[i] - xp[i - 1])) / (xp[i + 1] - xp[i - 1]));
        }
        A.changeelem(n - 1, n - 1, 2);
        A.changeelem(n, n, 2);
        A.changeelem(1, 2, 1);
        A.changeelem(n, n - 1, 1);
        d.changeelem(1, 6 / (xp[1] - xp[0]) * ((yp[1] - yp[0]) / (xp[1] - xp[0]) - cond1));
        d.changeelem(n, 6 / (xp[n - 1] - xp[n - 2]) * (cond2 - (yp[n - 1] - yp[n - 2]) / (xp[n - 1] - xp[n - 2])));
        vec<T> x(A.linear_f(d));
        for (int i = 0; i < n; i++)
            m[i] = x.getelem(i + 1);
    }
    else if (c == '2' || c == 'p' || c == 'P') //2nd or Periodic Boundary Condition
    {
        tridiagmatrix<T> A(n - 2);      //an (n-2)-dimensional tridiagonal matrix
        vec<T> d(n - 2);                //an (n-2)-dimensional vector
        for (int i = 2; i < n - 2; i++) //change the element of the matrix and the vector
        {
            T l = (xp[i + 1] - xp[i]) / (xp[i + 1] - xp[i - 1]);
            A.changeelem(i, i, 2);
            A.changeelem(i, i + 1, l);
            A.changeelem(i, i - 1, 1 - l);
            d.changeelem(i, 6 * ((yp[i + 1] - yp[i]) / (xp[i + 1] - xp[i]) - (yp[i] - yp[i - 1]) / (xp[i] - xp[i - 1])) / (xp[i + 1] - xp[i - 1]));
        }
        A.changeelem(1, 1, 2);
        if (n > 3)
        {
            A.changeelem(n - 2, n - 2, 2);
            A.changeelem(1, 2, (xp[2] - xp[1]) / (xp[2] - xp[0]));
            A.changeelem(n - 2, n - 3, (xp[n - 2] - xp[n - 3]) / (xp[n - 1] - xp[n - 3]));
        }
        d.changeelem(1, (6 * ((yp[2] - yp[1]) / (xp[2] - xp[1]) - (yp[1] - yp[0]) / (xp[1] - xp[0])) - cond1 * (xp[1] - xp[0])) / (xp[2] - xp[0]));
        if (n > 3)
            d.changeelem(n - 2, (6 * ((yp[n - 1] - yp[n - 2]) / (xp[n - 1] - xp[n - 2]) - (yp[n - 2] - yp[n - 3]) / (xp[n - 2] - xp[n - 3])) - cond2 * (xp[n - 1] - xp[n - 2])) / (xp[n - 1] - xp[n - 3]));
        // std::cout<<"-------------- A -------------------"<<std::endl;
        // A.print(6);
        // std::cout<<"-------------- d -------------------"<<std::endl;
        // d.print(6);
        if (c == '2') //2nd Boundary Condition
        {
            vec<T> x(A.linear_f(d));
            m[0] = cond1;
            m[n - 1] = cond2;
            for (int i = 1; i < n - 1; i++)
                m[i] = x.getelem(i);
        }
        else //Periodic Boundary Condition
        {
            lowtrimatrix<T> lm;
            uptrimatrix<T> um;
            vec<T> more(n - 2);
            A.LU(lm, um, zerolim);
            more.changeelem(1, (xp[1] - xp[0]) / (xp[2] - xp[0]));
            more.changeelem(n - 2, (xp[n - 1] - xp[n - 2]) / (xp[n - 1] - xp[n - 3]));
            vec<T> xd(um.linear_f(lm.linear_f(d)));
            vec<T> xmore(um.linear_f(lm.linear_f(more)));
            m[n - 1] = m[0] = (6 * ((yp[1] - yp[0]) / (xp[1] - xp[0]) - (yp[n - 1] - yp[n - 2]) / (xp[n - 1] - xp[n - 2])) - (xp[1] - xp[0]) * xd.getelem(1) - (xp[n - 1] - xp[n - 2]) * xd.getelem(n - 2)) / (2 * (xp[1] - xp[0] + xp[n - 1] - xp[n - 2]) - (xp[1] - xp[0]) * xmore.getelem(1) - (xp[n - 1] - xp[n - 2]) * xmore.getelem(n - 2));
            for (int i = 1; i < n - 1; i++)
                m[i] = xd.getelem(i) - m[0] * xmore.getelem(i);
        }
    }
    else
        throw std::runtime_error("The boundary condition does not exist!");
}

template <class T>
int cubic_spline<T>::bisearch(T somex) const // Binary search
{
    if (somex < x[0] || somex > x[nodenum - 1])
        throw std::range_error("The value is not in the domain of definition!");
    int begin = 0, end = nodenum - 1;
    while (begin < end - 1)
    {
        T middle = (begin + end) / 2;
        if (somex < x[middle])
            end = middle;
        else
            begin = middle;
    }
    return begin;
}

template <class T>
T cubic_spline<T>::f(T somex) const
{
    int lindex = bisearch(somex);
    T hj = x[lindex + 1] - x[lindex], B = y[lindex] - m[lindex] * hj * hj / 6;
    T A = (y[lindex + 1] - y[lindex]) / hj - hj / 6 * (m[lindex + 1] - m[lindex]);
    return m[lindex] / 6 / hj * (x[lindex + 1] - somex) * (x[lindex + 1] - somex) * (x[lindex + 1] - somex) + (m[lindex + 1] / 6 / hj * (somex - x[lindex]) * (somex - x[lindex]) + A) * (somex - x[lindex]) + B;
}

template <class T>
T cubic_spline<T>::getcoeff(int n, int deg) const
{
    if (n <= 0 || n >= nodenum)
        throw std::range_error("The range does not exist!");
    if (deg < 0 || deg > 3)
        throw std::range_error("The degree of the polynomial is not correct!");
    n--;
    T hj = x[n + 1] - x[n], B = y[n] - m[n] * hj * hj / 6;
    T A = (y[n + 1] - y[n]) / hj - hj / 6 * (m[n + 1] - m[n]);
    if (deg == 3) // Get the coefficients of the cubic term of this polynomial
        return (m[n + 1] - m[n]) / 6 / hj;
    if (deg == 2)
        return (m[n] * x[n + 1] - m[n + 1] * x[n]) / 2 / hj;
    if (deg == 1)
        return (-m[n] * x[n + 1] * x[n + 1] + m[n + 1] * x[n] * x[n]) / 2 / hj + A;
    else //deg==1
        return m[n] / 6 / hj * x[n + 1] * x[n + 1] * x[n + 1] + -(m[n + 1] / 6 / hj * x[n] * x[n] + A) * x[n] + B;
}

#endif //INTERPOLATION_H_