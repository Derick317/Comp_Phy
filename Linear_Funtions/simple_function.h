#ifndef SIMPLE_FUNCTION_H_
#define SIMPLE_FUNCTION_H_
#include <complex>


template <class T> inline
std::complex<T> getconj(std::complex<T> z)
{
    return std::conj(z);
}

template <class T> inline
T getconj(T x)
{
    return x;
}


#endif  //SIMPLE_FUNCTION_H_