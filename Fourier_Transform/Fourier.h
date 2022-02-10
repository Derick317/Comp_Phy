#ifndef FOURIER_H_
#define FOURIER_H_
#include <complex> //相当于已经 include<cmath>
#define PI 3.1415926535897932

//快速傅里叶变换, 默认数组长度 n 是 2 的整数次幂; 如果不是的话, 可能会出现神奇的 bug;
template <class T>                                                                            //我们直接认为输入的是复数
void fft(std::complex<T> *x, std::complex<T> *z, std::complex<T> *temp, int n, int delta = 1) //delta表示间隔, 初始值(默认值为 1)
{
    if (n == 1)
    {
        *z = *x;
        return;
    }
    fft(x, z, temp, n / 2, 2 * delta);
    fft(x + delta, z + delta, temp + delta, n / 2, 2 * delta);
    for (int i = 0; i < n; i += 2) //把所有的 z 挪到 temp 中, 并把奇数项乘以单位根的相应幂次
    {
        temp[i * delta] = z[i * delta];
        temp[(i + 1) * delta] = z[(i + 1) * delta] * std::complex<T>(std::cos(i * PI / n), -std::sin(i * PI / n));
    }
    for (int i = 0; i < n / 2; i++)
    {
        z[i * delta] = temp[2 * i * delta] + temp[(2 * i + 1) * delta];
        z[(i + n / 2) * delta] = temp[2 * i * delta] - temp[(2 * i + 1) * delta];
    }
}

#endif //FOURIER_H_