//用高斯消元法算法解以 Hilbert 矩阵为系数的方程
#include <iostream>
#include "../../Linear_Funtions/matrix.h"
using namespace std;

int main()
{
    for (int n = 2; n <= 15; n++)
    {
        vec<double> b(n);            //一个 n 维向量
        for (int i = 1; i <= n; i++) //修改向量 b 的元素
            b.changeelem(i, 1.0);
        // //高斯消元法
        // matrix<double> H_gem(n);     ////一个 n 维普通矩阵
        // for (int i = 1; i <= n; i++) //给矩阵元赋值
        //     for (int j = 1; j <= n; j++)
        //         H_gem.changeelem(i, j, 1.0 / (i + j - 1));
        // Cholesky 分解法
        sym_band_matrix<double> H(n, n - 1); //一个 n 维对称矩阵
        for (int i = 1; i <= n; i++)         //给矩阵元赋值
            for (int j = 1; j <= i; j++)
                H.changeelem(i, j, 1.0 / (n + n - i - j + 1));
        // cout << "-- n=" << n << " ---GEM---" << endl;
        // (H_gem.linear_f(b)).print(10, 13);
        cout << "-- n=" << n << " ---Cholesky---" << endl;
        (H.linear_f(b)).print(10, 10);
    }

    return 0;
}