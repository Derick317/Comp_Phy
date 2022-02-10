#include<iostream>
#include"../../Linear_Funtions/matrix.h"
using namespace std;

int main()
{
    //等式右端向量
    vec<double> b(4);
    b.changeelem(1,0.23);
    b.changeelem(2,0.32);
    b.changeelem(3,0.33);
    b.changeelem(4,0.31);
    //高斯消元法
    matrix<double> A(4);  //构造一个四维矩阵
    A.changeelem(1,1,0.05);
    A.changeelem(1,2,0.07);
    A.changeelem(1,3,0.06);
    A.changeelem(1,4,0.05);
    A.changeelem(2,1,0.07);
    A.changeelem(2,2,0.10);
    A.changeelem(2,3,0.08);
    A.changeelem(2,4,0.07);
    A.changeelem(3,1,0.06);
    A.changeelem(3,2,0.08);
    A.changeelem(3,3,0.10);
    A.changeelem(3,4,0.09);
    A.changeelem(4,1,0.05);
    A.changeelem(4,2,0.07);
    A.changeelem(4,3,0.09);
    A.changeelem(4,4,0.10);
    (A.linear_f(b)).print(20,17);
    //Cholesky 分解
    sym_band_matrix<double> AS(4,3);  //构造一个 4 维对称矩阵
    AS.changeelem(1,1,0.05);
    AS.changeelem(1,2,0.07);
    AS.changeelem(1,3,0.06);
    AS.changeelem(1,4,0.05);
    AS.changeelem(2,2,0.10);
    AS.changeelem(2,3,0.08);
    AS.changeelem(2,4,0.07);
    AS.changeelem(3,3,0.10);
    AS.changeelem(3,4,0.09);
    AS.changeelem(4,4,0.10);
    (AS.linear_f(b)).print(20,17);

    return 0;
}