#ifndef MATRIX_H_
#define MATRIX_H_
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <list>
#include <complex>
#include <stdexcept>
#include "simple_function.h"

const double ZREOLIM = 0.0000000000000001; //The default upper limit of 0

// T should not be set to ‘int’, because the division of int type is integer division, and there is no overloaded division operator /
template <class T>
class vec;
template <class T>
class sym_band_matrix; //FORWARD REFERENCE declaration
template <class T>
class lowtrimatrix; //FORWARD REFERENCE declaration
template <class T>
class uptrimatrix; //FORWARD REFERENCE declaration
template <class T>
class diagmatrix;
template <class T>
class sparse_matrix;
template <class T>
class matrix;
template <class Type>
bool cgm(vec<Type> &x, const vec<Type> &b, const sym_band_matrix<Type> &A, int stepnum, Type zerolim = ZREOLIM); //共轭梯度算法

//Nodes of a sparse matrix
template <class T>
struct spmnode
{
    T elem;    //实际存储元素
    int index; //元素真实列号
};

//a vector
template <class T>
class vec
{
    template <class Type>
    friend class tridiagmatrix;
    template <class Type>
    friend class lowtrimatrix;
    template <class Type>
    friend class uptrimatrix;
    template <class Type>
    friend class matrix;
    template <class Type>
    friend class diagmatrix;
    template <class Type>
    friend class sym_band_matrix;
    template <class Type>
    friend void swap(vec<Type> &lhs, vec<Type> &rhs); //交换两个向量
    template <class Type>
    friend Type f(const vec<Type> &v1, const vec<Type> &v2); //向量內积
    template <class Type>
    friend std::complex<Type> f(const vec<std::complex<Type>> &v1, const vec<std::complex<Type>> &v2); //Hermite 內积
    template <class Type>
    friend vec<Type> operator*(const sym_band_matrix<Type> &A, const vec<Type> &v); //矩阵乘向量
    template <class Type>
    friend vec<Type> operator*(const matrix<Type> &A, const vec<Type> &ve); //矩阵乘向量
    template <class Type>
    friend vec<Type> operator*(const sparse_matrix<Type> &spma, const vec<Type> &ve); //稀疏矩阵乘向量
    template <class Type>
    friend Type f(const sym_band_matrix<Type> &A, const vec<Type> &p); //內积 (p^T)Ap
    template <class Type>
    friend bool cgm(vec<Type> &x, const vec<Type> &b, const sym_band_matrix<Type> &A, int stepnum, Type zerolim); //共轭梯度算法
    template <class Type>
    friend vec<Type> operator*(const Type &x, vec<Type> ve); //overload operator向量数乘

private:
    int dim;
    bool IsCol; //1表示列向量，0表示行向量
    T *v;       //指向存储向量的指针
public:
    vec() : dim(0), IsCol(1), v(NULL) {}          //default constructor
    vec(int n, bool c = true);                    //构造函数, 构造 n 维 0 向量
    vec(const vec &ve);                           //copy constructor
    vec(vec &&ve) noexcept;                       //moving constructor
    vec(const matrix<T> &ma, int i, bool is_col); //构造函数，用矩阵第 i 行或第 i列构造
    ~vec() { delete[] v; }                        //Destructor 
    vec &operator=(vec rhs);                      //overload operator = 
    vec &operator+=(const vec &rhs);              //overload operator +=
    vec &operator-=(const vec &rhs);              //overload operator -=
    vec &operator*=(const T &num);                //数乘一个数

    void trans() { IsCol = !IsCol; }                //转置
    int getdim() const { return dim; }              //获得维数
    bool direction() const { return IsCol; }        //列向量还是行向量
    T getelem(int i) const;                         //返回第i个元素
    void changeelem(int i, T elem);                 //将第i个元素的值改为elem
    bool exch_elem(int i, int j);                   //swap two elements
    void print(int w = 4, int precision = 0) const; //打印此向量
};

//An ordinary matrix
template <class T>
class matrix
{
    template <class Type>
    friend class lowtrimatrix;
    template <class Type>
    friend class uptrimatrix;
    template <class Type>
    friend class vec;
    template <class Type>
    friend void swap(matrix<Type> &lhs, matrix<Type> &rhs); //交换两个矩阵
    template <class Type>
    friend vec<Type> operator*(const matrix<Type> &A, const vec<Type> &v); //矩阵乘向量
    template <class Type>
    friend matrix<Type> operator*(const matrix<Type> &ma, const uptrimatrix<Type> &up); //矩阵乘上三角矩阵
    template <class Type>
    friend matrix<Type> operator*(const matrix<Type> &A, const matrix<Type> &B); //矩阵乘矩阵
    template <class Type>
    friend matrix<Type> operator*(const sym_band_matrix<Type> &sm, const matrix<Type> &ma); //对称矩阵乘以普通矩阵
    template <class Type>
    friend matrix<Type> operator*(const matrix<Type> &ma, const sym_band_matrix<Type> &sm); //普通矩阵乘对称矩阵
    template <class Type>
    friend matrix<Type> operator*(const uptrimatrix<Type> &um, const matrix<Type> &ma); //左乘一个上三角矩阵

private:
    int dim; //矩阵级数
    T **m;   //指向储存矩阵的二维数组的指针

public:
    matrix();                     //default constructor
    matrix(int n);                //构造函数, 构造n级0矩阵
    matrix(const matrix &ma);     //copy constructor
    matrix(matrix &&ma) noexcept; //moving constructor
    ~matrix();                    //Destructor 
    matrix &operator=(matrix ma); //overload operator = 

    int getdim() const { return dim; }                                               //获得级数
    T getelem(int i, int j) const;                                                   //返回第i行第j列的值
    void changeelem(int i, int j, T elem);                                           //将第i行第i列的值改为elem
    T row_max_abs(int i, int b, int e) const;                                        //获取第i行第b个(包含)到第e个元素(不包含)中最大绝对值
    int row_maxabs_index(int i, int b, int e) const;                                 //获取第i行第b个(包含)到第e个元素(不包含)中最大绝对值元素的列数, 返回值在1~n之间
    T col_max_abs(int i, int b, int e) const;                                        //获取第i列第b个(包含)到第e个元素(不包含)中最大绝对值
    int col_maxabs_index(int i, int b, int e) const;                                 //获取第i列第b个(包含)到第e个元素(不包含)中最大绝对值元素的行数, 返回值在1~n之间
    void exch_row(int i, int j);                                                     //交换两行元素
    void exch_col(int i, int j);                                                     //交换两列元素
    void proportion_mainelem(int i);                                                 //比例因子主元, 找到第i~n行, 第i~n列的子矩阵中比例因子主元所在行，并把这一行与第i行交换.
    void partial_mainelem(int i);                                                    //部分主元, 找到第i~n行, 第i~n列的子矩阵中部分主元所在行，并把这一行与第i行交换
    bool elementary_transformation(int i1, T k, int i2, bool is_col, int begin = 1); //矩阵 i1 行(或列)的 k 倍加到 i2 行(列), 操作范围是这一行(列)的 begin 到末尾
    void print(int w = 4) const;                                                     //打印此矩阵

    vec<int> LU(int mainelem = 1, T zerolim = ZREOLIM);                //LU 分解，返回一个用向量表示的置换矩阵 P, 使得 PA=LU.
    uptrimatrix<T> Schmidt(T zerolim = ZREOLIM);                       //Schmidt 正交化方法
    uptrimatrix<T> Householder_QR(T zerolim = ZREOLIM);                // Householder方法进行 QR 分解
    uptrimatrix<T> Givens_QR(T zerolim = ZREOLIM);                     // Givens变换进行 QR 分解
    matrix<T> upper_Hessenberg(T zerolim = ZREOLIM);                   //正交相似与一个上海森堡矩阵 (*this=Q^THQ), 返回 Q
    vec<T> linear_f(vec<T> &b, int mainelem = 1, T zerolim = ZREOLIM); //返回以此矩阵为系数矩阵，以b为非齐次项的线性方程组的解, 列主元采用部分主元
};

//a sparse matrix
template <class T>
class sparse_matrix
{
    template <class Type>
    friend vec<Type> operator*(const sparse_matrix<Type> &spma, const vec<Type> &ve); //稀疏矩阵乘向量
    template <class Type>
    friend void swap(sparse_matrix<Type> &lhs, sparse_matrix<Type> &rhs); //交换稀疏两个矩阵

private:
    int dim;
    std::list<spmnode<T>> *spm;

public:
    sparse_matrix() : dim(0), spm(NULL) {}       //default constructor
    sparse_matrix(int n);                        //构造函数
    sparse_matrix(const sparse_matrix &sma);     //copy constructor
    sparse_matrix(sparse_matrix &&sma) noexcept; //moving constructor
    ~sparse_matrix() { delete[] spm; }           //Destructor 
    sparse_matrix &operator=(sparse_matrix sma); //overload operator = 

    int getdim() const { return dim; }     //获得级数
    T getelem(int i, int j) const;         //返回第i行第j列的值
    bool changeelem(int i, int j, T elem); //将第i行第i列的值改为elem
    void print(int w = 6) const;           //打印此稀疏矩阵
};

//symmetric banded matrix
template <class T>
class sym_band_matrix
{
    template <class Type>
    friend class matrix;
    template <class Type>
    friend void swap(sym_band_matrix<Type> &lhs, sym_band_matrix<Type> &rhs); //交换两个对称带状矩阵
    template <class Type>
    friend vec<Type> operator*(const sym_band_matrix<Type> &A, const vec<Type> &v); //矩阵乘向量
    template <class Type>
    friend matrix<Type> operator*(const sym_band_matrix<Type> &sm, const matrix<Type> &ma); //对称矩阵乘以普通矩阵
    template <class Type>
    friend matrix<Type> operator*(const matrix<Type> &ma, const sym_band_matrix<Type> &sm); //普通矩阵乘对称矩阵
    template <class Type>
    friend Type f(const sym_band_matrix<Type> &A, const vec<Type> &p); //內积 (p^T)Ap
    template <class Type>
    friend bool cgm(vec<Type> &x, const vec<Type> &b, const sym_band_matrix<Type> &A, int stepnum, Type zerolim); //共轭梯度算法

protected:
    int dim;             //矩阵维数
    int halfbandwidth;   //矩阵半带宽, 0矩阵用 halfbandwidth=-1 表示, 对角矩阵 halfbandwidth=0
    T **sbm;             //指向存储矩阵数组的指针，只存储左下一半, 第 i 行 j 列(i>=j)存在sbm[i-j][j-1]
    T **copysbm() const; //复制一个和sbm相同的二维数组，并返回其头指针
    void dltsbm();       //把数组 sbm 释放掉

public:
    sym_band_matrix() : dim(0), halfbandwidth(-1), sbm(NULL) {} //default constructor, 构造 0 级矩阵, 即什么都没有
    sym_band_matrix(int n, int hbw);                            //构造函数, 构造 n 级,半带宽 hbw 的 0 矩阵
    sym_band_matrix(const sym_band_matrix &sbma);               //copy constructor
    sym_band_matrix(sym_band_matrix &&sbma) noexcept;           //moving constructor
    ~sym_band_matrix();                                         //Destructor 
    sym_band_matrix &operator=(sym_band_matrix sbma);           //overload operator "=" 

    int getdim() const { return dim; }                     //获得级数
    int gethalfbandwidth() const { return halfbandwidth; } //获取半带宽
    T getelem(int i, int j) const;                         //返回第i行第j列的值
    void changeelem(int i, int j, T elem);                 //change the element in row i, column j into elem(同时也会改变对称元素)
    //bool elementary_transformation(int i1, T k, int i2, bool is_col, int begin = 1); //矩阵 i1 行(或列)的 k 倍加到 i2 行(列), 操作范围是这一行(列)的 begin 到末尾
    void print(int w = 4) const; //打印此对称带状矩阵

    bool Cholesky(lowtrimatrix<T> &lm, diagmatrix<T> &dm, T zerolim = ZREOLIM); // Cholesky 分解
    matrix<T> Jacobi(int i, int j, T zerolim = ZREOLIM);                        //一次 Jacobi 迭代, (*this)=QᵗAQ, 返回 Q, (*this)变为 A
    T Sturm(int n, int stepnum, T zerolim = ZREOLIM);                           // Sturm 方法求从小到大第 n 个本征值
    vec<T> linear_f(const vec<T> &b, T zerolim = ZREOLIM);                      //返回以此上三角矩阵为系数矩阵，以 b为非齐次项的线性方程组的解
};

//tridiagonal matrix
template <class T>
class tridiagmatrix
{
private:
    int dim; //矩阵维数
    T **tdm; //指向存储矩阵数组的指针, 第 i 行第 j 列存在 tdm[j-i+1][j-1]

public:
    tridiagmatrix() : dim(0), tdm(NULL) {}        //default constructor, 构造 0 级矩阵, 即什么都没有
    tridiagmatrix(int n);                         //构造函数, 构造 n 级 0 矩阵
    tridiagmatrix(const tridiagmatrix &tdma);     //copy constructor
    tridiagmatrix(tridiagmatrix &&tdma) noexcept; //moving constructor
    ~tridiagmatrix();                             //Destructor 
    tridiagmatrix &operator=(tridiagmatrix tdma); //overload operator "=" 

    int getdim() const { return dim; }     //获得级数
    T getelem(int i, int j) const;         //返回第i行第j列的值
    void changeelem(int i, int j, T elem); //change the element in row i, column j into elem
    void print(int w = 4) const;           //打印此三对角矩阵

    bool LU(lowtrimatrix<T> &lm, uptrimatrix<T> &um, T zerolim = ZREOLIM) const; //三对角矩阵 LU 分解
    vec<T> linear_f(const vec<T> &b, T zerolim = ZREOLIM) const;                 //返回以此三对角矩阵为系数矩阵，以 b 为非齐次项的线性方程组的解
};

//diagonal matrix, inheret from sym_band_matrix. when dim=0, hbw=-1; otherwise, hbw=0
template <class T>
class diagmatrix : public sym_band_matrix<T>
{
    template <class Type>
    friend uptrimatrix<Type> operator*(uptrimatrix<Type> um, const diagmatrix<Type> &dm); //上三角矩阵乘对角矩阵
    template <class Type>
    friend lowtrimatrix<Type> operator*(lowtrimatrix<Type> lm, const diagmatrix<Type> &dm); //下三角矩阵乘对角矩阵

public:
    diagmatrix() : sym_band_matrix<T>::sym_band_matrix(){};                                 //default constructor
    diagmatrix(int n) : sym_band_matrix<T>::sym_band_matrix(n, 0) {}                        //构造函数, 构造n级0矩阵
    diagmatrix(diagmatrix &&diagm) noexcept : sym_band_matrix<T>::sym_band_matrix(diagm) {} //moving constructor
    ~diagmatrix() {}                                                                        //Destructor 

    vec<T> linear_f(const vec<T> &b, T zerolim = ZREOLIM) const; //返回以此对角矩阵为系数矩阵，以 b为非齐次项的线性方程组的解
};

//upper triangular (banded) matrix, inheret from sym_band_matrix, which can be treated as its top right port
template <class T>
class uptrimatrix : public sym_band_matrix<T>
{
    template <class Type>
    friend class tridiagmatrix;
    template <class Type>
    friend class lowtrimatrix;
    template <class Type>
    friend void swap(uptrimatrix<Type> &lhs, uptrimatrix<Type> &rhs); //交换两个对称带状矩阵
    template <class Type>
    friend matrix<Type> operator*(const matrix<Type> &ma, const uptrimatrix<Type> &up); //矩阵乘上三角矩阵
    template <class Type>
    friend matrix<Type> operator*(const uptrimatrix<Type> &um, const matrix<Type> &ma); //上三角矩阵乘矩阵
    template <class Type>
    friend uptrimatrix<Type> operator*(const uptrimatrix<Type> &um, const diagmatrix<Type> &dm); //上三角矩阵乘对角矩阵

private:
    bool IsT; //是否是某个下三角矩阵的转至，如果是，则不能释放其存储矩阵
public:
    uptrimatrix() : IsT(0) {}                                                            //default constructor
    uptrimatrix(int n) : sym_band_matrix<T>::sym_band_matrix(n, n - 1), IsT(0) {}        //构造函数, 构造 n 级半带宽 (n-1) 的 0 矩阵
    uptrimatrix(int n, int hbw) : sym_band_matrix<T>::sym_band_matrix(n, hbw), IsT(0) {} //构造函数, 构造 n 级半带宽 hbw 的 0 矩阵
    //uptrimatrix(const lowtrimatrix<T> &lm);                                                 //构造函数，由下三角矩阵转置而来，且始终与之关联
    uptrimatrix(const matrix<T> &ma);                                                                //构造函数, 取矩阵 ma 的右上一半
    uptrimatrix(const uptrimatrix &um) : sym_band_matrix<T>::sym_band_matrix(um), IsT(0) {}          //copy constructor
    uptrimatrix(uptrimatrix &&um) noexcept : sym_band_matrix<T>::sym_band_matrix(um), IsT(um.IsT) {} //moving constructor
    ~uptrimatrix()                                                                                   //Destructor 
    {
        if (IsT)
        {
            sym_band_matrix<T>::halfbandwidth = -1;
            sym_band_matrix<T>::sbm = NULL;
        }
    }
    uptrimatrix &operator=(uptrimatrix um); //overload operator "=" , 保留转置关联

    T getelem(int i, int j) const;         //返回第i行第j列的值
    void changeelem(int i, int j, T elem); //将第i行第i列的值改为elem
    void print(int w = 5) const;           //打印此上三角矩阵
    lowtrimatrix<T> transpose() const;     //返回一个它的转置

    vec<T> linear_f(const vec<T> &b, T zerolim = ZREOLIM) const; //返回以此上三角矩阵为系数矩阵，以 b为非齐次项的线性方程组的解
};

//lower triangular (banded) matrix, inheret from sym_band_matrix, which can be treated as its bottom left port
template <class T>
class lowtrimatrix : public sym_band_matrix<T>
{
    template <class Type>
    friend class tridiagmatrix;
    template <class Type>
    friend class uptrimatrix;
    template <class Type>
    friend class sym_band_matrix;
    template <class Type>
    friend void swap(lowtrimatrix<Type> &lhs, lowtrimatrix<Type> &rhs); //交换两个对称带状矩阵
    template <class Type>
    friend lowtrimatrix<Type> operator*(lowtrimatrix<Type> lm, const diagmatrix<Type> &dm); //下三角矩阵乘对角矩阵

private:
    bool IsT; //是否是某个上三角矩阵的转至，如果是，则不能释放其存储矩阵
public:
    lowtrimatrix() : IsT(0) {}                                                            //default constructor
    lowtrimatrix(int n) : sym_band_matrix<T>::sym_band_matrix(n, n - 1), IsT(0) {}        //构造函数, 构造 n 级半带宽 (n-1) 的 0 矩阵
    lowtrimatrix(int n, int hbw) : sym_band_matrix<T>::sym_band_matrix(n, hbw), IsT(0) {} //构造函数, 构造 n 级半带宽 hbw 的 0 矩阵
    //lowtrimatrix(const uptrimatrix<T> &um);                                                   //构造函数，由一个上三角矩阵转置而来, 且始终与之关联
    lowtrimatrix(const matrix<T> &ma);                                                                 //构造函数, 取矩阵 ma 的左下一半
    lowtrimatrix(const lowtrimatrix &lm) : sym_band_matrix<T>::sym_band_matrix(lm), IsT(0) {}          //copy constructor
    lowtrimatrix(lowtrimatrix &&lm) noexcept : sym_band_matrix<T>::sym_band_matrix(lm), IsT(lm.IsT) {} //copy constructor
    ~lowtrimatrix()                                                                                    //Destructor 
    {
        if (IsT)
        {
            sym_band_matrix<T>::halfbandwidth = -1;
            sym_band_matrix<T>::sbm = NULL;
        }
    }
    lowtrimatrix &operator=(lowtrimatrix lm); //overload operator "=" , 保留转置关联

    T getelem(int i, int j) const;         //返回第i行第j列的值
    uptrimatrix<T> transpose() const;      //返回一个它的转置
    void changeelem(int i, int j, T elem); //将第i行第i列的值改为elem
    void print(int w = 6) const;           //打印此下三角矩阵

    vec<T> linear_f(const vec<T> &b, T zerolim = ZREOLIM) const; //返回以此下三角矩阵为系数矩阵，以b为非齐次项的线性方程组的解
};

//************** definition of functions in class 'vec' ******************
template <class T>
vec<T>::vec(int n, bool c) : dim(n), IsCol(c) //Constructor, construct an n-dimensional zero vector
{
    if (n <= 0)
    {
        dim = 0;
        v = NULL;
    }
    else
    {
        v = new T[n];
        for (int i = 0; i < n; i++)
            v[i] = 0;
    }
}

template <class T>
vec<T>::vec(const vec &ve) : dim(ve.dim), IsCol(ve.IsCol) //copy constructor
{
    if (dim == 0) //Initialize the pointer
        v = NULL;
    else
        v = new T[dim];
    for (int i = 0; i < dim; i++) //copy the elements
        v[i] = ve.v[i];
    //std::cout << "Calling a copy constructor of a vec." << std::endl;
}

template <class T>
vec<T>::vec(vec &&ve) noexcept : dim(ve.dim), v(ve.v), IsCol(ve.IsCol) //moving constructor
{
    ve.v = nullptr;
    ve.dim = 0;
    //std::cout << "Calling a moving constructor of a vec." << std::endl;
}

template <class T>
vec<T>::vec(const matrix<T> &ma, int i, bool is_col) //Constructor, construct from i-th row or i-th column of a matrix
{
    dim = ma.dim;
    if (dim == 0)
        v = NULL;
    else
        v = new T[dim];
    for (int k = 0; k < dim; k++)
        if (is_col)
            v[k] = ma.m[k][i - 1];
        else
            v[k] = ma.m[i - 1][k];
}

template <class T>
vec<T> &vec<T>::operator=(vec<T> rhs) //overload operator = 
{
    swap(*this, rhs);
    return *this;
}

template <class T>
vec<T> &vec<T>::operator+=(const vec<T> &rhs) //overload operator +=
{
    if (rhs.dim != dim) //if the dimension is difference, throw an exception
        throw std::invalid_argument("The dimensions of the two vectors are different!");
    else
        for (int i = 0; i < dim; i++)
            v[i] += rhs.v[i];
    return *this;
}

template <class T>
vec<T> &vec<T>::operator-=(const vec<T> &rhs) //overload operator -=
{
    if (rhs.dim != dim) //if the dimension is difference, throw an exception
        throw std::invalid_argument("The dimensions of the two vectors are different!");
    else //if the dimensions are the same
        for (int i = 0; i < dim; i++)
            v[i] -= rhs.v[i];
    return *this;
}

template <class T>
vec<T> &vec<T>::operator*=(const T &num) //overload operator number mutiplication
{
    for (int i = 0; i < dim; i++)
        v[i] *= num;
    return *this;
}

template <class T>
T vec<T>::getelem(int i) const //return i-th element
{
    if (i <= 0 || i > dim)
        throw std::out_of_range("The element wanted getting is out of bound!");
    return v[i - 1];
}

template <class T>
void vec<T>::changeelem(int i, T elem) //change the i-th element into elem
{
    if (i <= 0 || i > dim)
        throw std::out_of_range("The element you want to change is out of bound!");
    v[i - 1] = elem;
}

template <class T>
void vec<T>::print(int w, int precision) const //print this vector
{
    if (IsCol) // if it is a column vector
        for (int i = 0; i < dim; i++)
            if (precision <= 0)
                std::cout << std::setw(w) << v[i] << std::endl;
            else
            {
                std::cout << std::setprecision(precision) << std::setw(w) << v[i] << std::endl;
                std::cout << std::setprecision(6);
            }
    else
    {
        for (int i = 0; i < dim; i++)
            if (precision <= 0)
                std::cout << std::setw(w) << v[i];
            else
            {
                std::cout << std::setprecision(precision) << std::setw(w) << v[i] << std::endl;
                std::cout << std::setprecision(6);
            }

        std::cout << std::endl;
    }
}

template <class T>
bool vec<T>::exch_elem(int i, int j) //swap two elements
{
    if (i > dim || i <= 0 || j > dim || j <= 0)
        throw std::out_of_range("The elements you want to exchange is out of bound!");
    using std::swap;
    swap(v[i - 1], v[j - 1]);
    // T temp = v[i - 1];
    // v[i - 1] = v[j - 1];
    // v[j - 1] = temp;
    return true;
}

//*************** definition of functions in class 'sym_band_matrix' *****************
template <class T>
T **sym_band_matrix<T>::copysbm() const //复制一个和sbm相同的二维数组，并返回其头指针
{
    if (halfbandwidth == -1)
        return NULL;
    T **temp = new T *[halfbandwidth + 1];
    for (int i = 0; i <= halfbandwidth; i++)
    {
        temp[i] = new T[dim - i];
        for (int j = 0; j < dim - i; j++)
            temp[i][j] = sbm[i][j];
    }
    return temp;
}

template <class T>
void sym_band_matrix<T>::dltsbm() //把数组 sbm 释放掉
{
    for (int i = 0; i <= halfbandwidth; i++)
        delete[] sbm[i];
    delete[] sbm;
}

template <class T>
sym_band_matrix<T>::sym_band_matrix(int n, int hbw) : dim(n), halfbandwidth(hbw) //构造函数, 构造 n 级,半带宽 hbw 的 0 矩阵
{
    if (n <= 0)
    {
        dim = 0;
        halfbandwidth = -1;
        sbm = NULL;
    }
    else if (hbw < 0)
    {
        halfbandwidth = -1;
        sbm = NULL;
    }
    else
    {
        if (hbw > n - 1) //如果半带宽太大了，就至多设为 (dim-1), 把矩阵填满
            halfbandwidth = n - 1;
        sbm = new T *[halfbandwidth + 1];
        for (int i = 0; i <= halfbandwidth; i++)
        {
            sbm[i] = new T[dim - i];
            for (int j = 0; j < dim - i; j++)
                sbm[i][j] = 0;
        }
    }
}

template <class T>
sym_band_matrix<T>::sym_band_matrix(const sym_band_matrix &sbma) //copy constructor
{
    dim = sbma.dim;
    halfbandwidth = sbma.halfbandwidth;
    sbm = sbma.copysbm();
}

template <class T>
sym_band_matrix<T>::sym_band_matrix(sym_band_matrix &&sbma) noexcept : dim(sbma.dim), halfbandwidth(sbma.halfbandwidth), sbm(sbma.sbm) //moving constructor
{
    std::cout << "Moving Constructor" << std::endl;
    sbma.halfbandwidth = -1;
    sbma.sbm = nullptr;
}

template <class T>
sym_band_matrix<T>::~sym_band_matrix() //Destructor 
{
    dltsbm();
    //std::cout << "Calling the destructor of sym_band_matrix!" << std::endl;
}

template <class T>
sym_band_matrix<T> &sym_band_matrix<T>::operator=(sym_band_matrix<T> sbma) //overload operator "=" 
{
    swap(*this, sbma);
    return *this;
}

template <class T>
T sym_band_matrix<T>::getelem(int i, int j) const //return the element in row i and column j
{
    if (i > dim || i <= 0 || j > dim || j <= 0)
        throw std::out_of_range("The element you want to get is out of bound!");
    if (i - j > halfbandwidth || j - i > halfbandwidth)
        return 0;
    if (i >= j)
        return sbm[i - j][j - 1];
    return sbm[j - i][i - 1];
}

template <class T>
void sym_band_matrix<T>::changeelem(int i, int j, T elem) //change the element in row i, column j into elem(同时也会改变对称元素)
{
    if (i > dim || i <= 0 || j > dim || j <= 0)
        throw std::out_of_range("The element you want to change is out of bound!");
    if (i - j > halfbandwidth || j - i > halfbandwidth)
        throw std::out_of_range("The element cannot be changed because it is out of the band!");
    if (i >= j)
        sbm[i - j][j - 1] = elem;
    else
        sbm[j - i][i - 1] = elem;
}

template <class T>
void sym_band_matrix<T>::print(int w) const //print this matrix
{
    for (int i = 1; i <= dim; i++) // i表示矩阵真实行号，而不是存储指标
    {
        for (int j = 1; j <= dim; j++) // j表示矩阵真实列号，而不是存储指标
            std::cout << std::setw(w) << getelem(i, j);
        std::cout << std::endl;
    }
}

template <class T>
bool sym_band_matrix<T>::Cholesky(lowtrimatrix<T> &lm, diagmatrix<T> &dm, T zerolim) // Cholesky 分解
{
    if (dim == 0 || halfbandwidth < 0) // 如果是零维的, 或者没有宽度, 就不玩了
        return false;
    //构造存储数据的数组
    T *d = new T[dim];
    T **t = new T *[halfbandwidth + 1];
    T **l = new T *[halfbandwidth + 1];
    for (int i = 1; i <= halfbandwidth; i++) //t 的主对角元没用，为了指标方便仍然保留
        t[i] = new T[dim - i];
    for (int i = 0; i <= halfbandwidth; i++)
        l[i] = new T[dim - i];
    // 开始 Cholesky 分解
    for (int j = 0; j < dim; j++) //先把 l 主对角元变为 1
        l[0][j] = 1;
    d[0] = sbm[0][0];              //d_1=a_1_1
    for (int i = 2; i <= dim; i++) //i 是真实矩阵行号
    {
        if (std::abs(d[i - 2]) < std::abs(zerolim)) //d_(i-1), 即刚刚构造的 d 十分接近 0, 则退出
        {
            delete[] d;
            for (int i = 1; i <= halfbandwidth; i++)
                delete[] t[i];
            delete[] t;
            for (int i = 0; i <= halfbandwidth; i++)
                delete[] l[i];
            delete[] l;
            std::cout << "对称矩阵不可逆！" << std::endl;
            return false;
        }
        int beginindex = std::max(1, i - halfbandwidth); //第 i 行列的起始下标
        T tempsum = 0;
        for (int j = beginindex; j <= i - 1; j++) //j 是矩阵真实列号
        {
            for (int k = 1; k <= j - 1; k++)
                tempsum += t[i - k][k - 1] * l[j - k][k - 1];
            t[i - j][j - 1] = sbm[i - j][j - 1] - tempsum; //t_i_j=a_i_j-[t_i_1*l_j_1+t_i_2*l_j_2+...+t_i_(j-1)*l_j_(j-1)]
            l[i - j][j - 1] = t[i - j][j - 1] / d[j - 1];  //t_i_j = l_i_j * d_j
            tempsum = 0;
        }
        for (int j = beginindex; j <= i - 1; j++)
            tempsum += t[i - j][j - 1] * l[i - j][j - 1];
        d[i - 1] = sbm[0][i - 1] - tempsum; //d_i=a_i_i - [t_i_1*l_i_1+t_i_2*l_j_2+...+t_i_(i-1)*l_i_(i-1)]
    }
    if (std::abs(d[dim - 1]) < std::abs(zerolim)) //d_dim 十分接近 0, 则退出
    {
        delete[] d;
        for (int i = 1; i <= halfbandwidth; i++)
            delete[] t[i];
        delete[] t;
        for (int i = 0; i <= halfbandwidth; i++)
            delete[] l[i];
        delete[] l;
        std::cerr << "对称矩阵不可逆！" << std::endl;
        return false;
    }
    //生成需要的矩阵
    dm.dltsbm();
    dm.dim = lm.dim = dim;
    dm.sbm = new T *[1];
    dm.sbm[0] = d;
    dm.halfbandwidth = 0;
    lm.dltsbm();
    lm.halfbandwidth = halfbandwidth;
    lm.IsT = 0;
    lm.sbm = l;
    //将辅助数组 t 释放掉
    for (int i = 1; i <= halfbandwidth; i++)
        delete[] t[i];
    delete[] t;

    return true;
}

template <class T>
vec<T> sym_band_matrix<T>::linear_f(const vec<T> &b, T zerolim) //返回以此上三角矩阵为系数矩阵，以 b为非齐次项的线性方程组的解
{
    vec<T> so;
    if (sym_band_matrix<T>::dim != b.dim)
        throw std::invalid_argument("The dimensions of the vector and the matrix are different!");
    diagmatrix<T> dm;
    lowtrimatrix<T> lm;
    bool Chol = Cholesky(lm, dm, std::abs(zerolim));
    // if(dim>10)
    // {
    // std::cout << "------------ dm ---------- n=" << dim << " ---------" << std::endl;
    // dm.print(8);
    // std::cout << "------------ lm ---------- n=" << dim << " ---------" << std::endl;
    // lm.print(8);
    // }
    diagmatrix<T> dm2(dm.dim);
    for (int i = 1; i <= dm.dim; i++)
        dm2.changeelem(i, i, std::sqrt(dm.getelem(i, i)));
    std::cout << "---------- H -----------" << std::endl;
    (lm * dm2).print(11);
    if (!Chol)
        return so;
    vec<T> y(lm.linear_f(b));            // Ly=b
    vec<T> z(dm.linear_f(y));            //Dz=y
    return (lm.transpose()).linear_f(z); //L^T x=z
}

template <class T>
matrix<T> sym_band_matrix<T>::Jacobi(int i, int j, T zerolim) //一次 Jacobi 迭代, (*this)=QᵗAQ, 返回 Q, (*this)变为 A
{
    if (i > dim || i <= 0 || j > dim || j <= 0 || i - j > halfbandwidth || j - i > halfbandwidth || i == j)
        throw std::runtime_error("The element you want to rotate is illegal!");
    matrix<T> Q(dim); //待返回的正交矩阵
    if (i < j)        //使得 i>=j
        std::swap(i, j);
    for (int k = 1; k <= dim; k++)
        Q.changeelem(k, k, 1);
    if (std::abs(sbm[i - j][j - 1]) < zerolim)
        return Q;
    T s, c, t;
    T eta = (sbm[0][i - 1] - sbm[0][j - 1]) / (2 * sbm[i - j][j - 1]); // η=(a_ii-a_jj)/(2a_ij)
    if (eta >= 0)
        t = 1 / (eta + sqrt(1 + eta * eta));
    else
        t = 1 / (eta - sqrt(1 + eta * eta));
    c = 1 / sqrt(1 + t * t);
    s = c * t;
    Q.changeelem(i, i, c);
    Q.changeelem(j, j, c);
    Q.changeelem(i, j, s);
    Q.changeelem(j, i, -s);
    // for (int k = 1; k <= dim; k++) //先对两列做变换, 即先计算中间矩阵乘右边矩阵; 假设 k 比 i,j 小
    // {
    //     if (k < i && i - k <= halfbandwidth) //b_ki=sa_kj+ca_ki
    //     {
    //         sbm[i - k][k - 1] *= c;
    //         if (k < j)
    //             sbm[i - k][k - 1] += s * sbm[j - k][k - 1];
    //     }
    //     if (k < j && j - k <= halfbandwidth) //b_kj=ca_kj-sa_ki
    //     {
    //         sbm[j - k][k - 1] /= c;
    //         sbm[j - k][k - 1] -= t * sbm[i - k][k - 1];
    //     }
    // }
    // sbm[i - j][j - 1] = 0;
    // sbm[0][j - 1] *= c * c;
    // sbm[0][j - 1] -= s * s * sbm[0][i - 1];
    // sbm[0][i - 1] *= 1 - s * s / c * c;
    // sbm[0][i - 1] -= t * t * sbm[0][j - 1];
    T *ai = nullptr;
    ai = new T[dim]; //ai[k-1]=a_ki
    for (int k = std::max(1, i - halfbandwidth); k < i; k++)
        ai[k - 1] = sbm[i - k][k - 1];
    for (int k = i; k <= i + halfbandwidth; k++)
        ai[k - 1] = sbm[k - i][i - 1];
    for (int k = std::max(1, i - halfbandwidth); k <= std::min(dim, i + halfbandwidth); k++) //修改第 i 列(不含第 j 行)
    {
        if (k < i && k != j)
            sbm[i - k][k - 1] *= c;
        else if (k > i)
            sbm[k - i][i - 1] *= c;
        if (k < j)
            sbm[i - k][k - 1] += s * sbm[j - k][k - 1];
        else if (k > j && k < i)
            sbm[i - k][k - 1] += s * sbm[k - j][j - 1];
        else if (k > i && k <= j + halfbandwidth)
            sbm[k - i][i - 1] += s * sbm[k - j][j - 1];
    }
    for (int k = std::max(1, j - halfbandwidth); k <= std::min(dim, j + halfbandwidth); k++) //修改第 j 列
    {
        if (k < j)
            sbm[j - k][k - 1] *= c;
        else if (k > j && k != i)
            sbm[k - j][j - 1] *= c;
        if (k >= i - halfbandwidth && k < j)
            sbm[j - k][k - 1] -= s * ai[k - 1];
        else if (k > j && k != i)
            sbm[k - j][j - 1] -= s * ai[k - 1];
    }
    sbm[0][i - 1] *= c * c;
    sbm[0][i - 1] += 2 * s * c * sbm[i - j][j - 1] + s * s * sbm[0][j - 1];
    sbm[0][j - 1] *= c * c;
    sbm[0][j - 1] += s * s * ai[i - 1] - 2 * s * c * sbm[i - j][j - 1];
    sbm[i - j][j - 1] = 0;

    delete[] ai;
    return Q;
}

template <class T>
T sym_band_matrix<T>::Sturm(int n, int stepnum, T zerolim) // Sturm 方法求从小到大第 n 个本征值
{
    if (n <= 0 || n > dim)
        throw std::runtime_error("The eigenstate does not exist!");
    if (halfbandwidth != 1)
        throw std::runtime_error("Sturm method cannot be used!");
    T lowlim = sbm[0][0] - std::abs(sbm[1][0]), uplim = sbm[0][0] + std::abs(sbm[1][0]); //特征值下限和上限
    for (int i = 1; i < dim - 1; i++)
    {
        if (lowlim > sbm[0][i] - std::abs(sbm[1][i]) - std::abs(sbm[1][i - 1]))
            lowlim = sbm[0][i] - std::abs(sbm[1][i]) - std::abs(sbm[1][i - 1]);
        if (uplim < sbm[0][i] + std::abs(sbm[1][i]) + std::abs(sbm[1][i - 1]))
            uplim = sbm[0][i] + std::abs(sbm[1][i]) + std::abs(sbm[1][i - 1]);
    }
    if (lowlim > sbm[0][dim] - std::abs(sbm[1][dim - 1]))
        lowlim = sbm[0][dim] - std::abs(sbm[1][dim - 1]);
    if (uplim < sbm[0][dim] + std::abs(sbm[1][dim - 1]))
        uplim = sbm[0][dim] + std::abs(sbm[1][dim - 1]);
    T mid = (lowlim + uplim) / 2;
    for (int i = 0; i < stepnum; i++)
    {
        if (uplim - lowlim < std::abs(zerolim))
        {
            //std::cout << "i=" << i << std::endl;
            return mid;
        }
        if ((1 + i) % 5 == 0)
            std::cout << mid << std::endl;
        int s = 0;                    //s 表示 {p_1(mid), p_2(mid), ... , p_dim(mid)}的变号次数, 我们把 0 视为正
        T q = sbm[0][0] - mid;        // q_i=p_i/p_(i-1), q_1=p_1(mid)
        for (int j = 1; j < dim; j++) //计算变号数 s
        {
            if (q < 0)
                s++;
            // if (std::abs(q) < std::abs(zerolim)) //如果 q 很接近 0
            q = sbm[0][j] - mid - sbm[1][j - 1] * sbm[1][j - 1] / q;
        }
        if (q < 0)
            s++;
        if (s >= n)
            uplim = mid;
        else
            lowlim = mid;
        mid = (lowlim + uplim) / 2;
    }
    return mid;
}

//***************** definition of functions in class 'sparse_matrix' ***************
template <class T>
sparse_matrix<T>::sparse_matrix(int n) : dim(n), spm(NULL) //构造函数, 构造 n 行 0 列矩阵
{
    if (n <= 0)
    {
        std::cerr << "The " << n << "-dimision matrix fails to be founded!" << std::endl;
        dim = 0;
    }
    else
        spm = new std::list<spmnode<T>>[n];
}

template <class T>
sparse_matrix<T>::sparse_matrix(const sparse_matrix &sma) : dim(sma.dim) //copy constructor
{
    delete[] spm;
    if (dim == 0)
        spm = NULL;
    else
    {
        spm = new std::list<spmnode<T>>[dim];
        for (int i = 0; i < dim; i++)
            spm[i] = sma.spm[i];
    }
}

template <class T>
sparse_matrix<T>::sparse_matrix(sparse_matrix &&sma) noexcept : dim(sma.dim), spm(sma.spm) //moving constructor
{
    sma.dim = 0;
    sma.spm = nullptr;
}

template <class T>
bool sparse_matrix<T>::changeelem(int i, int j, T e) //change the element in row i, column i into elem
{
    if (i > dim || i <= 0 || j > dim || j <= 0)
        throw std::out_of_range("The element you want to change is out of bound!");
    typename std::list<spmnode<T>>::iterator it = spm[i - 1].end();
    if (spm[i - 1].empty() || (--it)->index < j) //如果最后一个元素列指标都小于 j 或者列表为空, 那直接加在后面就可以了
    {
        spmnode<T> node{e, j};
        spm[i - 1].push_back(node);
        return true;
    }
    for (it = spm[i - 1].begin(); it != spm[i - 1].end(); it++)
    {
        if (it->index == j) //看看 it 所指元素是不是在第 j 列
        {
            it->elem = e;
            return true;
        }
        if (it->index > j) //看看第 j 列在不在 it 和下一个结点之间
        {
            spmnode<T> node{e, j};
            spm[i - 1].insert(it, node);
            return true;
        }
    }
    return false;
}

template <class T>
void sparse_matrix<T>::print(int w) const
{
    for (int i = 1; i <= dim; i++) //打印第 i 行
    {
        int last_index = 0;
        for (typename std::list<spmnode<T>>::iterator it = spm[i - 1].begin(); it != spm[i - 1].end(); it++)
        {
            for (int k = 1; k < it->index - last_index; k++)
                std::cout << std::setw(w) << 0;
            last_index = it->index;
            std::cout << std::setw(w) << it->elem;
        }
        for (int k = 1; k <= dim - last_index; k++) //输入剩下的 0
            std::cout << std::setw(w) << 0;
        std::cout << std::endl;
    }
}

template <class T>
sparse_matrix<T> &sparse_matrix<T>::operator=(sparse_matrix<T> sma) //overload operator = 
{
    swap(*this, sma);
    return *this;
}

//***************** definition of functions in class 'tridiagmatrix' ***************
template <class T>
tridiagmatrix<T>::tridiagmatrix(int n) : dim(n) //constructor, construct an n-dimensional zero matrix
{
    if (n <= 0)
    {
        dim = 0;
        tdm = NULL;
    }
    else
    {
        tdm = new T *[3];
        if (dim >= 2) //存储对角线左下的次对角线
            tdm[0] = new T[dim - 1];
        else
            tdm[0] = NULL;
        for (int j = 0; j < dim - 1; j++)
            tdm[0][j] = 0;
        tdm[1] = new T[dim]; //存储对角线
        for (int j = 0; j < dim; j++)
            tdm[1][j] = 0;
        tdm[2] = new T[dim]; //存储对角线右上的次对角线, 这里多留了一个位置，是为了方便读取列号
        for (int j = 0; j < dim; j++)
            tdm[2][j] = 0;
    }
}

template <class T>
tridiagmatrix<T>::tridiagmatrix(const tridiagmatrix &tdma) //copy constructor
{
    if (tdma.dim <= 0)
    {
        dim = 0;
        tdm = NULL;
    }
    else
    {
        dim = tdma.dim;
        tdm = new T *[3];
        if (dim >= 2) //存储对角线左下的次对角线
            tdm[0] = new T[dim - 1];
        else
            tdm[0] = NULL;
        for (int j = 0; j < dim - 1; j++)
            tdm[0][j] = tdma.tdm[0][j];
        tdm[1] = new T[dim]; //存储对角线
        for (int j = 0; j < dim; j++)
            tdm[1][j] = tdma.tdm[1][j];
        tdm[2] = new T[dim]; //存储对角线右上的次对角线, 这里多留了一个位置，是为了方便读取列号
        for (int j = 0; j < dim; j++)
            tdm[2][j] = tdma.tdm[2][j];
    }
}

template <class T>
tridiagmatrix<T>::tridiagmatrix(tridiagmatrix &&tdma) noexcept : dim(tdma.dim), tdm(tdma.tdm) //moving constructor
{
    tdma.dim = 0;
    tdma.tdm = nullptr;
}

template <class T>
tridiagmatrix<T>::~tridiagmatrix() //Destructor 
{
    if (dim >= 1)
    {
        for (int i = 0; i < 3; i++)
            delete[] tdm[i];
        delete[] tdm;
    }
}

template <class T>
tridiagmatrix<T> &tridiagmatrix<T>::operator=(tridiagmatrix tdma) //overload operator "=" 
{
    swap(*this, tdma);
    return *this;
}

template <class T>
T tridiagmatrix<T>::getelem(int i, int j) const //return the element in row i and column j
{
    if (i > dim || i <= 0 || j > dim || j <= 0)
        throw std::out_of_range("The element you want to get is out of bound!");
    else if (i - j > 1 || j - i > 1)
        return 0;
    else
        return tdm[j - i + 1][j - 1];
}

template <class T>
void tridiagmatrix<T>::changeelem(int i, int j, T elem) //change the element in row i, column j into elem
{
    if (i > dim || i <= 0 || j > dim || j <= 0)
        throw std::out_of_range("The element you want to change is out of bound!");
    if (i - j > 1 || j - i > 1)
        std::cerr << "The element (" << i << ',' << j << ") cannot be changed!" << std::endl;
    else
        tdm[j - i + 1][j - 1] = elem;
}

template <class T>
void tridiagmatrix<T>::print(int w) const //print this tridiagonal matrix
{
    for (int i = 1; i <= dim; i++) //i表示矩阵真实行号，而不是存储指标
    {
        for (int j = 1; j <= dim; j++) //j表示矩阵真实列号，而不是存储指标
            std::cout << std::setw(w) << getelem(i, j);
        std::cout << std::endl;
    }
}

template <class T>
bool tridiagmatrix<T>::LU(lowtrimatrix<T> &lm, uptrimatrix<T> &um, T zerolim) const //三对角矩阵 LU 分解
{
    uptrimatrix<T> um1(dim, 1);  //LU 分解后的上三角矩阵
    lowtrimatrix<T> lm1(dim, 1); //LU 分解后的下三角矩阵
    if (dim == 1)
    {
        um1.sbm[0][0] = 1.0;
        lm1.sbm[0][0] = tdm[1][0];
    }
    else if (dim > 1)
    {
        for (int j = 2; j <= dim; j++) //j 是矩阵真实列号
            lm1.sbm[1][j - 2] = tdm[0][j - 2];
        lm1.sbm[0][0] = tdm[1][0];
        for (int i = 1; i <= dim - 1; i++) //i 是矩阵真实行号
        {
            if (std::abs(lm1.sbm[0][i - 1]) < std::abs(zerolim))
            {
                std::cout << "三对角矩阵不可逆!" << std::endl;
                return false;
            }
            else
            {
                um1.sbm[1][i - 1] = tdm[2][i] / lm1.sbm[0][i - 1];             //u_i=c_i / l_i
                lm1.sbm[0][i] = tdm[1][i] - tdm[0][i - 1] * um1.sbm[1][i - 1]; //l_(i+1)=b_(i+1)-a_(i+1)u_i
            }
        }
        for (int i = 1; i <= dim; i++)
            um1.sbm[0][i - 1] = 1;
    }
    swap(lm, lm1);
    swap(um, um1);
    return true;
}

template <class T>
vec<T> tridiagmatrix<T>::linear_f(const vec<T> &b, T zerolim) const //返回以此三对角矩阵为系数矩阵，以b为非齐次项的线性方程组的解
{
    vec<T> so;
    if (dim != b.dim)
        throw std::invalid_argument("The dimensions of the vector and the matrix are different!");
    uptrimatrix<T> um1;    //LU 分解后的上三角矩阵
    lowtrimatrix<T> lm1;   //LU 分解后的下三角矩阵
    LU(lm1, um1, zerolim); //LU 分解
    vec<T> y(lm1.linear_f(b));
    if (y.dim == 0) //说明L不可逆, 返回了一个零维向量
        return so;
    vec<T> x(um1.linear_f(y));
    return x;
    // if (dim == 1)
    // {
    //     if (std::abs(tdm[1][0]) < std::abs(zerolim))
    //     {
    //         std::cout << "系数矩阵不可逆!" << std::endl;
    //         return so;
    //     }
    //     vec<T> z(1);
    //     z.changeelem(1, b.v[0] / tdm[1][0]);
    //     return z;
    // }

    // uptrimatrix<T> um1(dim, 1);  //LU 分解后的上三角矩阵
    // lowtrimatrix<T> lm1(dim, 1); //LU 分解后的下三角矩阵

    // for (int j = 2; j <= dim; j++) //j 是矩阵真实列号
    //     lm1.sbm[1][j - 2] = tdm[0][j - 2];
    // lm1.sbm[0][0] = tdm[1][0];
    // for (int i = 1; i <= dim - 1; i++) //i 是矩阵真实行号
    // {
    //     if (std::abs(lm1.sbm[0][i - 1]) < std::abs(zerolim))
    //     {
    //         std::cout << "系数矩阵不可逆!" << std::endl;
    //         return so;
    //     }
    //     else
    //     {
    //         um1.sbm[1][i - 1] = tdm[2][i] / lm1.sbm[0][i - 1];             //u_i=c_i / l_i
    //         lm1.sbm[0][i] = tdm[1][i] - tdm[0][i - 1] * um1.sbm[1][i - 1]; //l_(i+1)=b_(i+1)-a_(i+1)u_i
    //     }
    // }
    // for (int i = 1; i <= dim; i++)
    //     um1.sbm[0][i - 1] = 1;

    // // std::cout << "---------------lm1----------------" << std::endl;
    // // lm1.print(15);
    // // std::cout << "---------------um1----------------" << std::endl;
    // // um1.print(15);

    // vec<T> y(lm1.linear_f(b));
    // if (y.dim == 0) //说明L不可逆, 返回了一个零维向量
    //     return so;
    // // cout << "---------------y----------------" << endl;
    // // y.print();
    // vec<T> x(um1.linear_f(y));
    // return x;
}

//***************** definition of functions in class 'diagmatrix' ***************
template <class T>
vec<T> diagmatrix<T>::linear_f(const vec<T> &b, T zerolim) const //返回以此对角矩阵为系数矩阵，以 b为非齐次项的线性方程组的解
{
    vec<T> so;
    if (sym_band_matrix<T>::dim != b.dim)
    {
        std::cerr << "系数矩阵和非齐次项的维数不同！" << std::endl;
        return so;
    }
    vec<T> x(b.dim);
    for (int i = 1; i <= sym_band_matrix<T>::dim; i++)
    {
        if (std::abs(sym_band_matrix<T>::sbm[0][i - 1]) < std::abs(zerolim))
        {
            std::cout << "系数矩阵不可逆！" << std::endl;
            return so;
        }
        x.v[i - 1] = b.v[i - 1] / sym_band_matrix<T>::sbm[0][i - 1];
    }
    return x;
}

//*************** definition of functions in class 'uptrimatrix' *****************
template <class T>
uptrimatrix<T>::uptrimatrix(const matrix<T> &ma) : IsT(0) //构造函数, 取矩阵 ma 的右上一半
{
    sym_band_matrix<T>::dim = ma.dim;
    sym_band_matrix<T>::halfbandwidth = ma.dim - 1;
    if (ma.dim >= 1)
        sym_band_matrix<T>::sbm = new T *[ma.dim];
    else
        sym_band_matrix<T>::sbm = NULL;
    for (int i = 0; i < ma.dim; i++)
    {
        sym_band_matrix<T>::sbm[i] = new T[ma.dim - i];
        for (int j = 0; j < ma.dim - i; j++)
            sym_band_matrix<T>::sbm[i][j] = ma.m[j][i + j];
    }
}

template <class T>
uptrimatrix<T> &uptrimatrix<T>::operator=(uptrimatrix<T> um) //overload operator "=" 
{
    swap(*this, um);
    return *this;
}

template <class T>
T uptrimatrix<T>::getelem(int i, int j) const //return the element in row i and column j
{
    if (i > sym_band_matrix<T>::dim || i <= 0 || j > sym_band_matrix<T>::dim || j <= 0)
    {
        std::cerr << "The element (" << i << ',' << j << ") is not exist!" << std::endl;
        return 0;
    }
    if (i > j || j - i > sym_band_matrix<T>::halfbandwidth)
        return 0;
    return sym_band_matrix<T>::sbm[j - i][i - 1];
}

template <class T>
void uptrimatrix<T>::changeelem(int i, int j, T elem) //change the element in row i, column i into elem
{
    if (i > sym_band_matrix<T>::dim || i <= 0 || j > sym_band_matrix<T>::dim || j <= 0)
    {
        std::cerr << "The element (" << i << ',' << j << ") is not exist!" << std::endl;
        return;
    }
    if (i > j || j - i > sym_band_matrix<T>::halfbandwidth)
    {
        std::cerr << "The element (" << i << ',' << j << ") cannot be changed!" << std::endl;
        return;
    }
    sym_band_matrix<T>::sbm[j - i][i - 1] = elem;
}

template <class T>
lowtrimatrix<T> uptrimatrix<T>::transpose() const //return the transpose
{
    lowtrimatrix<T> l;
    l.IsT = true;
    l.dim = sym_band_matrix<T>::dim;
    l.halfbandwidth = sym_band_matrix<T>::halfbandwidth;
    l.sbm = sym_band_matrix<T>::sbm;
    return l;
}

template <class T>
void uptrimatrix<T>::print(int w) const //print this upper triangular matrix
{
    for (int i = 1; i <= sym_band_matrix<T>::dim; i++) //i是矩阵真实行号
    {
        for (int j = 1; j <= sym_band_matrix<T>::dim; j++) //j是矩阵真实列号
            std::cout << std::setw(w) << getelem(i, j);
        std::cout << std::endl;
    }
}

template <class T>
vec<T> uptrimatrix<T>::linear_f(const vec<T> &b, T zerolim) const //返回以此上三角矩阵为系数矩阵，以b为非齐次项的线性方程组的解
{
    vec<T> so;
    if (sym_band_matrix<T>::dim != b.dim)
        throw std::invalid_argument("The dimensions of the vector and the matrix are different!");
    vec<T> x(sym_band_matrix<T>::dim);
    for (int i = sym_band_matrix<T>::dim; i > 0; i--) //i是真实矩阵行号, 不是数组存储指标
    {
        if (std::abs(sym_band_matrix<T>::sbm[0][i - 1]) < std::abs(zerolim)) //如果第行 i 对角元近似为 0
        {
            std::cout << i << std::endl;
            std::cout << "The upper triangle matrix is not invertible!" << std::endl;
            return so;
        }
        T tempsum = 0;
        if (i >= sym_band_matrix<T>::dim - sym_band_matrix<T>::halfbandwidth) //如果i比较大, 还没有到达带宽之外
            for (int j = i + 1; j <= sym_band_matrix<T>::dim; j++)            //j是真实矩阵行号, 不是数组存储指标
                tempsum += sym_band_matrix<T>::sbm[j - i][i - 1] * x.v[j - 1];
        else
            for (int j = i + 1; j <= i + sym_band_matrix<T>::halfbandwidth; j++) //j是真实矩阵行号, 不是数组存储指标
                tempsum += sym_band_matrix<T>::sbm[j - i][i - 1] * x.v[j - 1];
        x.v[i - 1] = (b.v[i - 1] - tempsum) / sym_band_matrix<T>::sbm[0][i - 1];
    }
    return x;
}

// ************** definition of functions in class 'lowtrimatrix' *******************
template <class T>
lowtrimatrix<T>::lowtrimatrix(const matrix<T> &ma) : IsT(false) //构造函数, 取矩阵 ma 的左下一半
{
    sym_band_matrix<T>::dim = ma.dim;
    sym_band_matrix<T>::halfbandwidth = ma.dim - 1;
    if (ma.dim >= 1)
        sym_band_matrix<T>::sbm = new T *[ma.dim];
    else
        sym_band_matrix<T>::sbm = NULL;
    for (int i = 0; i < ma.dim; i++)
    {
        sym_band_matrix<T>::sbm[i] = new T[ma.dim - i];
        for (int j = 0; j < ma.dim - i; j++)
            sym_band_matrix<T>::sbm[i][j] = ma.m[i + j][j];
    }
}

template <class T>
lowtrimatrix<T> &lowtrimatrix<T>::operator=(lowtrimatrix lm) //overload operator "=" 
{
    swap(*this, lm);
    return *this;
}

template <class T>
T lowtrimatrix<T>::getelem(int i, int j) const //return the element in row i and column j
{
    if (i > sym_band_matrix<T>::dim || i <= 0 || j > sym_band_matrix<T>::dim || j <= 0)
    {
        std::cerr << "The element (" << i << ',' << j << ") is not exist!" << std::endl;
        return 0;
    }
    if (i < j || i - j > sym_band_matrix<T>::halfbandwidth)
        return 0;
    return sym_band_matrix<T>::sbm[i - j][j - 1];
}

template <class T>
uptrimatrix<T> lowtrimatrix<T>::transpose() const //返回一个转置
{
    uptrimatrix<T> u;
    u.IsT = true;
    u.dim = sym_band_matrix<T>::dim;
    u.halfbandwidth = sym_band_matrix<T>::halfbandwidth;
    u.sbm = sym_band_matrix<T>::sbm;
    return u;
}

template <class T>
void lowtrimatrix<T>::changeelem(int i, int j, T elem) //change the element in row i, column i into elem
{
    if (i > sym_band_matrix<T>::dim || i <= 0 || j > sym_band_matrix<T>::dim || j <= 0)
    {
        std::cerr << "The element (" << i << ',' << j << ") is not exist!" << std::endl;
        return;
    }
    if (i < j || i - j > sym_band_matrix<T>::halfbandwidth)
    {
        std::cerr << "The element (" << i << ',' << j << ") cannot be changed!" << std::endl;
        return;
    }
    sym_band_matrix<T>::sbm[i - j][j - 1] = elem;
}

template <class T>
void lowtrimatrix<T>::print(int w) const //print this lower triangular matrix
{
    for (int i = 1; i <= sym_band_matrix<T>::dim; i++) //i 是矩阵真实行号
    {
        for (int j = 1; j <= sym_band_matrix<T>::dim; j++) //j 是矩阵真实列号
            std::cout << std::setw(w) << getelem(i, j);
        std::cout << std::endl;
    }
}

template <class T>
vec<T> lowtrimatrix<T>::linear_f(const vec<T> &b, T zerolim) const //返回以此下三角矩阵为系数矩阵，以 b 为非齐次项的线性方程组的解
{
    vec<T> so;
    if (sym_band_matrix<T>::dim != b.dim)
    {
        std::cerr << "系数矩阵和非齐次项的维数不同！" << std::endl;
        return so;
    }
    vec<T> x(sym_band_matrix<T>::dim);
    for (int i = 1; i <= sym_band_matrix<T>::dim; i++) //i 是矩阵真实行号，不是存储指标
    {
        if (std::abs(sym_band_matrix<T>::sbm[0][i - 1]) < std::abs(zerolim))
        {
            std::cout << "The lower triangle matrix is not invertible!" << std::endl;
            return so;
        }
        T tempsum = 0;
        if (i <= sym_band_matrix<T>::halfbandwidth + 1)
            for (int j = 1; j < i; j++) //j 是矩阵真实列号
                tempsum += sym_band_matrix<T>::sbm[i - j][j - 1] * x.v[j - 1];
        else
            for (int j = i - sym_band_matrix<T>::halfbandwidth; j < i; j++) //j 是矩阵真实列号
                tempsum += sym_band_matrix<T>::sbm[i - j][j - 1] * x.v[j - 1];
        x.v[i - 1] = (b.v[i - 1] - tempsum) / sym_band_matrix<T>::sbm[0][i - 1];
    }
    return x;
}

//************* definition of functions in class 'matrix' *******************
template <class T>
matrix<T>::matrix() //default constructor
{
    std::cout << " Default construct function" << std::endl;
    matrix::dim = 0;
    matrix::m = NULL;
}

template <class T>
matrix<T>::matrix(int n) //构造函数, 构造n级0矩阵
{
    //std::cout << "Construct function (int " << n << ")" << std::endl;
    if (n <= 0)
    {
        matrix::dim = 0;
        matrix::m = NULL;
    }
    else
    {
        matrix::dim = n;
        matrix::m = new T *[n]; //申请行向量数组, 共n行
        for (int i = 0; i < n; i++)
        {
            matrix::m[i] = new T[n]; //申请列, 第(i+1)存于m[i], i=0,1,...,n-1.
            for (int j = 0; j < n; j++)
                matrix::m[i][j] = 0;
        }
    }
}

template <class T>
matrix<T>::matrix(const matrix &ma) //copy constructor
{
    dim = ma.dim;
    matrix::m = new T *[dim]; //申请行向量数组, 共dim行
    for (int i = 0; i < dim; i++)
    {
        matrix::m[i] = new T[dim]; //申请列, 第(i+1)存于m[i], i=0,1,...,dim-1.
        for (int j = 0; j < dim; j++)
            matrix::m[i][j] = ma.m[i][j];
    }
    std::cout << "Calling the copy construncor of a matrix!" << std::endl;
}

template <class T>
matrix<T>::matrix(matrix &&ma) noexcept : dim(ma.dim), m(ma.m) //moving constructor
{
    ma.dim = 0;
    ma.m = nullptr;
}

template <class T>
matrix<T>::~matrix() //Destructor 
{
    for (int i = 0; i < matrix::dim; i++)
        delete[] matrix::m[i];
    delete[] matrix::m;
    //std::cout << "Destructor " << std::endl;
}

template <class T>
matrix<T> &matrix<T>::operator=(matrix<T> ma) //overload operator = 
{
    swap(ma, *this);
    return *this;
}

template <class T>
T matrix<T>::getelem(int i, int j) const //return the element in row i and column j
{
    if (i > matrix::dim || i <= 0 || j > matrix::dim || j <= 0)
    {
        std::cerr << "The element (" << i << ',' << j << ") is not exist!" << std::endl;
        return 0;
    }
    else
        return matrix::m[i - 1][j - 1];
}

template <class T>
void matrix<T>::changeelem(int i, int j, T elem) //change the element in row i, column i into elem
{
    if (i > matrix::dim || i <= 0 || j > matrix::dim || j <= 0)
    {
        std::cerr << "The element (" << i << ',' << j << ") does not exist!" << std::endl;
        return;
    }
    matrix::m[i - 1][j - 1] = elem;
}

template <class T>
T matrix<T>::row_max_abs(int i, int b, int e) const //获取第i行最大绝对值
{
    int maxabs_index = row_maxabs_index(i, b, e);
    if (maxabs_index == 0)
        return 0;
    return m[i - 1][maxabs_index - 1];
}

template <class T>
int matrix<T>::row_maxabs_index(int i, int b, int e) const //获取第i行第b个(包含)到第e个元素(不包含)中最大绝对值元素的列数, 返回值在1~n之间
{
    if (i > matrix::dim || i <= 0)
    {
        std::cerr << "Row " << i << " does not exist!" << std::endl;
        return 0;
    }
    else if (b < 1 || e > dim + 1 || b >= e)
    {
        std::cerr << "There's no element between " << b << " and " << e << "!" << std::endl;
        return 0;
    }
    int maxabs_index = b;
    T maxabs = std::abs(m[i - 1][b - 1]);
    for (int k = b + 1; k < e; k++) // k 表示矩阵真实列号，而非存储指标
        if (std::abs(m[i - 1][k - 1]) > maxabs)
        {
            maxabs_index = k;
            maxabs = std::abs(m[i - 1][k - 1]);
        }
    return maxabs_index;
}

template <class T>
T matrix<T>::col_max_abs(int i, int b, int e) const //获取第i列第b个(包含)到第e个元素(不包含)中最大绝对值
{
    int maxabs_index = col_maxabs_index(i, b, e);
    if (maxabs_index == 0)
        return 0;
    return m[maxabs_index - 1][i - 1];
}

template <class T>
int matrix<T>::col_maxabs_index(int i, int b, int e) const //获取第i列第b个(包含)到第e个元素(不包含)中最大绝对值元素的行数, 返回值在1~n之间
{
    if (i > matrix::dim || i <= 0)
        throw std::out_of_range("The column you select is out of range!");
    else if (b < 1 || e > dim + 1 || b >= e)
    {
        std::cerr << "There's no element between " << b << " and " << e << "!" << std::endl;
        return 0;
    }
    int maxabs_index = b;
    T maxabs = std::abs(m[b - 1][i - 1]);
    for (int k = b + 1; k < e; k++) // k 表示矩阵真实列号，而非存储指标
        if (std::abs(m[k - 1][i - 1]) > maxabs)
        {
            maxabs_index = k;
            maxabs = std::abs(m[k - 1][i - 1]);
        }
    return maxabs_index;
}

template <class T>
void matrix<T>::exch_row(int i, int j) //swap elements in two rows
{
    if (i > matrix::dim || i <= 0 || j > matrix::dim || j <= 0)
        throw std::out_of_range("The row you select is out of range!");
    T *temp;
    temp = matrix::m[i - 1];
    m[i - 1] = m[j - 1];
    m[j - 1] = temp;
}

template <class T>
void matrix<T>::exch_col(int i, int j) //swap elements in two columns
{
    if (i > matrix::dim || i <= 0 || j > matrix::dim || j <= 0)
        throw std::out_of_range("The column you select is out of range!");
    T temp;
    for (int k = 0; k < dim; k++)
    {
        temp = m[k][i - 1];
        m[k][i - 1] = m[k][j - 1];
        m[k][j - 1] = temp;
    }
}

template <class T>
void matrix<T>::partial_mainelem(int i) //部分主元, 找到第i~n行, 第i~n列的子矩阵中部分主元所在行，并把这一行与第i行交换.
{
    if (i > dim || i <= 0)
        throw std::out_of_range("The row you select is out of range!");
    int row = col_maxabs_index(i, i, dim + 1);
    if (row > i)
        exch_row(row, i);
}

template <class T>
void matrix<T>::proportion_mainelem(int i) //比例因子主元, 找到第i~n行, 第i~n列的子矩阵中比例因子主元所在行，并把这一行与第i行交换.
{
    if (i > dim || i <= 0)
        throw std::out_of_range("The row you select is out of range!");
    if (i == dim)
        return;
    int row = i;                  //row是真实矩阵行号, 不是数组存储指标
    for (int k = i; k < dim; k++) //k是数组存储指标, 不是真实矩阵行号
        if (std::abs(m[k][i - 1]) * std::abs(m[row - 1][row_maxabs_index(row, i, dim + 1) - 1]) > std::abs(m[row - 1][i - 1]) * std::abs(m[k][row_maxabs_index(k + 1, i, dim + 1) - 1]))
            row = k + 1;
    exch_row(row, i);
}

template <class T>
bool matrix<T>::elementary_transformation(int i1, T k, int i2, bool is_col, int begin) //矩阵 i1 行(或列)的 k 倍加到 i2 行(列)
{
    if (i1 > dim || i1 <= 0 || i2 > dim || i2 <= 0 || begin > dim || begin <= 0)
        throw std::out_of_range("The index is out of range!");
    for (int j = begin; j <= dim; j++)
    {
        if (is_col)
            m[j - 1][i2 - 1] += k * m[j - 1][i1 - 1];
        else
            m[i2 - 1][j - 1] += k * m[i1 - 1][j - 1];
    }
    return true;
}

template <class T>
void matrix<T>::print(int w) const //print this matrix
{
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
            std::cout << std::setw(w) << m[i][j];
        std::cout << std::endl;
    }
}

template <class T>
vec<int> matrix<T>::LU(int mainelem, T zerolim) //LU 分解，返回一个用向量表示的置换矩阵 P, 使得 PA=LU.
{
    vec<int> p(dim);
    for (int i = 1; i <= dim; i++)
        p.v[i - 1] = i;
    for (int i = 1; i <= dim - 1; i++) //i 是矩阵真实行号, 而不是存储指标
    {
        if (mainelem == 1) //选主元, 部分主元
            partial_mainelem_vec(*this, i, p);
        else //一般设为 2
            proportion_mainelem_vec(*this, i, p);
        if (std::abs(m[i - 1][i - 1]) < std::abs(zerolim)) //如果主元为零，则系数矩阵不可逆
        {
            std::cout << "The matrix is not invertible!" << std::endl;
            vec<int> so;
            return so;
        }
        for (int j = i + 1; j <= dim; j++) //消元, 把第i列下面的都消掉
        {
            m[j - 1][i - 1] /= m[i - 1][i - 1]; //消元系数(倍数), 相当于 L 的元素
            for (int k = i + 1; k <= dim; k++)  //把第j行减掉第i行的一个倍数
                m[j - 1][k - 1] -= m[j - 1][i - 1] * m[i - 1][k - 1];
        }
    }
    return p;
}

template <class T>
uptrimatrix<T> matrix<T>::Schmidt(T zerolim) //用 Schmidt 正交化方法进行 QR 分解, 自己变成 Q, 返回 R
{
    uptrimatrix<T> R(dim);          //R to be returned
    vec<T> norm2(dim);              //用于存储各个 β 的模的平方
    vec<T> *beta = new vec<T>[dim]; //用于存储各个 beta
    for (int j = 1; j <= dim; j++)
    {
        vec<T> alpha_j(*this, j, 1);
        for (int k = 1; k < j; k++)                  //计算 beta_j
            if (norm2.v[k - 1] >= std::abs(zerolim)) //模方大于 0 才需要减
            {
                R.sbm[j - k][k - 1] = f(beta[k - 1], alpha_j) / norm2.v[k - 1];
                elementary_transformation(k, -R.sbm[j - k][k - 1], j, true);
            }
        vec<T> beta_j(*this, j, true);
        beta[j - 1] = beta_j;
        norm2.v[j - 1] = f(beta_j, beta_j);
    }
    for (int i = 1; i <= dim; i++)
        R.sbm[0][i - 1] = 1;
    delete[] beta;
    return R;
}

template <class T>
uptrimatrix<T> matrix<T>::Householder_QR(T zerolim) // QR decomposition by Householder method
{
    matrix<T> Q(dim);
    for (int i = 1; i <= dim; i++) //change Q into identity
        Q.m[i - 1][i - 1] = 1;
    for (int j = 1; j < dim; j++) //对第 j 列进行操作
    {
        T norm = 0; //矩阵第 j 列后面待消去向量(包括第 j 行)的模方
        for (int i = j; i <= dim; i++)
            norm += std::abs(m[i - 1][j - 1]) * std::abs(m[i - 1][j - 1]);
        norm = sqrt(norm);
        vec<T> w(dim);
        w.v[j - 1] = m[j - 1][j - 1] - norm;
        for (int i = j + 1; i <= dim; i++)
            w.v[i - 1] = m[i - 1][j - 1];
        // std::cout << "--------------------- w -----------------------" << std::endl;
        // w.print(15);
        T normw = sqrt(f(w, w));
        if (normw < std::abs(zerolim)) //如果这一列主对角线下面全是 0, 那就不用消了
            continue;
        for (int k = j; k <= dim; k++) //把 w 归一化
            w.v[k - 1] /= normw;
        // std::cout << "--------------------- w -----------------------" << std::endl;
        // w.print(15);
        for (int k = j + 1; k <= dim; k++)
        {
            vec<T> m_k(*this, k, true);
            T neiji = f(w, m_k);
            elementary_transformation(j, -2 * neiji / normw, k, true, j);
            m[j - 1][k - 1] += 2 * neiji / normw * norm; //当然第 j 行还少了一点
        }
        for (int i = j + 1; i <= dim; i++)
            m[i - 1][j - 1] = 0;
        m[j - 1][j - 1] = norm;
        //更新 Q
        vec<T> Qw(Q * w);
        for (int k = 1; k <= dim; k++)
            for (int l = j; l <= dim; l++)
                Q.m[k - 1][l - 1] -= 2 * Qw.v[k - 1] * getconj(w.v[l - 1]);
        // std::cout << "--------------------- m -----------------------" << std::endl;
        // print(15);
        // std::cout << "--------------------- Q -----------------------" << std::endl;
        // Q.print(15);
    }
    uptrimatrix<T> R(*this);
    swap(*this, Q);
    return R;
}

template <class T>
uptrimatrix<T> matrix<T>::Givens_QR(T zerolim) // QR decomposition by Givens transform
{
    matrix<T> Q(dim);
    for (int i = 1; i <= dim; i++) //Change Q into indentity
        Q.m[i - 1][i - 1] = 1;
    for (int j = 1; j < dim; j++) //对第 j 列进行操作, 最后一列就不用了
    {
        for (int i = j + 1; i <= dim; i++) //对第 i 列进行操作
        {
            if (std::abs(m[j - 1][j - 1]) * std::abs(zerolim) > std::abs(m[i - 1][j - 1])) //如果两个元素差太多，就可以认为小的那个是 0
                continue;
            T s = 0, c = 0;
            T norm = sqrt(std::abs(m[j - 1][j - 1]) * std::abs(m[j - 1][j - 1]) + std::abs(m[i - 1][j - 1]) * std::abs(m[i - 1][j - 1]));
            s = getconj(m[i - 1][j - 1]) / norm;
            c = getconj(m[j - 1][j - 1]) / norm;
            if (std::abs(norm) < std::abs(zerolim)) //如果模太小，也不用算了
                continue;
            // = sqrt(std::abs(m[j - 1][j - 1]) * std::abs(m[j - 1][j - 1]) + std::abs(m[i - 1][j - 1]) * std::abs(m[i - 1][j - 1]));
            // s = std::conj(m[i - 1][j - 1]) / norm;
            // c = std::conj(m[j - 1][j - 1]) / norm;
            m[j - 1][j - 1] = norm;
            m[i - 1][j - 1] = 0;
            for (int k = j + 1; k <= dim; k++) //同时修改第 i,j 行第 k 列
                //Givens_change(m[j - 1][k - 1], m[i - 1][k - 1], s, c);
            {
                T m_jk = m[j - 1][k - 1];
                m[j - 1][k - 1] = c * m_jk + s * m[i - 1][k - 1];
                m[i - 1][k - 1] = -getconj(s) * m_jk + getconj(c) * m[i - 1][k - 1];
            }
            for (int k = 1; k <= dim; k++) //同时修改正交(酉)矩阵 Q 的第 i,j 列的第 k 行
                //Givens_change(Q.m[k - 1][j - 1], Q.m[k - 1][i - 1], s, c, false);
            {
                T Q_kj = Q.m[k - 1][j - 1];
                Q.m[k - 1][j - 1] = getconj(c) * Q.m[k - 1][j - 1] + getconj(s) * Q.m[k - 1][i - 1];
                Q.m[k - 1][i - 1] = c * Q.m[k - 1][i - 1] - s * Q_kj;
            }
        }
        // std::cout << "--------------------- m -----------------------" << std::endl;
        // print(15);
        // std::cout << "--------------------- Q -----------------------" << std::endl;
        // Q.print(15);
    }
    uptrimatrix<T> um(*this);
    swap(*this, Q);
    return um;
}

template <class T>
matrix<T> matrix<T>::upper_Hessenberg(T zerolim) //正交相似于一个上海森堡矩阵 (*this=Q^THQ), 返回 Q
{
    matrix<T> Q(dim);
    for (int i = 1; i <= dim; i++) //Change Q into indentity
        Q.m[i - 1][i - 1] = 1;
    for (int j = 1; j < dim - 1; j++) //消第 j 列，最后两列就不用消了
    {
        T norm = 0; //矩阵第 j 列后面待消去向量(包括第 j+1 行)的模
        for (int i = j + 1; i <= dim; i++)
            norm += std::abs(m[i - 1][j - 1]) * std::abs(m[i - 1][j - 1]);
        norm = sqrt(norm);
        vec<T> w(dim);
        w.v[j] = m[j][j - 1] - norm;
        for (int i = j + 2; i <= dim; i++)
            w.v[i - 1] = m[i - 1][j - 1];
        T normw = sqrt(f(w, w));
        if (normw < std::abs(zerolim)) //如果这一列主对角线下面全是 0, 那就不用消了
            continue;
        for (int k = j + 1; k <= dim; k++) //把 w 归一化
            w.v[k - 1] /= normw;
        // std::cout << "--------------------- w -----------------------" << std::endl;
        // w.print(15);
        for (int k = j + 1; k <= dim; k++) //对第 k 列进行操作
        {
            vec<T> m_k(*this, k, true); //原矩阵的第 k 列
            T inner_product = f(w, m_k);
            elementary_transformation(j, -2 * inner_product / normw, k, true, j + 1);
            m[j][k - 1] += 2 * inner_product / normw * norm; //当然第 (j+1) 行还少了一点
        }
        for (int i = j + 2; i <= dim; i++) //把第 j 列第 (j+2) 行及以后的行都改为 0
            m[i - 1][j - 1] = 0;
        m[j][j - 1] = norm;
        for (int i = 1; i <= dim; i++) //右乘一个 Householder 矩阵, 对第 i 行变换
        {
            T inner_product = 0;
            for (int k = j + 1; k <= dim; k++)
                inner_product += m[i - 1][k - 1] * w.v[k - 1];
            for (int k = j + 1; k <= dim; k++)
                m[i - 1][k - 1] -= 2 * inner_product * getconj(w.v[k - 1]);
        }
        //更新 Q
        for (int k = 1; k <= dim; k++)
        {
            T inner_product = 0;
            for (int i = j + 1; i <= dim; i++)
                inner_product += getconj(w.v[i - 1]) * Q.m[i - 1][k - 1];
            for (int i = j + 1; i <= dim; i++)
                Q.m[i - 1][k - 1] -= 2 * inner_product * w.v[i - 1];
        }
    }
    return Q;
}

template <class T>
vec<T> matrix<T>::linear_f(vec<T> &b, int mainelem, T zerolim) //返回以此矩阵为系数矩阵，以b为非齐次项的线性方程组的解
{
    vec<int> p(LU());
    if (p.dim == 0)
    {
        vec<T> so;
        return so;
    }
    //求出对应的上三角矩阵和下三角矩阵
    lowtrimatrix<T> l(*this);
    uptrimatrix<T> u(*this);
    for (int i = 0; i < dim; i++)
        l.sbm[0][i] = 1;
    vec<T> b2(dim);
    for (int i = 1; i <= dim; i++)
        b2.v[i - 1] = b.v[p.v[i - 1] - 1];
    vec<T> y(l.linear_f(b2));
    return u.linear_f(y);
}

//*************** Other Functions ********* SWAP ************
template <class Type>
inline void swap(vec<Type> &lhs, vec<Type> &rhs) //swap two vectors
{
    if (&lhs != &rhs)
    {
        using std::swap;
        swap(lhs.dim, rhs.dim);
        swap(lhs.v, rhs.v);
        swap(lhs.IsCol, rhs.IsCol);
    }
}

template <class Type>
inline void swap(matrix<Type> &lhs, matrix<Type> &rhs) //swap two matrices
{
    if (&lhs != &rhs)
    {
        using std::swap;
        swap(lhs.dim, rhs.dim);
        swap(lhs.m, rhs.m);
    }
}

template <class Type>
inline void swap(sparse_matrix<Type> &lhs, sparse_matrix<Type> &rhs) //swap two matrices
{
    if (&lhs != &rhs)
    {
        using std::swap;
        swap(lhs.dim, rhs.dim);
        swap(lhs.spm, rhs.spm);
    }
}

template <class Type>
inline void swap(sym_band_matrix<Type> &lhs, sym_band_matrix<Type> &rhs) //swap two symmetric banded matrices
{
    if (&lhs != &rhs)
    {
        using std::swap;
        swap(lhs.dim, rhs.dim);
        swap(lhs.sbm, rhs.sbm);
        swap(lhs.halfbandwidth, rhs.halfbandwidth);
    }
}

template <class Type>
inline void swap(uptrimatrix<Type> &lhs, uptrimatrix<Type> &rhs) //swap two upper triangular matrices
{
    if (&lhs != &rhs)
    {
        using std::swap;
        swap(lhs.dim, rhs.dim);
        swap(lhs.sbm, rhs.sbm);
        swap(lhs.halfbandwidth, rhs.halfbandwidth);
        swap(lhs.IsT, rhs.IsT);
    }
}

template <class Type>
inline void swap(lowtrimatrix<Type> &lhs, lowtrimatrix<Type> &rhs) //swap two lower triangular matrices
{
    if (&lhs != &rhs)
    {
        using std::swap;
        swap(lhs.dim, rhs.dim);
        swap(lhs.sbm, rhs.sbm);
        swap(lhs.halfbandwidth, rhs.halfbandwidth);
        swap(lhs.IsT, rhs.IsT);
    }
}

template <class Type>
inline void swap(tridiagmatrix<Type> &lhs, tridiagmatrix<Type> &rhs) //swap two tridiagonal matrices
{
    if (&lhs != &rhs)
    {
        using std::swap;
        swap(lhs.dim, rhs.dim);
        swap(lhs.tdm, rhs.tdm);
    }
}

//*************** Other Functions ********* operator overloading ************
template <class T>
vec<T> operator+(vec<T> lhs, const vec<T> &rhs) //overload operator vec + 
{
    lhs += rhs;
    return lhs;
}

template <class T>
vec<T> operator-(vec<T> lhs, const vec<T> &rhs) //overload operator vec - 
{
    lhs -= rhs;
    return lhs;
}

template <class T>
vec<T> operator*(const T &x, vec<T> ve) //overload operator number multiplication of vec
{
    for (int i = 0; i < ve.dim; i++)
        ve.v[i] *= x;
    return ve;
}

template <class Type>
vec<Type> operator*(const sym_band_matrix<Type> &A, const vec<Type> &v) //a symmetric matrix multiplied by a vector
{
    if (v.dim != A.dim)
        throw std::invalid_argument("The dimensions of the vector and the matrix are different!");
    vec<Type> vso(v.dim);
    for (int i = 1; i <= v.dim; i++)
    {
        Type tempsum = 0;
        if (i <= A.halfbandwidth + 1)
            for (int j = 1; j <= i; j++)
                tempsum += A.sbm[i - j][j - 1] * v.v[j - 1];
        else
            for (int j = i - A.halfbandwidth; j <= i; j++)
                tempsum += A.sbm[i - j][j - 1] * v.v[j - 1];
        if (i <= v.dim - A.halfbandwidth)
            for (int j = i + 1; j <= i + A.halfbandwidth; j++)
                tempsum += A.sbm[j - i][i - 1] * v.v[j - 1];
        else
            for (int j = i + 1; j <= v.dim; j++)
                tempsum += A.sbm[j - i][i - 1] * v.v[j - 1];
        vso.v[i - 1] = tempsum;
    }
    return vso;
}

template <class Type>
uptrimatrix<Type> operator*(uptrimatrix<Type> um, const diagmatrix<Type> &dm) //a upper triangular matrix multiplied by a diagonal matrix
{
    if (um.dim != dm.dim)
        throw std::invalid_argument("The dimensions of the matrices are different!");
    for (int i = 1; i <= um.dim; i++)
        for (int j = i; j <= um.dim; j++)
            um.sbm[j - i][i - 1] *= dm.sbm[0][j - 1];
    return um;
}

template <class Type>
lowtrimatrix<Type> operator*(lowtrimatrix<Type> lm, const diagmatrix<Type> &dm) //a lower triangular matrix multiplied by a diagonal matrix
{
    if (lm.dim != dm.dim)
        throw std::invalid_argument("The dimensions of the matrices are different!");
    for (int i = 1; i <= lm.dim; i++)
        for (int j = 1; j <= i; j++)
            lm.sbm[i - j][j - 1] *= dm.sbm[0][j - 1];
    return lm;
}

template <class Type>
matrix<Type> operator*(const matrix<Type> &ma, const uptrimatrix<Type> &up) //a matrix multiplied by a upper triangular matrix
{
    if (ma.dim != up.dim)
        throw std::invalid_argument("The dimensions of the two matrices are different!");
    matrix<Type> so(ma.dim);
    for (int i = 1; i <= ma.dim; i++)
        for (int j = 1; j <= ma.dim; j++)
            for (int k = std::max(1, j - up.halfbandwidth); k <= j; k++)
                so.m[i - 1][j - 1] += ma.m[i - 1][k - 1] * up.sbm[j - k][k - 1];
    return so;
}

template <class T>
vec<T> operator*(const sparse_matrix<T> &A, const vec<T> &ve) //a sparse matrix multiplied by a vector
{
    if (ve.dim != A.dim)
        throw std::invalid_argument("The dimensions of the vector and the matrix are different!");
    vec<T> b(A.dim);
    for (int i = 1; i <= A.dim; i++)
        for (typename std::list<spmnode<T>>::iterator it = A.spm[i - 1].begin(); it != A.spm[i - 1].end(); it++)
            b.v[i - 1] += it->elem * ve.v[it->index - 1];
    return b;
}

template <class T>
matrix<T> operator*(const matrix<T> &A, const matrix<T> &B) //a matrix multiplied by a matrix
{
    if (B.dim != A.dim)
        throw std::invalid_argument("The dimensions of the two matrices are different!");
    matrix<T> product(A.dim);
    for (int i = 1; i <= A.dim; i++)
        for (int j = 1; j <= A.dim; j++)
            for (int k = 1; k <= A.dim; k++)
                product.m[i - 1][j - 1] += A.m[i - 1][k - 1] * B.m[k - 1][j - 1];
    return product;
}

template <class Type>
matrix<Type> operator*(const uptrimatrix<Type> &um, const matrix<Type> &ma) //a upper triangular matrix multiplied by a matrix
{
    if (ma.dim != um.dim)
        throw std::invalid_argument("The dimensions of the two matrices are different!");
    matrix<Type> product(ma.dim);
    for (int i = 1; i <= ma.dim; i++)
        for (int j = 1; j <= ma.dim; j++)
            for (int k = i; k <= std::min(ma.dim, i + um.halfbandwidth); k++)
                product.m[i - 1][j - 1] += um.sbm[k - i][i - 1] * ma.m[k - 1][j - 1];
    return product;
}

template <class Type>
matrix<Type> operator*(const sym_band_matrix<Type> &sm, const matrix<Type> &ma) //a symmetric matrix multiplied by an ordinary matrix
{
    if (ma.dim != sm.dim)
        throw std::invalid_argument("The dimensions of the two matrices are different!");
    matrix<Type> product(ma.dim);
    for (int k = 0; k <= sm.halfbandwidth; k++) //先把对称矩阵上半(包括对角线)算完
        for (int l = 0; l < sm.dim - k; l++)
            for (int t = 0; t < sm.dim; t++) //product 的第(l+1)行第(t+1)列
                product.m[l][t] += sm.sbm[k][l] * ma.m[k + l][t];
    for (int k = 1; k <= sm.halfbandwidth; k++) //再算左下半(不包括对角线)
        for (int l = 0; l < sm.dim - k; l++)
            for (int t = 0; t < sm.dim; t++) //product 的第(k+l+1)行第(t+1)列
                product.m[k + l][t] += sm.sbm[k][l] * ma.m[l][t];
    return product;
}

template <class Type>
matrix<Type> operator*(const matrix<Type> &ma, const sym_band_matrix<Type> &sm) //a matrix multiplied by a symmetric matrix
{
    if (ma.dim != sm.dim)
        throw std::invalid_argument("The dimensions of the two matrices are different!");
    matrix<Type> product(ma.dim);
    for (int k = 0; k <= sm.halfbandwidth; k++) //先把对称矩阵上半(包括对角线)算完
        for (int l = 0; l < sm.dim - k; l++)
            for (int t = 0; t < sm.dim; t++) //product 的第(t+1)行第(k+l+1)列
                product.m[t][l + k] += sm.sbm[k][l] * ma.m[t][l];
    for (int k = 1; k <= sm.halfbandwidth; k++) //再算左下半(不包括对角线)
        for (int l = 0; l < sm.dim - k; l++)
            for (int t = 0; t < sm.dim; t++) //product 的第(t+1)行第(l+1)列
                product.m[t][l] += ma.m[t][k + l] * sm.sbm[k][l];
    return product;
}

template <class Type>
vec<Type> operator*(const matrix<Type> &A, const vec<Type> &ve) //a matrix multiplied by a vector
{
    if (ve.dim != A.dim)
        throw std::invalid_argument("The dimensions of the vector and the matrix are different!");
    vec<Type> Av(A.dim);
    for (int i = 1; i <= A.dim; i++)
        for (int j = 1; j <= A.dim; j++)
            Av.v[i - 1] += A.m[i - 1][j - 1] * ve.v[j - 1];
    return Av;
}

//*************** Other Functions ********* Inner Product ************
template <class Type>
Type f(const vec<Type> &v1, const vec<Type> &v2) //Inner product of two vectors
{
    if (v1.dim != v2.dim)
        throw std::invalid_argument("The dimensions of the two vectors are different!");
    Type sum = 0;
    for (int i = 0; i < v1.dim; i++)
        sum += getconj(v1.v[i]) * v2.v[i];
    return sum;
}

template <class Type>
Type f(const sym_band_matrix<Type> &A, const vec<Type> &p) //Inner product: (p^T)Ap
{
    if (p.dim != A.dim)
        throw std::invalid_argument("The dimensions of the vector and the matrix are different!");
    Type tempsum = 0;
    for (int j = 0; j < A.dim; j++)
        tempsum += p.v[j] * p.v[j] * A.sbm[0][j];
    for (int i = 1; i <= A.halfbandwidth; i++)
        for (int j = 0; j < A.dim - i; j++)
            tempsum += 2 * A.sbm[i][j] * p.v[j] * p.v[i + j];
    return tempsum;
}

template <class Type>
bool cgm(vec<Type> &x, const vec<Type> &b, const sym_band_matrix<Type> &A, int stepnum, Type zerolim) //Conjugate Gradient Algorithm
{
    if (x.dim != A.dim || x.dim != b.dim || A.dim != b.dim) //if the dimensions are different
        return false;
    if (stepnum > A.dim)
        stepnum = A.dim;    //最多迭代 A.dim 步
    vec<Type> r(A * x - b); //残差
    vec<Type> p(r);         //搜索方向
    Type rTrlast = f(r, r);
    Type alpha, beta, rTr;
    Type pTAp = f(A, p);
    for (int stepdone = 0; stepdone <= stepnum; stepdone++)
    {
        if (std::abs(rTrlast) < std::abs(zerolim) || std::abs(pTAp) < std::abs(zerolim))
            return true;
        alpha = rTrlast / pTAp;
        x -= alpha * p; //更新解
        r = A * x - b;  //更新残差
        rTr = f(r, r);
        beta = rTr / rTrlast;
        p *= beta;
        p += r; //更新搜索方向
        rTrlast = rTr;
        pTAp = f(A, p);
    }
    return false;
}

template <class T1, class T2>
bool partial_mainelem_vec(matrix<T1> &ma, int i, vec<T2> &b) //部分主元, 找到第i~n行, 第i~n列的子矩阵中部分主元所在行，并把这一行与第i行交换, 同时交换向量 b 相应元素
{
    if (i > ma.getdim() || i <= 0)
        throw std::out_of_range("The row you select is out of range!");
    if (ma.getdim() != b.getdim())
        throw std::invalid_argument("The dimensions of the vector and the matrix are different!");
    int row = ma.col_maxabs_index(i, i, ma.getdim() + 1);
    if (row > i)
    {
        ma.exch_row(row, i);
        b.exch_elem(row, i);
    }
    return true;
}

template <class T1, class T2>
bool proportion_mainelem_vec(matrix<T1> &ma, int i, vec<T2> &b) //部分主元, 找到第i~n行, 第i~n列的子矩阵中部分主元所在行，并把这一行与第i行交换, 同时交换向量 b 相应元素
{
    if (i > ma.getdim() || i <= 0)
    {
        std::cerr << "Row " << i << " does not exist!" << std::endl;
        return false;
    }
    if (ma.getdim() != b.getdim())
    {
        std::cerr << "The dimension of the matrix and the vector is  different!" << std::endl;
        return false;
    }
    if (i == ma.getdim())
        return true;
    int row = i;                          //row是真实矩阵行号, 不是数组存储指标
    for (int k = i; k < ma.getdim(); k++) //k是数组存储指标, 不是真实矩阵行号
        if (std::abs(ma.getelem(k + 1, i)) * std::abs(ma.getelem(row, ma.row_maxabs_index(row, i, ma.getdim() + 1))) > std::abs(ma.getelem(row, i)) * std::abs(ma.getelem(k + 1, ma.row_maxabs_index(k + 1, i, ma.getdim() + 1))))
            row = k + 1;
    ma.exch_row(row, i);
    b.exch_elem(row, i);
    return true;
}

#endif //MATRIX_H_