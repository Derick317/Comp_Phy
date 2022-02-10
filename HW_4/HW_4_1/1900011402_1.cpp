#include "../../Linear_Funtions/matrix.h"
using namespace std;

int main()
{
    matrix<double> m(4), Q(4);
    sym_band_matrix<double> sm1(4, 3), sm2(4, 1);
    uptrimatrix<double> R(4, 3);
    for (int i = 1; i < 5; i++)
    {
        m.changeelem(i, i, i);
        sm1.changeelem(i, i, i);
        sm2.changeelem(i, i, i);
    }
    for (int i = 1; i < 4; i++)
    {
        m.changeelem(i, i + 1, -1);
        m.changeelem(i + 1, i, -1);
        sm1.changeelem(i + 1, i, -1);
        sm2.changeelem(i + 1, i, -1);
    }
    m.print(15);
    // QR 算法
    cout << "---------- QR -----------" << endl;
    for (int i = 0; i < 50; i++)
    {
        R = m.Givens_QR();
        m = R * m;
        if (!((i + 1) % 5))
        {
            cout << "---------- " << i + 1 << " -----------" << endl;
            m.print(15);
        }
    }
    //Jacobi 算法, 每次找绝对值最大非对角元旋转
    cout << "---------- Jacobi -----------" << endl;
    for (int i = 0; i < 30; i++)
    {
        int r = 2, c = 1; // 最大非对角元行, 列指标
        for (int k = 2; k < 5; k++)
            for (int l = 1; l < k; l++)
                if (abs(sm1.getelem(k, l)) > abs(sm1.getelem(r, c)))
                {
                    r = k;
                    c = l;
                }
        sm1.Jacobi(r, c);
        if (!((i + 1) % 5))
        {
            cout << "---------- " << i + 1 << " -----------" << endl;
            sm1.print(15);
        }
    }
    // Strum 算法
    cout << "---------- Strum -----------" << endl;
    for (int i = 1; i < 5; i++)
        cout << sm2.Sturm(i, 30,0.00001) << endl;
    return 0;
}

// int main()
// {
//     matrix<double> Q, QT(4), m(4);
//     sym_band_matrix<double> sm(4, 3);
//     for (int i = 1; i < 5; i++)d
//     {
//         sm.changeelem(i, i, i);
//     }
//     for (int i = 1; i < 4; i++)
//     {
//         sm.changeelem(i + 1, i, -1);
//     }

//     Q = sm.Jacobi(2, 1);
//     for (int i = 1; i < 5; i++)
//         for (int j = 1; j < 5; j++)
//             QT.changeelem(i, j, Q.getelem(j, i));
//     cout << "---------- Q -----------" << endl;
//     Q.print(10);
//     cout << "---------- sm -----------" << endl;
//     sm.print(10);
//     cout << "---------- sm0 -----------" << endl;
//     (QT * sm * Q).print(10);
//     return 0;
// }