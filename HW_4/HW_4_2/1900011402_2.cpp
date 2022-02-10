#include "../../Linear_Funtions/matrix.h"
using namespace std;

int main()
{
    sym_band_matrix<double> A(10, 9);
    vec<double> q(10), z(10);
    double nu = 0;
    for (int i = 1; i < 10; i++)
    {
        A.changeelem(i, i, 2);
        A.changeelem(i, i + 1, -1);
        //A.changeelem(i + 1, i, -1);
    }
    A.changeelem(10, 10, 2);
    A.changeelem(1, 10, -1);
    //A.changeelem(10, 1, -1);
    q.changeelem(1, 1);

    //cout << f(q, q) << endl;

    for (int i = 1; i <= 60; i++)
    {
        z = A * q;
        q = (1.0 / sqrt(f(z, z))) * z;
        if(i%5==0)
        cout << i << " " << f(q, A * q) << endl;
    }
    nu = f(q, A * q);
    q.trans();

    cout << "------------- q ----------------" << endl;
    q.print(10);
    return 0;
}