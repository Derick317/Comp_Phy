//直接展开法求 e^(-x)
#include <iostream>
using namespace std;

int main()
{
    double pow[11]; //用来存放指数结果的数组
    for (int i = 0; i < 11; i++)
    {
        pow[i] = 1.0;
        double x = i*10 ;

        for (int n = 1; n < 3.0 * x; n++) //为了保证收敛, n 至少大于 e*x
        {
            double addnum = 1.0;
            if (n % 2) // 如果 n 是奇数, 结果为负
                addnum = -1.0;
            for (int k = 1; k <= n; k++)
            {
                addnum /= k;
                addnum *= x;
            }
            pow[i]+=addnum;
        }
        cout<<pow[i]<<endl;
    }

    return 0;
}