#include <iostream>
using namespace std;

int main()
{
    double pow[11]; //用来存放指数结果的数组
    for (int i = 0; i < 11; i++)  //直接展开法求 e^(-x)
    {
        pow[i] = 1.0;
        double x = i*10 ;

        for (int n = 1; n < 3.0 * x; n++) //为了保证收敛, n 至少大于 e*x
        {
            double addnum = 1.0;
            if (n % 2) // 如果 n 是奇数, 结果为负
                addnum = -1.0;
            for (int k = 1; k <= n; k++)  //考虑到先把方幂求出来再除以阶乘可能会越界, 所以我乘一个之后就除以一个
            {
                addnum /= k;
                addnum *= x;
            }
            pow[i]+=addnum;
        }
        cout<<pow[i]<<" ";
    }
    cout<<endl;

    for (int i = 0; i <= 10; i++) //递推计算 e^(-x)
    {
        pow[i]=1.0;
        double x=i*10.0;
        double s=1.0;
        double s_last=1.0;
        for(int n=1;n<3.0*x;n++)  
        {
            s=-s_last/n*x;
            pow[i]+=s;
            s_last=s;
        }
        cout<<pow[i]<<" ";
    }
    cout<<endl;

     for (int i = 0; i <= 10; i++)//先计算e^x, 再取倒数
    {
        pow[i]=1.0;
        double x=i*10.0;
        double s=1.0;
        double s_last=1.0;
        for(int n=1;n<3.0*x;n++)  
        {
            s=s_last/n*x;
            pow[i]+=s;
            s_last=s;
        }
        pow[i]=1.0/pow[i];
        cout<<pow[i]<<" ";
    }
    cout<<endl;
    return 0;
}
