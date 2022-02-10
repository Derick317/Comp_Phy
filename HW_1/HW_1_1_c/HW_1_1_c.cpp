//先计算e^x, 再取倒数
#include <iostream>
using namespace std;

int main()
{
    double pow[11];
    for (int i = 0; i <= 10; i++)
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
        cout<<pow[i]<<endl;
    }

    return 0;
}