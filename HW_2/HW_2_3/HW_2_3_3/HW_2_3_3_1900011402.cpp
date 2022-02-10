#include"../../../Interpolation/interpolation.h"
#include<fstream>
using namespace std;

inline double f(double x)
{
    return 1.0/(1.0+25.0*x*x);
}

int main()
{
    double x[21],y[21];
    for(int i=0;i<21;i++)
    {
        x[i]=-1.0+i*0.1;
        y[i]=f(x[i]);
    }
    cubic_spline<double> s(x,y,21,'2');
    // for(int i=0;i<41;i++) //打表
    // {
    //     double x1=-1.0+i*0.05;
    //     cout<<x1<<' '<<f(x1)<<' '<<setprecision(8)<<s.f(x1)<<endl;
    // }
    ofstream fout;
    fout.open("Runge_Spline.txt");
    if(!fout)
    return 1;
    fout.clear();
    for(int i=0;i<501;i++) //把数据记录下来
    {
        double x1=-1.0+i*0.004;
        fout<<x1<<' '<<f(x1)<<' '<<setprecision(8)<<s.f(x1)<<endl;
    }
    fout.close();
    return 0;
}