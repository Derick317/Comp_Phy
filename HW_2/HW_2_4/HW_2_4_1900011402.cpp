#include "../../Interpolation/interpolation.h"
#include <fstream>
using namespace std;

int main()
{
    double t[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    double xt[9] = {0, 0.207106781, 0, -1.207106781, -2, -1.207106781, 0, 0.207106781, 0};
    double yt[9] = {0, 0.207106781, 1, 1.207106781, 0, -1.207106781, -1, -0.207106781, 0};
    cubic_spline<double> sx(t, xt, 9, 'p');
    cubic_spline<double> sy(t, yt, 9, 'p');
    //输出数据到一个文件, 每隔 0.05 读一个点
    ofstream fout;                //
    fout.open("Heart_spline_date.txt"); //打开文件
    if (!fout)
        return 1;
    fout.clear();
    for (int i = 0; i <= 240; i++)
    {
        double t0 = i / 30.0;
        fout << sx.f(t0) << ' ' << sy.f(t0) << endl; //把数存进文件
    }
    fout.close(); //关闭文件
    cout<<"Sx(t)"<<endl;
    for (int i = 1; i < 9; i++)//输出 Sx(t)的系数
        cout << sx.getcoeff(i, 3) << "t^3+" << sx.getcoeff(i, 2) << "t^2+" << sx.getcoeff(i, 1) << "t+" << sx.getcoeff(i, 3) << " ,t\\in [" << i - 1 << ',' << i << "]\\\\" << endl;
    cout<<"Sy(t)"<<endl;
    for (int i = 1; i < 9; i++)//输出 Sx(t)的系数
        cout << sy.getcoeff(i, 3) << "t^3+" << sy.getcoeff(i, 2) << "t^2+" << sy.getcoeff(i, 1) << "t+" << sy.getcoeff(i, 3) << " ,t\\in [" << i - 1 << ',' << i << "]\\\\" << endl;

    return 0;
}