#include "../../../Root_Finding/root.h"
using namespace std;

double f1(double x)
{
    return x - 2 * sin(x);
}

double df1(double x)
{
    return 1 - 2 * cos(x);
}

double f2(double x)
{
    return (x - 2 * sin(x)) * (x - 2 * sin(x));
}

double df2(double x)
{
    return 2 * (x - 2 * sin(x)) * (1 - 2 * cos(x));
}

int main()
{
    // double root = Secant(f1, 1.5, 1.5-f1(1.5)/df1(1.5), 100, 0.00001);
    // cout << std::setprecision(8) << root << endl;
    double root = Newton(f1,df1, 1.5,100, 0.00001);
     cout << std::setprecision(8) << root << endl;
    return 0;
}
