#include "../../Integral/Integral.h"
#define _E_ 2.718281828459045
using namespace std;

double f(double x)
{
    return pow(_E_, -x) / x;
}

int main()
{
    cout << trapezoid(f, 1.0, 100.0, 10) << endl;
    cout << trapezoid(f, 1.0, 100.0, 100) << endl;
    cout << trapezoid(f, 1.0, 100.0, 1000) << endl;
    cout << Simpson(f, 1.0, 100.0, 10) << endl;
    cout << Simpson(f, 1.0, 100.0, 100) << endl;
    cout << Simpson(f, 1.0, 100.0, 1000) << endl;
    return 0;
}