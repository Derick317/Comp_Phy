#include "Integral.h"
using namespace std;

double f(double x)
{
    return x * x;
}

int main()
{
    cout<<trapezoid(f,0.0,3.0,20)<<endl;
    return 0;
}