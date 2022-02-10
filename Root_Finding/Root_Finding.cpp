#include <iostream>
#include "root.h"
//#include "..\IsZero.h"
#include <cmath>
#include <complex>
using namespace std;

double f(double x)
{
    return (x+3)*(x-1)*(x-1);
}
double phi(double x)
{
    return -1+5/x-3/x/x;
}

int main()
{
    double x=2.0;
    cout<<Iterative_Aitken(phi,x,100,f,0.00000001)<<endl;

    return 0;
}
