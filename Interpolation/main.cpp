#include "interpolation.h"
using namespace std;

int main()
{
    double x[] = {0, 0.6, 0.9};
    double y[] = {1, 0.9358968236779, 0.689498432951747 };
    cubic_spline<double> s(x, y, 3, '2');
    for(int i=1;i<3;i++)
    for(int d=3;d>=0;d--)
    cout<<s.getcoeff(i,d)<<endl;
    
    return 0;
}