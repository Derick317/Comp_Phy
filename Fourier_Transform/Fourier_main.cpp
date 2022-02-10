#include<iostream>
#include"Fourier.h"
using namespace std;

int main()
{
    complex<double> x[4]={1.0,2.0,-1.0,3.0};
    complex<double> z[4];
    complex<double> temp[4];
    fft(x,z,temp,4);
    for(int i=0;i<4;i++)
    cout<<z[i]<<endl;
    return 0;
}