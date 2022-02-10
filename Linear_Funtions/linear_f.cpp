#include <iostream>
#include "matrix.h"
using namespace std;

int main()
{
	vec<complex<double>> q(5);
	vec<double> a(5);
	q.changeelem(1,2);
	a.changeelem(3,2);
	q.print();
	cout<<f(q,q)<<endl;
	return 0;
}