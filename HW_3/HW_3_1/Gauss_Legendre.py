import scipy.special
import math

def f(x):
    return (math.exp(-x))/x

N=10
[xi,wi] = scipy.special.roots_legendre(N)
xi=49.5*xi+50.5
sum=0
for i in range(N):
    sum+=f(xi[i])*wi[i]
sum*=49.5
print(sum)

N=100
[xi,wi] = scipy.special.roots_legendre(N)
xi=49.5*xi+50.5
sum=0
for i in range(N):
    sum+=f(xi[i])*wi[i]
sum*=49.5
print(sum)
