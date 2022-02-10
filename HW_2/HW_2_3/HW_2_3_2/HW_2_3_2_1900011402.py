import numpy as np
import matplotlib.pyplot as plt
import math

def f(x):
    return 1/(1+25*x*x)

def T(k,x):
    return math.cos(k*math.acos(x))

def C(c,x): #c是系数列表
    y=0
    for i in range(len(c)):
        y+=c[i]*T(i,x)
    return y

x=[] #Chebyshev节点
for i in range (20):
    x.append(math.cos(math.pi*(i+0.5)/20))

c=[] #系数
c0=0
for i in range(20):
    c0+=f(x[i])
c0/=20
c.append(c0)
for k in range (1,20): #对剩下19个系数操作, 因为它们都只除以10
    ck=0
    for i in range(20):
        ck+=T(k,x[i])*f(x[i])
    ck/=10
    c.append(ck)

#列表
print(-1.0,' ',f(-1),' ',C(c,-1))
print(x[19],' ',f(x[19]),' ',C(c,x[19]))
for i in range (19,0,-1):
    # print(i)
    mid=(x[i]+x[i-1])/2
    print(mid,' ',f(mid),' ',C(c,mid))
    print(x[i-1],' ',f(x[i-1]),' ',C(c,x[i-1]))
print(1.0,' ',f(1),' ',C(c,1))

# x1=[]  #用于画图的自变量
# y1=[]  #用于画图的真实函数值
# t1=[]  #用于画图的拟合值

# #画图
# for i in range(501):
#     x1.append(-1+0.004*i)
#     y1.append(f(x1[i]))
#     t1.append(C(c,x1[i]))
# plt.title('Chebyshev Approxiamtin of Runge Function: f(x)=1/(1+25x²)') #写上图题
# plt.xlabel('x') #为x轴命名为“x”
# plt.ylabel('y') #为y轴命名为“y”
# plt.xlim(-1,1) #设置x轴的范围为[0,1]
# plt.ylim(0,1.2) #同上
# plt.xticks([-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1]) #设置x轴刻度
# plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2]) #设置y轴刻度
# plt.tick_params(labelsize = 5) #设置刻度字号
# plt.plot(x1,y1,linewidth=1) #第一个data表示选取data为数据集，第二个是函数，data的平方
# plt.plot(x1,t1,linewidth=1,alpha=0.8) #同上
# plt.legend(['Runge Function','Chebyshev Approximation']) #打出图例
# plt.show() #显示图形