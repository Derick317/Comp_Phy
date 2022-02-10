import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return 1.0 / (1.0 + 25.0 * x * x)

def l(x, x0, y0, nudenum):
    d = []  # 分母
    nu = []  # 分子
    for k in range(nudenum):
        d.append(1.0)
        nu.append(1.0)
        for j in range(nudenum):
            if (j != k):
                d[k] = d[k]*(x0[k]-x0[j])
                nu[k] = nu[k]*(x-x0[j])
    lx = 0.0
    for k in range(nudenum):
        lx += y0[k]*nu[k]/d[k]
    return lx

x0 = []
y0 = []
for i in range(21):
    x0.append(-1+i*0.1)
    y0.append(f(x0[i]))

#打表
x1 = []
y1 = []
p1 = []
for i in range(41):
    x1.append(-1+0.05*i)
    y1.append(f(x1[i]))
    p1.append(l(x1[i], x0, y0, 21))
    print(x1[i], " ", y1[i], " ", p1[i])

#画图
x2=[]
y2=[]
p2=[]
for i in range(501):
    x2.append(-1+0.004*i)
    y2.append(f(x2[i]))
    p2.append(l(x2[i], x0, y0, 21))
plt.title('Runge Function: f(x)=1/(1+25x²)') #写上图题
plt.xlabel('x') #为x轴命名为“x”
plt.ylabel('y') #为y轴命名为“y”
plt.xlim(-1,1) #设置x轴的范围为[0,1]
plt.ylim(-2,2) #同上
plt.xticks([-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1]) #设置x轴刻度
plt.yticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]) #设置y轴刻度
plt.tick_params(labelsize = 5) #设置刻度字号
plt.plot(x2,y2) #第一个data表示选取data为数据集，第二个是函数，data的平方
plt.plot(x2,p2) #同上
plt.legend(['y=f(x)','y=P₂₀(x)']) #打出图例
plt.show() #显示图形
