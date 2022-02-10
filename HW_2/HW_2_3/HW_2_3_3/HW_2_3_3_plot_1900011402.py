import numpy as np
import matplotlib.pyplot as plt

with open("Runge_Spline.txt") as f:
    data=f.readlines()
x=[]
y=[]  #真实函数值
s=[]  #拟合函数值
for i in range(len(data)):
    item=data[i].split()
    x.append(float(item[0]))
    y.append(float(item[1]))
    s.append(float(item[2]))

plt.title('Spline Interpolation of Runge Function: f(x)=1/(1+25x²)') #写上图题
plt.xlabel('x') #为x轴命名为“x”
plt.ylabel('y') #为y轴命名为“y”
plt.xlim(-1,1) #设置x轴的范围为[0,1]
plt.ylim(0,1.2) #同上
plt.xticks([-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1]) #设置x轴刻度
plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2]) #设置y轴刻度
plt.tick_params(labelsize = 5) #设置刻度字号
plt.plot(x,y,linewidth=0.5) #
plt.plot(x,s,linewidth=0.5,alpha=0.8) #
plt.legend(['Runge Function','Spline Interpolation']) #打出图例
plt.show() #显示图形