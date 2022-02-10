import numpy as np
import matplotlib.pyplot as plt
import math

with open("Heart_spline_date.txt") as f:
    data = f.readlines()
x = []  # 拟合心形线
y = []
for i in range(len(data)):
    item = data[i].split()
    x.append(float(item[0]))
    y.append(float(item[1]))

x0 = []  # 用来存放真实心形线, 共分 400 个点
y0 = []
for i in range(401):
    phi = math.pi*i/200
    x0.append((1-math.cos(phi))*math.cos(phi))
    y0.append((1-math.cos(phi))*math.sin(phi))

y1 = [0, 0.207106781, 1, 1.207106781, 0, -1.207106781, -1, -0.207106781, 0]  # 插值节点
x1 = [0, 0.207106781, 0, -1.207106781, -2, -1.207106781, 0, 0.207106781, 0]

plt.title('Heart Line')  # 写上图题
plt.xlabel('x')  # 为x轴命名为“x”
plt.ylabel('y')  # 为y轴命名为“y”
plt.xlim(-2.5, 0.5)  # 设置x轴的范围为[0,1]
plt.ylim(-1.5, 1.5)  # 同上
plt.xticks([-2.5, -2.0,-1.5, -1, -0.5,0,0.5])  # 设置x轴刻度
plt.yticks([-1.5, -1, -0.5, 0, 0.5, 1, 1.5])  # 设置y轴刻度
plt.tick_params(labelsize=5)  # 设置刻度字号
plt.plot(x0, y0)  # 第一个data表示选取data为数据集，第二个是函数，data的平方
plt.plot(x, y,linewidth=1)  # 同上
plt.scatter(x1,y1,marker='x',s=20) #插值节点
plt.legend(['Heart Line', 'Cubic Spline','Insert Nude'])  # 打出图例
plt.show()  # 显示图形
