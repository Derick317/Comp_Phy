#import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3d

with open("Lorenz_attractor_data.txt") as f:
    data = f.readlines()
t = []
y1 = []
y2 = []
y3 = []
for i in range(len(data)):
    item = data[i].split()
    t.append(float(item[0]))
    y1.append(float(item[1]))
    y2.append(float(item[2]))
    y3.append(float(item[3]))

fig = plt.figure()
ax = p3d.Axes3D(fig)
ax.plot(y1, y2, y3, 'green')
ax.plot(y1[0], y2[0],y3[0], 'r*')
plt.show()
