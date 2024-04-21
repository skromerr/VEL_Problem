import matplotlib.pyplot as plt
import numpy as np

x, y, Az, u = [], [], [], []

with open("results.txt") as file:
    for line in file:
        xT, yT, AzT = line.split()
        x.append(float(xT))
        y.append(float(yT))
        Az.append(float(AzT))

with open("kek.txt") as file:
    for line in file:
        u.append(float(line))

x = list(dict.fromkeys(x))
y = list(dict.fromkeys(y))

levels = np.linspace(min(u), max(u), 8)



# plot
fig, ax = plt.subplots()

X, Y = np.meshgrid(x, y)
u = np.reshape(u, (len(y), len(x)))
colorBar = plt.contourf(X, Y, u, levels=100, cmap='jet')
plt.colorbar(colorBar, ax=ax, format='%.0e')
plt.show()
