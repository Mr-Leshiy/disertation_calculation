import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
import numpy as np
import argparse
import ast

parser = argparse.ArgumentParser(description='Surface plot.')
parser.add_argument("-p", "--path", type=str, required=True)

file = open(parser.parse_args().path)
data = file.read().split("|")

x = ast.literal_eval(data[0])
y = ast.literal_eval(data[1])
z = ast.literal_eval(data[2])

z_min = min(min(values) for values in [min(values) for values in z])
z_max = max(max(values) for values in [max(values) for values in z])

X, Y = np.meshgrid(x, y)

Z = [np.array(values) for values in z]

fig1 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')

def update_plot(frame_nuber, Z, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(X, Y, Z[frame_nuber], cmap=cm.coolwarm, linewidth=0, antialiased=False)

plot = [ax.plot_surface(X, Y, Z[0], cmap=cm.coolwarm, linewidth=0, antialiased=False)]

ani =  animation.FuncAnimation(fig1, update_plot, len(z), fargs=(Z, plot), interval=50, blit=False)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("function")
ax.set_zlim([z_min, z_max])

plt.show()