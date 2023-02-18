import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
import numpy as np
import argparse
import ast

parser = argparse.ArgumentParser(description='Surface plot.')
parser.add_argument("-x", "--list_x", type=str, required=True)
parser.add_argument("-y", "--list_y", type=str, required=True)
parser.add_argument("-z", "--list_z", type=str, required=True)

args = parser.parse_args()

x = ast.literal_eval(args.list_x)
y = ast.literal_eval(args.list_y)
z = ast.literal_eval(args.list_z)

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
# ax.set_zlim([-13.5, -14.7])

plt.show()