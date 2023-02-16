import matplotlib.pyplot as plt
from matplotlib import cm
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

Z = np.array(z)

fig1 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')

ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

plt.legend()
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("function")

plt.show()