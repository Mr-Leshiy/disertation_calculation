import matplotlib.pyplot as plt
from matplotlib import cm
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

X, Y = np.meshgrid(x, y)

Z = np.array(z)

fig1 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')

ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

plt.legend()
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("sigma_y(x,y)")

plt.show()