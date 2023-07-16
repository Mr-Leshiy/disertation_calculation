import matplotlib.pyplot as plt
from matplotlib import cm
import argparse
import ast

parser = argparse.ArgumentParser(description='Function plot.')
parser.add_argument("-p", "--path", type=str, required=True)

file = open(parser.parse_args().path)
data = file.read().split("|")

omega = ast.literal_eval(data[0])
z = ast.literal_eval(data[1])

plt.plot(omega, z,)

plt.legend()
plt.xlabel("omega")
plt.ylabel("Function")

plt.show()