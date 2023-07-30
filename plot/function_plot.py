import matplotlib.pyplot as plt
from matplotlib import cm
import argparse
import ast

parser = argparse.ArgumentParser(description='Function plot.')
parser.add_argument("-p", "--path", type=str, required=True)

file = open(parser.parse_args().path)
data = file.read().split("|")
print(data[0])
print(data[1])

ox_name = data[0]
oy_name = data[1]
x = ast.literal_eval(data[2])
y = ast.literal_eval(data[3])

plt.plot(x, y,)

plt.legend()
plt.xlabel(ox_name)
plt.ylabel(oy_name)

plt.show()