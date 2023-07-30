import matplotlib.pyplot as plt
from matplotlib import cm
import argparse
import ast

parser = argparse.ArgumentParser(description='Function plot.')
parser.add_argument("-p", "--path", type=str, required=True)

file = open(parser.parse_args().path)
data = file.read().split("|")

ox_name = data[0]
oy_name = data[1]

i = 2
legend = []
while i < len(data):
    x = ast.literal_eval(data[i])
    y = ast.literal_eval(data[i + 1])
    legend.append(data[i + 2])
    plt.plot(x, y,)
    i += 3

plt.legend(legend)
plt.xlabel(ox_name)
plt.ylabel(oy_name)

plt.show()