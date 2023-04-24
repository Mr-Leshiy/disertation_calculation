import sympy as sp
import numpy as np

g, b, m, l, a, e = sp.symbols("g, b, m, l, a, e")

y1 = e*b*a*m - e + e**(-1)*b*a*m + e**(-1)
y2 = -e**(-1)*b*a*m + e**(-1) + e**(-1)*m - e*b*a*m - e - e*m
y3 = -e*g*b*a*m + e*2*g + e*l - e**(-1)*g*b*a*m + e**(-1)*2*g + e**(-1)*l
y4 = -e**(-1)*g*b*a*m + e**(-1)*g*m - e**(-1)*l + e*g*b*a*m + e*g*m - e*l

a = sp.expand(y2*y3)
print("a: {}".format(a))
