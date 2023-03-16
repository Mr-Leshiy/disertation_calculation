import sympy as sp
import numpy as np

g, b, m, l, a, e = sp.symbols("g, b, m, l, a, e")

x1 = 2*b*m + 2*a + 2*a*m
x2 = 2*b*m - 2*a
x3 = -2*b*m + 2*a + 2*a*m
x4 = 2*b*m + 2*a
x5 = -2*g*b*m - 2*g*a*m + 2*a*l
x6 = -2*g*b*m + 4*g*a + 2*a*l
x7 = -2*g*b*m + 2*g*a*m - 2*a*l
x8 = 2*g*b*m + 4*g*a + 2*a*l
x9 = 2 + 2*m
x10 = -2
x11 = 2 + m
x12 = -2 - m

y1 = x2 + e*x4
y2 = e*x3 - x1
y3 = x6 + e*x8
y4 = e*x7 - x5

y_tmp = sp.simplify(y2*y3 - y1*y4)

d_1_0 = -(1 + m)*4 * sp.simplify(y4*x10*x12 + y4*x10*x11 + y3*x9*x11) / sp.simplify(x9*x11*y_tmp)
d_1_0 = sp.simplify(d_1_0)
d_3_0 = (1 + m)*4 * sp.simplify(x12 * y4) / sp.simplify(x11 * y_tmp)
d_3_0 = sp.simplify(d_3_0)
d_2_1 = -(1 + m)*4 * sp.simplify(sp.simplify(x1*x10 - x2*x9)*y3 - sp.simplify(x5*x10 - x6*x9)*y1) / sp.simplify(x9*x11*y_tmp) - (1 + m)*4*e*sp.simplify(x10/x9/x11)
d_2_1 = sp.simplify(d_2_1)
d_4_1 = (1 + m)*4 * sp.simplify(sp.simplify(x5*x10 - x6*x9)*y2 - sp.simplify(x1*x10 - x2*x9)*y4) / sp.simplify(x9*x11*y_tmp) + (1 + m)*4*e/x11
d_4_1 = sp.simplify(d_4_1)

# (d_1_0 + d_3_0)
sum1 = sp.simplify(d_1_0 + d_3_0)
# (d_2_1 + d_4_1)
sum2 = sp.simplify(d_2_1 + d_4_1)
# (d_2_1 + d_4_1)*(d_1_0 + d_3_0)
res = sp.simplify(sum1 * sum2)

print("(d_2_1 + d_4_1)*(d_1_0 + d_3_0): [ {} ]".format(res))


