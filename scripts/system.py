from sympy import solve
from sympy import symbols, var, simplify, latex, collect

d_1, d_3, f_1, f_3 = var('d_1^0, d_3^0, f_1^0, f_3^0')
b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16 = var('b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16')

system = [
    b1*d_1 + b2*d_3 + b3*f_1 + b4*f_3 - 1,
    b5*d_1 + b6*d_3 + b7*f_1 + b8*f_3,
    b9*d_1 + b10*d_3 + b11*f_1 + b12*f_3,
    b14*d_3 + b16*f_3
]
find = [d_1, d_3, f_1, f_3]
res = solve(system, find)