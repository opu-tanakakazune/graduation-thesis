import math

v = math.pi * math.pi
l = 1
x1 = 0
x2 = 2

#f(x)を計算・定義する関数
def calculateF(x, v, l):
    f = (x**2)*(1+(math.tan(x))**2) - 2*v*(l**2)
    return f

for i in range(16):
    xm = (x1+x2)/2
    value = calculateF(xm, v, l)
    if value < 0:
        x1 = xm
    else:
        x2 = xm        
print(x1)
print(x2)

value = x1 * math.tan(x1)
print(value)
    