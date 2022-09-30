import math

old_a = []
old_b = []
a = []
b = []
v = []

dx = 0.05
expected_x = -0.5
expected_p = 100.0
h = 0.01
h2 = 0.0000001

# 初期値の設定
for i in range(400):
    x = -2.0 + h*i
    a_value = ((2.0*math.pi*dx*dx)**(-0.25)) * math.exp(-(x-expected_x)*(x-expected_x)/4.0/dx/dx) * math.cos(expected_p*x)
    b_value = ((2.0*math.pi*dx*dx)**(-0.25)) * math.exp(-(x-expected_x)*(x-expected_x)/4.0/dx/dx) * math.sin(expected_p*x)
    v_value = 0
    if x<0.5 and x>0:
        v_value = 8000
    old_a += [a_value]
    old_b += [b_value]
    v += [v_value]

#t=0のdataファイル作成
with open('schrodinger1.dat', 'w') as f1:
    for i in range(400):
        x = -2.0 + h*i
        expected_value = (old_a[i]**2.0) + (old_b[i]**2.0)
        print(x, expected_value, v[i], file=f1)

def calculate_a(a, b1, b, b2, v):
    answer = a - (h2/(h*h))*(b1 -(2*b) + b2) + h2*v*b
    return answer

def calculate_b(b, a1, a, a2, v):
    answer = b + (h2/(h*h))*(a1 -(2*a) + a2) - h2*v*a
    return answer

for i in range(20000):
    a = []
    b = []
    v = []
    for j in range(400):
        x = -2.0 + h*j
        v_value = 0
        if x<0.5 and x>0:
            v_value = 8000
        if j == 0:
            a += [calculate_a(old_a[j], old_b[j+1], old_b[j], 0, v_value)]
            b += [calculate_b(old_b[j], old_a[j+1], old_b[j], 0, v_value)]
            v += [v_value]
        elif j == 399:
            a += [calculate_a(old_a[j], 0, old_b[j], old_b[j-1], v_value)]
            b += [calculate_b(old_b[j], 0, old_b[j], old_b[j-1], v_value)]
            v += [v_value]
        else:
            a += [calculate_a(old_a[j], old_b[j+1], old_b[j], old_b[j-1], v_value)]
            b += [calculate_b(old_b[j], old_a[j+1], old_b[j], old_b[j-1], v_value)]
            v += [v_value]
    old_a = a
    old_b = b
#t=20000のdataファイル作成
with open('schrodinger2.dat', 'w') as f2:
    for i in range(400):
        x = -2.0 + h*i
        expected_value = a[i]**2.0 + b[i]**2.0
        print(x, expected_value, v[i], file=f2)

for i in range(20000):
    a = []
    b = []
    v = []
    for j in range(400):
        x = -2.0 + h*j
        v_value = 0
        if x<0.5 and x>0:
            v_value = 8000
        if j == 0:
            a += [calculate_a(old_a[j], old_b[j+1], old_b[j], 0, v_value)]
            b += [calculate_b(old_b[j], old_a[j+1], old_b[j], 0, v_value)]
            v += [v_value]
        elif j == 399:
            a += [calculate_a(old_a[j], 0, old_b[j], old_b[j-1], v_value)]
            b += [calculate_b(old_b[j], 0, old_b[j], old_b[j-1], v_value)]
            v += [v_value]
        else:
            a += [calculate_a(old_a[j], old_b[j+1], old_b[j], old_b[j-1], v_value)]
            b += [calculate_b(old_b[j], old_a[j+1], old_b[j], old_b[j-1], v_value)]
            v += [v_value]
    old_a = a
    old_b = b
#t=40000のdataファイル作成
with open('schrodinger3.dat', 'w') as f3:
    for i in range(400):
        x = -2.0 + h*i
        expected_value = a[i]**2.0 + b[i]**2.0
        print(x, expected_value, v[i], file=f3)

print(9**(-1/4))