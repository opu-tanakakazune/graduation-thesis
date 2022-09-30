import cmath

old_phi = []
phi = []
v = []

dx = 0.05
expected_x = -0.5
expected_p = 100.0
h = 0.01
h2 = 0.0000001

# 初期値の設定
for i in range(400):
    x = -2.0 + h*i
    a = -(x-expected_x)*(x-expected_x)/4.0/dx/dx
    b = (expected_p * x)
    phi_value = ((2.0*cmath.pi*dx*dx)**(-0.25) * cmath.exp(a+b*1j))
    v_value = 0
    if x<0.5 and x>0:
        v_value = 8000
    old_phi += [phi_value]
    v += [v_value]
s0 = 0
#  t=0のdataファイル作成
with open('schrodinger1.dat', 'w') as f1:
    for i in range(400):
        x = -2.0 + h*i
        expected_value = (abs(old_phi[i]))**2
        s0 += expected_value*h
        print(x, expected_value, v[i]/4000, file=f1)

def calculate_phi(phi_p, phi_0, phi_m, v):
    a = (h2/(h*h)) * (phi_p -2*phi_0 + phi_m)
    b = h2*v*phi_0
    answer = phi_0 + (a*1j) - (b*1j)
    return answer

for i in range(20000):
    phi =[]
    v = []
    for j in range(400):
        x = -2.0 + h*j
        v_value = 0
        if x<0.5 and x>0:
            v_value = 8000
        if j == 0:
            phi += [calculate_phi(old_phi[j+1], old_phi[j], 0, v_value)]
            v += [v_value]
        elif j == 399:
            phi += [calculate_phi(0, old_phi[j], old_phi[j-1], v_value)]
            v += [v_value]
        else:
            phi += [calculate_phi(old_phi[j+1], old_phi[j], old_phi[j-1], v_value)]
            v += [v_value]
    old_phi = phi
#  t=20000のdataファイル作成
s1 = 0
with open('schrodinger2.dat', 'w') as f2:
    for i in range(400):
        x = -2.0 + h*i
        expected_value = (abs(phi[i]))**2
        s1 += expected_value*h
        print(x, expected_value, v[i]/4000, file=f2)

for i in range(20000):
    phi =[]
    v = []
    for j in range(400):
        x = -2.0 + h*j
        v_value = 0
        if x<0.5 and x>0:
            v_value = 8000
        if j == 0:
            phi += [calculate_phi(old_phi[j+1], old_phi[j], 0, v_value)]
            v += [v_value]
        elif j == 399:
            phi += [calculate_phi(0, old_phi[j], old_phi[j-1], v_value)]
            v += [v_value]
        else:
            phi += [calculate_phi(old_phi[j+1], old_phi[j], old_phi[j-1], v_value)]
            v += [v_value]
    old_phi = phi
s2 = 0
#  t=40000のdataファイル作成
with open('schrodinger3.dat', 'w') as f3:
    for i in range(400):
        x = -2.0 + h*i
        expected_value = (abs(phi[i]))**2
        s2 += expected_value*h
        print(x, expected_value, v[i]/4000, file=f3)

for i in range(20000):
    phi =[]
    v = []
    for j in range(400):
        x = -2.0 + h*j
        v_value = 0
        if x<0.5 and x>0:
            v_value = 8000
        if j == 0:
            phi += [calculate_phi(old_phi[j+1], old_phi[j], 0, v_value)]
            v += [v_value]
        elif j == 399:
            phi += [calculate_phi(0, old_phi[j], old_phi[j-1], v_value)]
            v += [v_value]
        else:
            phi += [calculate_phi(old_phi[j+1], old_phi[j], old_phi[j-1], v_value)]
            v += [v_value]
    old_phi = phi
s3 = 0
#  t=60000のdataファイル作成
with open('schrodinger4.dat', 'w') as f4:
    for i in range(400):
        x = -2.0 + h*i
        expected_value = (abs(phi[i]))**2
        s3 += expected_value*h
        print(x, expected_value, v[i]/4000, file=f4)

for i in range(20000):
    phi =[]
    v = []
    for j in range(400):
        x = -2.0 + h*j
        v_value = 0
        if x<0.5 and x>0:
            v_value = 8000
        if j == 0:
            phi += [calculate_phi(old_phi[j+1], old_phi[j], 0, v_value)]
            v += [v_value]
        elif j == 399:
            phi += [calculate_phi(0, old_phi[j], old_phi[j-1], v_value)]
            v += [v_value]
        else:
            phi += [calculate_phi(old_phi[j+1], old_phi[j], old_phi[j-1], v_value)]
            v += [v_value]
    old_phi = phi
s4 = 0
#  t=80000のdataファイル作成
with open('schrodinger5.dat', 'w') as f5:
    for i in range(400):
        x = -2.0 + h*i
        expected_value = (abs(phi[i]))**2
        s4 += expected_value*h
        print(x, expected_value, v[i]/4000, file=f5)

for i in range(20000):
    phi =[]
    v = []
    for j in range(400):
        x = -2.0 + h*j
        v_value = 0
        if x<0.5 and x>0:
            v_value = 8000
        if j == 0:
            phi += [calculate_phi(old_phi[j+1], old_phi[j], 0, v_value)]
            v += [v_value]
        elif j == 399:
            phi += [calculate_phi(0, old_phi[j], old_phi[j-1], v_value)]
            v += [v_value]
        else:
            phi += [calculate_phi(old_phi[j+1], old_phi[j], old_phi[j-1], v_value)]
            v += [v_value]
    old_phi = phi
s5 = 0
#  t=100000のdataファイル作成
with open('schrodinger6.dat', 'w') as f6:
    for i in range(400):
        x = -2.0 + h*i
        expected_value = (abs(phi[i]))**2
        s5 += expected_value*h
        print(x, expected_value, v[i]/4000, file=f6)

for i in range(20000):
    phi =[]
    v = []
    for j in range(400):
        x = -2.0 + h*j
        v_value = 0
        if x<0.5 and x>0:
            v_value = 8000
        if j == 0:
            phi += [calculate_phi(old_phi[j+1], old_phi[j], 0, v_value)]
            v += [v_value]
        elif j == 399:
            phi += [calculate_phi(0, old_phi[j], old_phi[j-1], v_value)]
            v += [v_value]
        else:
            phi += [calculate_phi(old_phi[j+1], old_phi[j], old_phi[j-1], v_value)]
            v += [v_value]
    old_phi = phi
s6 = 0
#  t=120000のdataファイル作成
with open('schrodinger7.dat', 'w') as f7:
    for i in range(400):
        x = -2.0 + h*i
        expected_value = (abs(phi[i]))**2
        s6 += expected_value*h
        print(x, expected_value, v[i]/4000, file=f7)

with open('schrodingerS.dat', 'w') as f8:
    print(0, s0, file=f8)
    print(1, s1, file=f8)
    print(2, s2, file=f8)
    print(3, s3, file=f8)
    print(4, s4, file=f8)
    print(5, s5, file=f8)
    print(6, s6, file=f8)
