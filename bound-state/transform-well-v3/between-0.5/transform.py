import cmath
import math

old_phi = []
phi = []
v = []
s = []
alphas = []
hamis = []

l = 1
h = 0.01
h2 = 0.000001
n = 2000
xmin = 0.01 * (n / 2)
xay = 1.278839111328125
eta = 4.2544867725581845
constant_v = cmath.pi * cmath.pi
alpha = h2/(2*h*h)
beta = h2/2
for i in range(n-1):
    alphas += [alpha]

# 初期値の設定
for i in range(n):
    x = -xmin + h*i
    k = xay / l
    ka = eta / l
    C = 1/((l + (1/ka))**(1/2))
   
    if -l <= x and x <= l:
        phi_value = C * cmath.cos(k * x)
    else:
        phi_value = C * cmath.cos(k * l) * cmath.exp(-1 * (ka * (abs(x) - l)))
    v_value = -constant_v
    # if -l - 0.5 <= x and x < -l:
    #     v_value = 0
    # if l < x and x <= l + 0.5:
    #     v_value = 0
    if x < -l:
        v_value = 0
    if l < x:
        v_value = 0
    old_phi += [phi_value]
    v += [v_value]
#  t=0のdataファイル作成
with open('bound-state/transform-well-v3/between-0.5/schrodingerV_W.dat', 'w') as f1:
    hami = 0
    for i in range(n):
        x = -xmin + h*i
        expected_value = (abs(old_phi[i]))**2
        if -l - 0.5 <= x and x <= l + 0.5:
            hami += expected_value*h
        print(x, expected_value, v[i]/20 + constant_v/20, math.log(expected_value), file=f1)
    hamis += [hami]

#3重対格行列計算用メソッド
#代入方法は参考画像参照
def TDMA(a, b, c, d):
    e = [c[0]/b[0]]
    f = [d[0]/b[0]]
    for i in range(1, len(d)-1):
        e.append(c[i]/(b[i]-a[i-1]*e[i-1]))
        f.append((d[i]-a[i-1]*f[i-1])/(b[i]-a[0]*e[0]))
    u = [(d[len(d)-1]-a[len(d)-2]*f[len(d)-2]) /
         (b[len(d)-1]-a[len(d)-2]*e[len(d)-2])]
    for i in range(len(d)-2, -1, -1):
        u.append(f[i]-e[i]*u[len(d)-2-i])
    u.reverse()
    return u

def calculate_b(alpha, beta, v):
    b = 1j -2*alpha -beta*v
    return b

def calculate_d(alpha, beta, v, phi1, phi2, phi3):
    d = -alpha*phi3 +(1j +2*alpha +beta*v)*phi2 -alpha*phi1
    return d

for m in range(2,2):
    file_name = 'bound-state/transform-well-v3/between-0.5/schrodinger' + str(m) + '.dat'
    f_number = 'f' + str(m)
    for i in range(20000):
        phi =[]
        v = []
        b = []
        d = []
        for j in range(n):
            x = -xmin + h*j
            v_value = -constant_v
            if -l -0.5 <= x and x < -l:
                v_value = 0
            if l < x and x <= l + 0.5:
                v_value = 0
            b += [calculate_b(alpha, beta, v_value)]
            if j != 0 and j != n-1:
                d += [calculate_d(alpha, beta, v_value, old_phi[j-1], old_phi[j], old_phi[j+1])]
            if j == 0:
                d += [-alpha*old_phi[j+1] +(1j +2*alpha +beta*v_value)*old_phi[j]]
            if j == n-1:
                d += [calculate_d(alpha, beta, v_value, old_phi[j-1], old_phi[j], 0)]
            v += [v_value]
        phi = TDMA(alphas, b, alphas, d)
        old_phi = phi
    with open(file_name, 'w') as f_number:
        hami = 0
        for i2 in range(n):
            x = -xmin + h*i2
            expected_value = (abs(phi[i2]))**2
            if -l - 0.5 <= x and x <= l + 0.5:
                hami += expected_value*h
            print(x, expected_value, v[i2]/20 + constant_v/20, math.log(expected_value), file=f_number)
        hamis += [hami]
    print(m)

count = 0
with open('bound-state/transform-well-v3/between-0.5/schrodingerHami.dat', 'w') as f100:
    for item in hamis:
        print(count+1, math.log(hamis[count]), file=f100)
        count = count +1
