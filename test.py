h = 0.01
h2 = 0.0000001

def calculate_phi(phi_p, phi_0, phi_m, v):
    a = (h2/(h*h)) * (phi_p -2*phi_0 + phi_m)
    b = h2*v*phi_0
    answer = phi_0 + (a*1j) - (b*1j)
    return answer
phi = calculate_phi(0.32108751086095555+0j, 0.31807154798425885+0j, 0.3150839139330875+0j, 0)
print(phi)
print((abs(phi))**2)

0.31807154798425885+0j
0.31807154798425885+2.832882552533666e-08j