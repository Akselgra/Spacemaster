import numpy as np
import matplotlib.pyplot as plt


cv0 = 650
cv1 = 1200
ck = 100000
m = 2.25
K = 273.15


E1 = cv0*m*(5)
E2 = E1 + m*ck
E3 = cv1*(65)*m + E2

Es = np.array([E1, E2, E3])

def T(E):
    cv0 = 650
    cv1 = 1200
    ck = 100000
    m = 2.25
    K = 273.15
    
    
    E1 = cv0*m*(5)
    E2 = E1 + m*ck
    E3 = cv1*(65)*m + E2
    
    Es = np.array([E1, E2, E3])
    
    if E < E1:
        return 273.15 - 30 + E/(cv0*m)
    if E >= E1 and E < E2:
        return 273.15 - 25
    if E > E2:
        return(273.15 - 25 + (E - E2)/(cv1*m))
    
Es = np.linspace(0, E3, 1000)
Ts = np.zeros_like(Es)
for i in range(len(Ts)):
    Ts[i] = T(Es[i])
    
Ts = Ts - 273.15
Es = Es/1000
plt.plot(Es, Ts)
plt.xticks([E1/1000, E2/1000, E3/1000])
plt.yticks([-30, -25, 40])
plt.xlabel("Energy [kJ]")
plt.ylabel("Temperature [$^\circ$C]")
plt.grid("on")
plt.show()
