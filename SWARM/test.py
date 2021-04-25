import numpy as np
import matplotlib.pyplot as plt
import SWARMprocess
pro = SWARMprocess.SWARMprocess()


def Br(theta, r, B0 = 3.12e-5, Re = 6370000):
    return(-2*B0*(Re/r)**3 * np.cos(theta))

def Bt(theta, r, B0 = 3.12e-5, Re = 6370000):
    return(-B0*(Re/r)**3 * np.sin(theta))

def B(theta, r, B0 = 3.12e-5, Re = 6370000):
    return(B0*(Re/r)**3*np.sqrt(1 + 3*np.cos(theta)**2))


Re = 6370000
theta = 0



xs = np.concatenate((np.linspace(-3*Re, -Re, 500), np.linspace(Re, 3*Re, 500)))
ys = np.copy(xs)

xs = np.linspace(-10*Re, 10*Re)
ys = np.linspace(-10*Re, 10*Re)
Xs, Ys = np.meshgrid(xs, ys)

Rs = np.sqrt(Xs**2 + Ys**2)
thetas = np.arctan2(Xs, Ys)

Bs = B(thetas, Rs)
Brs = Br(thetas, Rs)/Bs
Bts = Bt(thetas, Rs)/Bs

xdirs = Brs*np.cos(Bts)
ydirs = Brs*np.sin(Bts)

xvecs = xdirs*Bs
yvecs = ydirs*Bs

plt.streamplot(Xs, Ys, xvecs, yvecs)
plt.show()



