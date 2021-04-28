import numpy as np
import matplotlib.pyplot as plt
import SWARMprocess
pro = SWARMprocess.SWARMprocess()


fig_width_pt = 360  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inches
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*golden_mean       # height in inches
fig_height = fig_width
fig_size = [fig_width,fig_height]

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Roman",
#   "font.family": "sans-serif",
    # "font.sans-serif": ["Helvetica"],
    'font.size' : 10,
    'axes.labelsize' : 10,
    'font.size' : 10,
#    'text.fontsize' : 10,
    'legend.fontsize': 10,
    'xtick.labelsize' : 10,
    'ytick.labelsize' : 10,
    'figure.figsize': fig_size
})
#matplotlib.use('pgf')

Re = 6370000
xs = np.linspace(-3*Re, 3*Re, 1000)
zs = np.copy(xs)
# zs = np.linspace(-5*Re, 5*Re, 1000)


Xs, Zs = np.meshgrid(xs, zs)
Ys = np.zeros_like(Xs)

vecs = pro.dipol(Xs, Ys, Zs)
xvecs = vecs[0]
zvecs = vecs[2]

xvecs[np.where(np.sqrt(Xs**2 + Zs**2) < Re)] = 0
zvecs[np.where(np.sqrt(Xs**2 + Zs**2) < Re)] = 0

thetas = np.linspace(0, 2*np.pi, 1000)
circ_x = Re*np.cos(thetas)
circ_y = Re*np.sin(thetas)

plt.streamplot(Xs/Re, Zs/Re, xvecs, zvecs, color = "k")
plt.plot(circ_x/Re, circ_y/Re, "k")
plt.axis("equal")
plt.xlabel("[Earth radii]")
plt.ylabel("[Earth radii]")
plt.title("Dipole model of the geomagnetic field")
plt.savefig("Figures/dipole_model.png")
plt.show()



