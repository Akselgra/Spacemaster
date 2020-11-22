import numpy as np
import matplotlib.pyplot as plt

class SynthDat():
    """
    Class for creating synthetic 2d data that mimics SWARM data.
    """

    def __init__(self, fs = 2, v = 7615, n = 1000):
        """
        Constructor.
        ----------------------------------------
        Takes:
            fs - sampling frequency [Hz]
            v - velocity of satellite [m/s]
            n - size of one grid dimension
        """
        self.fs = fs
        self.v = v
        self.dx = v/fs #grid point distance
        self.n = n


    def gauss_curve2d(self, amplitude, x, y, meanx, meany, stdx, stdy):
        """
        2D Gaussian bell function
        """
        first = (x - meanx)**2/(2*stdx**2)
        second = (y - meany)**2/(2*stdy**2)
        return(amplitude*np.exp(-(first + second)))

    def flatfunc(self,A, x, y, width, stepwidth , xpos, ypos):
        """
        Defines a flat cylinder with continous borders
        """
        #A = self.amplitude
        B = self.background

        width_circ = width + stepwidth
        incline = (A-B)/(stepwidth)
        r = np.sqrt((x - xpos)**2 + (y - ypos)**2)

        if r <= width:
            return A

        elif width < r < width_circ:
            return(-(r- width)*incline + A)

        else:
            return B


    def flatfunc_vec(self, A, satpos, width, stepwidth, bobpos):
        """
        Generalized flatfunc for any dimension
        -----------------------------------------------
        takes:
            A - scalar max amplitude
            satpos - array with [x, y, ...] positions of satellite
            width - width of bubble
            stepwidth - width of incline
            bobpos - array with [x, y, ...] positions of bubble

        returns:
            Electron density at satellite position
        """

        width_circ = width + stepwidth
        incline = (A)/stepwidth
        r = np.sqrt(np.sum((satpos - bobpos)**2))

        if r <= width:
            return A

        elif width < r < width_circ:
            return((width - r)*incline + A)

        else:
            return 0

    def add_bubble(self, A, width, stepwidth, bobpos):
        """
        Adds bubble to data from flatfunc_vec method
        """
        if self.data == True:
            self.data = np.zeros((n, n))
        xs = np.arange(self.n)*obj.dx
        ys = np.arange(self.n)*obj.dx
        Xs, Ys = np.meshgrid(xs, ys)
        for i in range(len(self.data)):
            for j in range(len(self.data[i])):
                self.data[i,j] += obj.flatfunc_vec(A, np.array([Xs[i, j], Ys[i, j]]), width, stepwidth, bobpos)


    def density(self, r, C = 6.022e23):
        """
        The density of a bubble with a constant electron count.
        standard of C is set to Avogadros number
        """
        area = 4/3*np.pi*r**3
        return(C/area)


    def multibob(self, satpos, bobpos, width, stepwidth, A):
        """
        Function for summing over multiple flatfunc bubbles.
        takes:
            satpos - array with [x, y, ...] positions of satellite
            bobpos - list of bubble position arrays
            width - list of widths
            stepwidth - list of stepwidths
            A - list of amplitudes

        returns
            ne - electron density at position
        """
        test = len(bobpos) == len(width) == len(stepwidth) == len(A)
        msg = "Variables in multibob were not of equal length"
        assert test, msg

        ne = 0
        for i in range(len(bobpos)):
            ne = ne + self.flatfunc_vec(A[i], satpos, width[i],\
                                    stepwidth[i], bobpos[i])
        return(ne)

    def satrun(self,t0, t1, satpos, satdir, bobpos, bobvel, width, C, diffus):
        """
        Generates synthetic data for one satellite.
        takes:
            t0 - initial sampling time [s]
            t1 - final sampling time [s]
            satpos - array with initial satellite positions [m]
            satdir - direction of satellite velocity
            bobpos - list of initial bubble positions [m]
            bobvel - list of bubble velocities [m/s]
            width - list of initial bubble widths [m]
            C - list of total electron content in bubble
            diffus - time derivative of width [m/s]
        returns:
        nes - array with electron density values
        """

        n = int((t1 - t0)*self.fs) #nr of samples
        times = np.linspace(t0, t1, n)
        nes = np.zeros(n)
        stepwidth = np.zeros_like(width)
        A = np.zeros_like(width)
        dt = times[1] - times[0]

        for i in range(len(bobpos)):
            bobpos[i] = bobpos[i] + t0*bobvel[i]
            width[i] = width[i] + t0*diffus[i]

        for i in range(n):
            satpos = satpos + satdir*self.v*dt#update satellite position

            for j in range(len(bobpos)):
                bobpos[j] = bobpos[j] + bobvel[j]*dt #update bubble positions
                width[j] = width[j] + diffus[j]*dt #update bubble width
                stepwidth[j] = width[j]/10
                A[j] = self.density(width[j], C[j])
            nes[i] = self.multibob(satpos, bobpos, width, stepwidth, A)

        return(nes)


    def satrun2(self,t0, t1, satpos, satdir, bobpos, bobvel, width, A, growth):
        """
        Generates synthetic data for one satellite.
        takes:
            t0 - initial sampling time [s]
            t1 - final sampling time [s]
            satpos - array with initial satellite positions [m]
            satdir - direction of satellite velocity
            bobpos - list of initial bubble positions [m]
            bobvel - list of bubble velocities [m/s]
            width - list of initial bubble widths [m]
            A - list of electron density in bubble
            growth - time derivative of width [m/s]
        returns:
        nes - array with electron density values
        """

        n = int((t1 - t0)*self.fs) #nr of samples
        times = np.linspace(t0, t1, n)
        nes = np.zeros(n)
        stepwidth = np.zeros_like(width)
        dt = times[1] - times[0]

        for i in range(n):
            satpos = satpos + satdir*self.v*dt #update satellite position

            for j in range(len(bobpos)):
                bobpos[j] = bobpos[j] + bobvel[j]*dt #update bubble positions
                width[j] = width[j] + growth[j]*dt #update bubble width
                stepwidth[j] = width[j]/10
            nes[i] = self.multibob(satpos, bobpos, width, stepwidth, A)

        return(nes)


if __name__ == "__main__":
    t0 = 0
    t1 = 400
    fs = 2
    n = int((t1 - t0)*fs)
    v = 7615
    obj = SynthDat(fs = fs, v = v, n = n)


    width = np.array([17,18,19,20])*v

    c = 6.022e23
    C = np.zeros_like(width)+c

    diffus = np.array([0.01,0.01,0.01,0.01])*v

    bobpos1 = np.array([50, 0])*v
    bobpos2 = np.array([100, 0])*v
    bobpos3 = np.array([150, 0])*v
    bobpos4 = np.array([250, 0])*v
    bobpos = np.array([bobpos1, bobpos2, bobpos3, bobpos4])

    bobvel1 = np.array([0, 0])*v
    bobvel2 = np.array([0, 0])*v
    bobvel3 = np.array([0, 0])*v
    bobvel4 = np.array([0, 0])*v
    bobvel = np.array([bobvel1, bobvel2, bobvel3, bobvel4])

    satpos = np.array([0, 0])
    satdir = np.array([1, 0])

    width = np.array([50])*v
    C = np.array([c])
    diffus = np.array([0.01])*v
    bobpos = np.array([np.array([150, 0])])*v
    bobvel = np.array([np.array([-0.05, 0])])*v

    satpossies = np.linspace(0, t1*v, n)
    xs = np.linspace(-500*v, 500*v, 500)
    ys = np.linspace(-500*v, 500*v, 500)
    Xs, Ys = np.meshgrid(xs, ys)

    # data = np.zeros_like(Xs)
    # for i in range(len(data)):
    #     for j in range(len(data[i])):
    #         possies = np.array([Xs[i, j], Ys[i, j]])
    #         data[i, j] = obj.multibob(possies, bobpos, width, width/5, obj.density(width))
    #
    # plt.figure()
    # plt.pcolormesh(Xs, Ys, data, cmap = "magma")
    # plt.colorbar()
    # plt.plot(satpossies, np.zeros_like(satpossies), "r.")
    # plt.legend(["Satellite path"])
    # plt.show()

    nes = obj.satrun2(t0, t1, satpos, satdir, bobpos, bobvel, width, C, diffus)

    t0 += 50
    t1 += 50
    n = int((t1 - t0)*fs)
    nes2 = obj.satrun2(t0, t1, satpos, satdir, bobpos, bobvel, width, C, diffus)

    times = np.linspace(t0, t1, n)
    xs = np.linspace(0, v*(t1 - t0), n)
    plt.plot(xs/v, nes)
    plt.plot(xs/v, nes2)
    plt.xlabel("x/v (basically time)")
    plt.ylabel("Electron density")
    plt.legend(["t0 = 0", "t0 = 50"])
    plt.show()
