import numpy as np
import matplotlib.pyplot as plt

class SynthDat():
    """
    Class for creating synthetic 2d data that mimics SWARM data.
    """

    def __init__(self, fs = 2, v = 7615, n = 1000, m = 1, amplitude = 2.5e5,\
    background = 0):
        """
        Constructor.
        ----------------------------------------
        Takes:
            fs - sampling frequency [Hz]
            v - velocity of satellite [m/s]
            n - size of one grid dimension
            m - number of tops in data
        """
        self.fs = fs
        self.v = v
        self.dx = v/fs #grid point distance
        self.n = n
        self.m = m
        self.amplitude = amplitude
        self.background = background


    def gauss_curve2d(self, amplitude, x, y, meanx, meany, stdx, stdy):
        """
        2D Gaussian bell function
        """
        first = (x - meanx)**2/(2*stdx**2)
        second = (y - meany)**2/(2*stdy**2)
        return(amplitude*np.exp(-(first + second)))

    def gaussgrid(self, std = False):
        """
        Initializes a grid
        """
        A = self.amplitude

        xs = np.arange(n)*self.dx
        ys = np.arange(n)*self.dx
        Xs, Ys = np.meshgrid(xs, ys)

        data = np.zeros_like(Xs)

        partsize = int(self.n/self.m) #grid points in each top
        meanx = partsize/2*self.dx # position of tops in  x direction
        meany = self.n/2*self.dx
        if std == False:
            std = meanx/2

        for m in range(self.m): #filling data from gauss_curve2d
            curr_mean = meanx + partsize*m*self.dx
            ind0 = partsize*m
            for i in range(self.n):
                for j in range(partsize):
                    x_ind = i
                    y_ind = ind0 + j
                    x = Xs[x_ind, y_ind]
                    y = Ys[x_ind, y_ind]
                    data[x_ind, y_ind] = self.gauss_curve2d(A, x, y, curr_mean,\
                                                            meany, std, std)
        self.Xs = Xs
        self.Ys = Ys
        self.data = data


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
        
        B = self.background
        width_circ = width + stepwidth
        incline = (A-B)/stepwidth
        r = np.sqrt(np.sum((satpos - bobpos)**2))
        
        if r <= width:
            return A

        elif width < r < width_circ:
            return(-(r- width)*incline + A)

        else:
            return B

    def flatgrid(self, width = True, stepwidth = True):
        """
        Initializes a grid with flat bubbles
        """
        A = self.amplitude  
        xs = np.arange(n)*self.dx
        ys = np.arange(n)*self.dx
        Xs, Ys = np.meshgrid(xs, ys)

        data = np.zeros_like(Xs)

        partsize = int(self.n/self.m) #grid points in each top
        meanx = partsize/2*self.dx # position of tops in  x direction
        meany = self.n/2*self.dx
        if width:
            width = meanx/2
        if stepwidth:
            stepwidth = width/4

        for m in range(self.m): #filling data from flatfunc
            curr_mean = meanx + partsize*m*self.dx
            ind0 = partsize*m
            for i in range(self.n):
                for j in range(partsize):
                    x_ind = i
                    y_ind = ind0 + j
                    x = Xs[x_ind, y_ind]
                    y = Ys[x_ind, y_ind]
                    data[x_ind, y_ind] = self.flatfunc(A, x, y, width, stepwidth,\
                                                       curr_mean, meany)
        self.Xs = Xs
        self.Ys = Ys
        self.data = data



    def rand_flat_grid(self):
        """
        Initializes a grid with flat bubbles
        with random amplitude, width and stepwidth
        """

        A_base = self.amplitude  
        xs = np.arange(n)*self.dx
        ys = np.arange(n)*self.dx
        Xs, Ys = np.meshgrid(xs, ys)

        data = np.zeros_like(Xs)

        partsize = int(self.n/self.m) #grid points in each top
        meanx = partsize/2*self.dx # position of tops in  x direction
        meany = self.n/2*self.dx
    
        width_base = meanx/2
        

        for m in range(self.m): #filling data from flatfunc
            curr_mean = meanx + partsize*m*self.dx
            ind0 = partsize*m
            
            width = width_base * np.random.uniform(0.2, 1)#randomize variables
            stepwidth = width/4
            A = A_base*np.random.uniform(0.2, 1)
            
            for i in range(self.n):
                for j in range(partsize):
                    x_ind = i
                    y_ind = ind0 + j
                    x = Xs[x_ind, y_ind]
                    y = Ys[x_ind, y_ind]
                    satpos = np.array([x, y])
                    bubpos = np.array([curr_mean, meany])
                    data[x_ind, y_ind] = self.flatfunc_vec(A, satpos, width, stepwidth,\
                                                       bubpos)
        self.Xs = Xs
        self.Ys = Ys
        self.data = data




if __name__ == "__main__":
    fs = 2
    v = 7615
    n = 500
    m = 1
    amplitude = 225000
    obj = SynthDat(fs = fs, v = v, n = n, m = m, amplitude = amplitude)
    obj.rand_flat_grid()
    data = obj.data
    Xs = obj.Xs
    Ys = obj.Ys

    plt.pcolormesh(Xs, Ys, data, cmap = "magma")
    plt.colorbar()
    plt.show()

    plt.plot(Xs[int(n/2), :], data[int(n/2), :])
    plt.show()
