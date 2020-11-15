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
        self.data = np.zeros((n, n))


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
        
        width_circ = width + stepwidth
        incline = (A)/stepwidth
        r = np.sqrt(np.sum((satpos - bobpos)**2))
        
        if r <= width:
            return A

        elif width < r < width_circ:
            return(-(r- width)*incline + A)

        else:
            return 0

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


    def add_bubble(self, A, width, stepwidth, bobpos):
        """
        Adds bubble to data from flatfunc_vec method
        """
        xs = np.arange(n)*obj.dx
        ys = np.arange(n)*obj.dx
        Xs, Ys = np.meshgrid(xs, ys)
        for i in range(len(self.data)):
            for j in range(len(self.data[i])):
                self.data[i,j] += obj.flatfunc_vec(A, np.array([Xs[i, j], Ys[i, j]]), width, stepwidth, bobpos)
                
    
    def density(self, r, C = 6.022e23):
        """
        The density of a bubble with a constant electron count.
        standard of C is set to Avogadros number/10
        """
        return(3/4 * C/np.pi * 1/(r**3))
        


if __name__ == "__main__":
    fs = 2
    v = 7615
    n = 500
    m = 9
    amplitude = 225000
    obj = SynthDat(fs = fs, v = v, n = n, m = m, amplitude = amplitude)
    '''obj.rand_flat_grid()
    data = obj.data
    Xs = obj.Xs
    Ys = obj.Ys'''
    
    xs = np.arange(n)*obj.dx
    ys = np.arange(n)*obj.dx
    Xs, Ys = np.meshgrid(xs, ys)
    
    pos1 = np.array([obj.dx*250, obj.dx*250])
    width1 = 100*obj.dx
    stepwidth1 = width1/10
    A1 = obj.density(width1)
    
    
    pos2 = np.array([obj.dx*69, obj.dx*100])
    width2 = 40*obj.dx
    stepwidth2 = width1/5
    A2 = obj.density(width2)
    
    pos3 = np.array([obj.dx*50, obj.dx*360])
    width3 = 50*obj.dx
    stepwidth3 = width3/5
    A3 = obj.density(width3)
    
    obj.add_bubble(A1, width1, stepwidth1, pos1)
    obj.add_bubble(A2, width2, stepwidth2, pos2)
    obj.add_bubble(A3, width3, stepwidth3, pos3)
    data = obj.data

    plt.pcolormesh(Xs, Ys, data, cmap = "magma")
    plt.colorbar()
    plt.show()

    plt.plot(Xs[int(n/2), :], data[int(n/2), :])
    plt.show()
