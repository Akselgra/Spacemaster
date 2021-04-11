import numpy as np
import matplotlib.pyplot as plt
import SWARMprocess
pro = SWARMprocess.SWARMprocess()

def gradient(x, start, stop, val = 100):
    """
    input:
        x - current position [float]
        start - position where the gradient starts [float]
        stop - position where gradient stops [float]
        val - value at stop [float]
    output:
        value if start < x < stop:
        0 else
    """

    X = (x - start)
    if start < x < stop:
        return(X/(stop - start)*val)
    else:
        return(0)

def triangle(x, mid, width, val = 100):
    """
    input:
        x - current position [float]
        mid - position of triangle [float]
        width - width of triangle [float]
        val - max value of triangle [float]
    returns
        triangle value at x
    """

    temp_arr = np.copy(x)
    X = np.abs(np.abs(x - mid) - width)
    inds1 = np.nonzero(np.abs(x - mid) <= width)[0]
    inds2 = np.nonzero(np.abs(x - mid) > width)[0]
    temp_arr[inds1] = X[inds1]/width*val
    temp_arr[inds2] = 0
    return(temp_arr)

def limited_sin(x, mid, width, freq, amp = 100):
    """
    takes time series x and adds a sine between mid-width and mid+width
    input:
        x - current position [float]
        mid - position of sine[float]
        width - width of sine [float]
        freq - frequency of sine signal
        amp - amplitude of sine signal
    returns
        temp_arr - array with sine signal.
    """

    temp_arr = np.zeros_like(x)
    X = np.abs(np.abs(x - mid) - width)
    inds1 = np.nonzero(np.abs(x - mid) <= width)[0]
    inds2 = np.nonzero(np.abs(x - mid) > width)[0]
    sine = np.sin(x[inds1]*2*np.pi*freq)*amp
    temp_arr[inds1] = sine
    temp_arr[inds2] = 0
    return(temp_arr)

def triangle_combiner(x, widths, mids, vals):
    """
    Takes lists of triangle parameters and
    adds them together in a time series
    """
    temp_arr = np.zeros_like(x)
    for i in range(len(widths)):
        curr_triang = triangle(x, mids[i], widths[i], vals[i])
        curr_triang = pro.meanie(curr_triang, 10)
        temp_arr = temp_arr +  curr_triang
    return(temp_arr)

def triple_test():
    minfreq = 0.025
    maxfreq = 0.05
    fs = 2
    t = 5000
    n = int(fs*t)
    times = np.linspace(0, t, n)
    freqs = np.linspace(-fs/2, fs/2, n)
    widths = np.array([2, 3, 6, 9, 12])*25
    mids = np.array([100, 300, 500, 700, 900])*5
    vals = np.array([10, 1, 1, 1, 1])
    # widths = np.array([1])
    # mids = np.array([50])
    # vals = np.array([1])
    triangles = triangle_combiner(times, widths, mids, vals)

    plt.figure(0)
    plt.plot(times, triangles)
    plt.xlabel("Time [s]")
    plt.ylabel("Data")
    plt.title("Synthetic data")
    Freqs, Times, ffts = pro.fft_time_holes(triangles, times, n = 200, fs = fs)
    ffts += 1e-16

    plt.figure(1)
    plt.xlabel("Time [s]")
    plt.ylabel("Frequency [Hz]")
    plt.title("Time-Frequency fourier coefficients")
    plt.pcolormesh(Times, Freqs, np.log10(np.abs(ffts)), cmap = "gist_ncar")
    plt.colorbar()

    times, fourier_int = pro.fft_time_holes_integral(triangles, times, n = 200, fs = fs,\
                                                    minfreq = minfreq, maxfreq = maxfreq)

    plt.figure(2)
    plt.xlabel("Time [s]")
    plt.ylabel("Integrated fourier coeffiencts")
    plt.title("integral from f = %g to f = %g" % (minfreq, maxfreq))
    plt.plot(times, fourier_int)
    plt.show()

def triangle_plot():
    fs = 2
    t = 100
    n = int(fs*t)
    times = np.linspace(0, t, n)
    freqs = np.linspace(-fs/2, fs/2, n)
    width = 10
    mid = 50
    val = 1
    ys = triangle(times, mid, width, val)


    plt.plot(times, ys)
    plt.title("Triangle pulse")
    plt.ylabel("Measured unit")
    plt.xlabel("Time [s]")
    plt.show()

    fft = np.fft.fft(ys)/n
    fft = np.roll(fft, int(n/2))
    plt.plot(freqs, np.log10(np.abs(fft)))
    plt.xlabel("frequency [Hz]")
    plt.ylabel("log10(PSD)")
    plt.title("Power spectrum density")
    plt.show()
    

def index_test():
    fs = 2
    t = 3000
    n = int(fs*t)
    times = np.linspace(0, t, n)
    freqs = np.linspace(-fs/2, fs/2, n)
    widths = np.arange(5, 1500)
    indices = np.zeros(len(widths))
    for i in range(len(widths)):
        width1 = widths[i]
        width2 = 100
        mid = 1500
        val = 1
        ys1 = triangle(times, mid, width1, val)
        ys2 = triangle(times, mid, width2, val)#*(width1/widths[-1]))
    
        """
        plt.plot(times, ys1)
        plt.plot(times, ys2)
        plt.title("Triangle pulse")
        plt.ylabel("Measured unit")
        plt.xlabel("Time [s]")
        plt.legend(["Triangle A", "Triangle B"])
        plt.show()
        """
        l = 0
        
        freqs = np.linspace(-fs/2, fs/2, n)[int(n/2):]
        minfreq = 0.1 + l
        maxfreq = 1 + l
        fft1 = np.fft.fft(ys1)[:int(n/2)]/n*2
        fft2 = np.fft.fft(ys2)[:int(n/2)]/n*2
        fft1 = np.abs(fft1)
        fft2 = np.abs(fft2)
        """
        plt.plot(freqs, np.log10(fft1))
        plt.plot(freqs, np.log10(fft2))
        plt.title("Logarithmic PSD")
        plt.legend(["Triangle A", "Triangle B"])
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Log10(PSD)")
        plt.show()
        """        
        
        df = freqs[1] - freqs[0]
        N_maxfreq = int(maxfreq/df)
        N_minfreq = int(minfreq/df)
        
        sum1 = np.sum(fft1[N_minfreq:N_maxfreq]*df)
        sum2 = np.sum(fft2[N_minfreq:N_maxfreq]*df)
        
        index = (sum1 - sum2)/np.max([sum1, sum2])
        indices[i] = index

    witty = np.linspace(1, widths[-1]/width2)
    plt.plot(widths/width2, indices)
    #plt.plot(widths, 1 - widths/width2)
    #plt.plot([width2, width2], [-1, 1], "k")
    plt.plot(witty, -1 + 1/witty)
    plt.xlabel("Width of triangle A")
    plt.ylabel("$I_{A - B}$")
    plt.grid("on")
    #plt.legend(["Comparison Index", "1 - width$_1$/width$_2$"])
    plt.title("Comparison indices for varying triangle sizes")
    # plt.savefig("Figures/comp_ind_example.pdf")
    plt.show()
    
if __name__ == "__main__":
    index_test()