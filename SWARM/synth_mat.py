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

if __name__ == "__main__":
    fs = 2
    t = 5000
    n = int(fs*t)
    times = np.linspace(0, t, n)
    freqs = np.linspace(-fs/2, fs/2, n)
    widths = np.array([1, 3, 6, 9, 12])*25
    mids = np.array([100, 300, 500, 700, 900])*5
    vals = np.array([1, 1, 1, 1, 1])
    # widths = np.array([1])
    # mids = np.array([50])
    # vals = np.array([1])
    triangles = triangle_combiner(times, widths, mids, vals)

    plt.figure(0)
    plt.plot(times, triangles)
    Freqs, Times, ffts = pro.fft_time_holes(triangles, times, n = 200, fs = fs)
    ffts += 1e-16

    plt.figure(1)
    plt.pcolormesh(Times, Freqs, np.log10(np.abs(ffts)), cmap = "gist_ncar")
    plt.colorbar()

    times, fourier_int = pro.fft_time_holes_integral(triangles, times, n = 200, fs = fs,\
                                                    minfreq = 0.05, maxfreq = 0.1)

    plt.figure(2)
    plt.plot(times, fourier_int)
    plt.show()