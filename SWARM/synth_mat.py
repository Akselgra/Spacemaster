import numpy as np
import matplotlib.pyplot as plt
import SWARMprocess
pro = SWARMprocess.SWARMprocess()

fig_width_pt = 418.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inches
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*golden_mean       # height in inches
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

def square(x, mid, width, val = 100):
    """
    defines a square pulse signal
    """
    temp_arr = np.copy(x)
    inds1 = np.nonzero(np.abs(x-mid) <= width)[0]
    inds2 = np.nonzero(np.abs(x-mid) > width)[0]
    temp_arr[inds1] = val
    temp_arr[inds2] = 0
    return temp_arr
    
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

def trapez(x, mid, width, val = 100, grad = 0.1):
    """
    Defines a trapez
    """
    temp_arr = np.copy(x)
    X = np.abs(np.abs(x-mid)-width)
    
    gradind = int(val/grad)
    inds1 = np.nonzero(np.abs(x - mid) <= width)[0]
    inds2 = np.nonzero(np.logical_and(np.abs(x - mid) > width, np.abs(x-mid) < (width + gradind)))[0]
    inds3 = np.nonzero(np.abs(x - mid) > width + gradind)[0]
    
    X = np.abs(np.abs(x - mid) - (width + gradind))
    
    temp_arr[inds1] = val
    temp_arr[inds2] = X[inds2]*grad
    temp_arr[inds3] = 0

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
    
    fs = 1
    t = 5000
    f = 0.5
    n = int(fs*t)
    N = 200*fs
    times = np.linspace(0, t, n)
    freqs = np.linspace(-fs/2, fs/2, n)
    widths = np.array([2, 4, 6, 8, 10])*25
    mids = np.array([100, 300, 500, 700, 900])*5
    vals = np.array([1, 1, 1, 1, 1])
    minfreq = 0.025
    minfreq = 1/np.min(widths)
    maxfreq = 0.5
    # widths = np.array([1])
    # mids = np.array([50])
    # vals = np.array([1])
    triangles = triangle_combiner(times, widths, mids, vals)

    plt.figure(0)
    plt.plot(times, triangles)
    plt.xlabel("Time [s]")
    plt.ylabel("Data")
    plt.title("Synthetic data")
    plt.savefig("Figures/matfigs/synth/multi_triangle_pulse.pdf")
    Freqs, Times, ffts = pro.fft_time_holes(triangles, times, n = N, fs = fs)
    ffts += 1e-16

    ind = int(f*N/fs)
    print(ind)
    print(Freqs[:, :ind])
    plt.figure(1)
    plt.xlabel("Time [s]")
    plt.ylabel("Frequency [Hz]")
    plt.title("Time-Frequency fourier coefficients")
    plt.pcolormesh(Times[:, :ind], Freqs[:, :ind], np.log10(np.abs(ffts[:, :ind])), cmap = "gist_ncar")
    plt.colorbar(label = "Log10(PSD) [dB]")
    plt.savefig("Figures/matfigs/synth/multi_triangle_pulse_fft_time.pdf")

    times, fourier_int = pro.fft_time_holes_integral(triangles, times, n = N, fs = fs,\
                                                    minfreq = minfreq, maxfreq = maxfreq)

    plt.figure(2)
    plt.xlabel("Time [s]")
    plt.ylabel("Integrated fourier coeffiencts")
    plt.title("integral from f = %g to f = %g" % (minfreq, maxfreq))
    plt.plot(times, fourier_int)
    plt.savefig("Figures/matfigs/synth/multi_triangle_pulse_fft_time_integral.pdf")
    plt.show()

def triangle_plot():
    fs = 1000
    t = 100
    n = int(fs*t)
    times = np.linspace(0, t, n)
    df = fs/n
    freqs = np.arange(n)*df
    # freqs = np.linspace(-fs/2, fs/2, n)
    width = 10
    mid = 50
    val = 1
    ys = triangle(times, mid, width, val)


    plt.plot(times, ys)
    plt.title("Triangle pulse")
    plt.ylabel("Amplitude")
    plt.xlabel("Time [s]")
    plt.savefig("Figures/matfigs/synth/triangle_pulse.pdf")
    plt.show()

    fft = np.fft.fft(ys)/fs
    # fft = np.roll(fft, int(n/2))
    # plt.plot(freqs, np.log10(np.abs(fft)))
    # plt.xlabel("frequency [Hz]")
    # plt.ylabel("log10(PSD)")
    # plt.title("Power spectrum density")
    # plt.show()
    
    f = 0.5
    positive_freqs = freqs[:int(n/2)]
    positive_ffts = fft[:int(n/2)]
    positive_freqs = positive_freqs[:int(f*t)]
    positive_ffts = positive_ffts[:int(f*t)]
    
    

    plt.plot(positive_freqs, (np.abs(positive_ffts)))
    plt.plot(positive_freqs, (width*np.sinc(width*positive_freqs)**2), ".")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("PSD")
    plt.legend(["PSD of triangle", "Sinc$^2$"])
    plt.title("PSD of triangle pulse")
    plt.savefig("Figures/matfigs/synth/triangle_pulse_fft.pdf")
    plt.show()
    
    print(freqs[1] - freqs[0])
    print(fs/(n-1))

def triangle_plot2():
    fs = 1000
    t = 100
    n = int(fs*t)
    times = np.linspace(0, t, n)
    df = fs/n
    freqs = np.arange(n)*df
    # freqs = np.linspace(-fs/2, fs/2, n)
    width1 = 10
    width2 = 20
    mid = 50
    val = 1
    ys1 = triangle(times, mid, width1, val)
    ys2 = triangle(times, mid, width2, val)

    plt.plot(times, ys1)
    plt.plot(times, ys2)
    plt.title("Triangle pulses")
    plt.ylabel("Amplitude")
    plt.xlabel("Time [s]")
    plt.legend(["Triangle A", "Triangle B"])
    plt.savefig("Figures/matfigs/synth/double_triangle_pulse.pdf")
    plt.show()

    fft1 = np.fft.fft(ys1)/fs
    fft2 = np.fft.fft(ys2)/fs
    # fft = np.roll(fft, int(n/2))
    # plt.plot(freqs, np.log10(np.abs(fft)))
    # plt.xlabel("frequency [Hz]")
    # plt.ylabel("log10(PSD)")
    # plt.title("Power spectrum density")
    # plt.show()
    
    f = 0.5
    positive_freqs = freqs[:int(n/2)]
    positive_ffts1 = fft1[:int(n/2)]
    positive_ffts2 = fft2[:int(n/2)]
    positive_freqs = positive_freqs[:int(f*t)]
    positive_ffts1 = positive_ffts1[:int(f*t)]
    positive_ffts2 = positive_ffts2[:int(f*t)]
    
    

    plt.plot(positive_freqs, (np.abs(positive_ffts1)))
    plt.plot(positive_freqs, np.abs(positive_ffts2))
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("PSD")
    plt.legend(["PSD of triangle A", "PSD of triangle B"])
    plt.title("PSD of triangle pulses")
    plt.savefig("Figures/matfigs/synth/double_triangle_pulse_fft.pdf")
    plt.show()
    
    minfreq = 0.1
    minfreq = 1/width1
    print(minfreq)
    maxfreq = 0.5
    i_minfreq = int(minfreq/df)
    i_maxfreq = int(maxfreq/df)
    
    sumA = np.sum(np.abs(positive_ffts1[i_minfreq:i_maxfreq])*df)
    sumB = np.sum(np.abs(positive_ffts2[i_minfreq:i_maxfreq])*df)
    comp_ind = (sumA - sumB)/(np.max([sumA, sumB]))
    print(comp_ind)
    
def index_test():
    grad = 0.01
    fs = 2
    t = 3000
    n = int(fs*t)
    times = np.linspace(0, t, n)
    freqs = np.linspace(-fs/2, fs/2, n)
    widths = np.arange(5, 500)
    grads = np.copy(widths)/widths[-1]*5*grad
    # print(grads[-1])

    indices = np.zeros(len(widths))
    for i in range(len(widths)):
        width1 = widths[i]
        width2 = 100
        mid = 1500
        val = 1
        # grad1 = grads[i]
        # grad2 = grad
        ys1 = triangle(times, mid, width1, val)
        ys2 = triangle(times, mid, width2, val)#*(width1/widths[-1]))
        # ys1 = square(times, mid, width1, val)
        # ys2 = square(times, mid, width1, val)
        # ys1 = limited_sin(times, mid, width1, freq = 0.5, amp = val)
        # ys2 = limited_sin(times, mid, width2, freq = 0.5, amp = val)
        # ys1 = trapez(times, mid, width1, val, grad = grad1)
        # ys2 = trapez(times, mid, width2, val, grad = grad2)
        
    
        
        # plt.plot(times, ys1)
        # plt.plot(times, ys2)
        # plt.title("Triangle pulse")
        # plt.ylabel("Measured unit")
        # plt.xlabel("Time [s]")
        # plt.legend(["Triangle A", "Triangle B"])
        # plt.show()
        
        l = 0
        
        freqs = np.linspace(-fs/2, fs/2, n)[int(n/2):]
        minfreq = 0.1 + l
        maxfreq = 1 + l
        fft1 = np.fft.fft(ys1)[:int(n/2)]/n*2
        fft2 = np.fft.fft(ys2)[:int(n/2)]/n*2
        fft1 = np.abs(fft1)
        fft2 = np.abs(fft2)
        
        # plt.plot(freqs, np.log10(fft1))
        # plt.plot(freqs, np.log10(fft2))
        # plt.title("Logarithmic PSD")
        # plt.legend(["Triangle A", "Triangle B"])
        # plt.xlabel("Frequency [Hz]")
        # plt.ylabel("Log10(PSD)")
        # plt.show()
             
        
        df = freqs[1] - freqs[0]
        N_maxfreq = int(maxfreq/df)
        N_minfreq = int(minfreq/df)
        
        sum1 = np.sum(fft1[N_minfreq:N_maxfreq]*df)
        sum2 = np.sum(fft2[N_minfreq:N_maxfreq]*df)
        
        index = (sum1 - sum2)/np.max([sum1, sum2])
        indices[i] = index
    
    plt.plot(times, ys1)
    plt.plot(times, ys2)
    plt.show()

    witty = np.linspace(1, widths[-1]/width2)
    plt.plot(widths/width2, indices)
    #plt.plot(widths, 1 - widths/width2)
    #plt.plot([width2, width2], [-1, 1], "k")
    # plt.plot(witty, -1 + 1/witty)
    plt.xlabel("Width of triangle A")
    plt.ylabel("$I_{A - B}$")
    plt.grid("on")
    #plt.legend(["Comparison Index", "1 - width$_1$/width$_2$"])
    plt.title("Comparison indices for varying triangle sizes")
    # plt.savefig("Figures/comp_ind_example.pdf")
    plt.show()

def index_test_trapez():
    grad = 0.01
    fs = 2
    t = 3000
    n = int(fs*t)
    times = np.linspace(0, t, n)
    freqs = np.linspace(-fs/2, fs/2, n)
    widths = np.arange(5, 500)
    grads = np.copy(widths)/widths[-1]*5*grad
    # print(grads[-1])

    indices = np.zeros(len(widths))
    for i in range(len(widths)):
        width2 = widths[i]
        # width2 = 100
        width1 = 100
        mid = 1500
        val = 1
        # grad2 = grads[i]
        grad2 = grads[i]
        # grad2 = grad
        grad1 = grad
        
        ys1 = trapez(times, mid, width1, val, grad = grad1)
        ys2 = trapez(times, mid, width2, val, grad = grad2)
        # ys1 = square(times, mid, width1)
        # ys2 = square(times, mid, width2)
        
    
        
        # plt.plot(times, ys1)
        # plt.plot(times, ys2)
        # plt.title("Trapez pulse")
        # plt.ylabel("Measured unit")
        # plt.xlabel("Time [s]")
        # plt.legend(["Trapez A", "Trapez B"])
        # plt.show()
        
        l = 0
        
        freqs = np.linspace(-fs/2, fs/2, n)[int(n/2):]
        minfreq = 0.1 + l
        maxfreq = 1 + l
        fft1 = np.fft.fft(ys1)[:int(n/2)]/n*2
        fft2 = np.fft.fft(ys2)[:int(n/2)]/n*2
        fft1 = np.abs(fft1)
        fft2 = np.abs(fft2)
        
        # plt.plot(freqs, np.log10(fft1))
        # plt.plot(freqs, np.log10(fft2))
        # plt.title("Logarithmic PSD")
        # plt.legend(["Triangle A", "Triangle B"])
        # plt.xlabel("Frequency [Hz]")
        # plt.ylabel("Log10(PSD)")
        # plt.show()
             
        
        df = freqs[1] - freqs[0]
        N_maxfreq = int(maxfreq/df)
        N_minfreq = int(minfreq/df)
        
        sum1 = np.sum(fft1[N_minfreq:N_maxfreq]*df)
        sum2 = np.sum(fft2[N_minfreq:N_maxfreq]*df)
        
        index = (sum1 - sum2)/np.max([sum1, sum2])
        indices[i] = index
    
    plt.plot(times, ys1)
    plt.plot(times, ys2)
    plt.show()
    
    plt.plot(freqs, np.log10(fft1))
    plt.plot(freqs, np.log10(fft2))
    plt.title("Logarithmic PSD")
    plt.legend(["Trapez A", "Trapez B"])
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Log10(PSD)")
    plt.show()

    witty = np.linspace(1, widths[-1]/width2)
    plt.plot(grads/grad1, indices)
    #plt.plot(widths, 1 - widths/width2)
    #plt.plot([width2, width2], [-1, 1], "k")
    # plt.plot(witty, -1 + 1/witty)
    plt.xlabel("slope of B / slope of A")
    plt.ylabel("$I_{A - B}$")
    plt.grid("on")
    #plt.legend(["Comparison Index", "1 - width$_1$/width$_2$"])
    plt.title("Comparison indices for varying trapez slopes")
    # plt.savefig("Figures/comp_ind_example.pdf")
    plt.show()
def sine_plot1():
    fs = 1000
    t = 1/0.04
    f1 = 0.04
    f_lim = 0.5
    n = int(fs*t)
    df = fs/n
    times = np.linspace(0, t, n)
    signal = np.sin(2*np.pi*f1*times)
    
    plt.plot(times, signal)
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.title("Sine signal")
    plt.savefig("Figures/signal_example.pdf")
    plt.show()
    
    # freqs = np.linspace(-fs/2, fs/2, n)
    freqs = np.arange(n)*df
    ffts = np.fft.fft(signal)/n*2
    # ffts = np.roll(ffts, int(n/2))

    
    
    freqs = freqs[:int(n/2)]
    ffts = ffts[:int(n/2)]
    
    freqs = freqs[:int(f_lim*t)]
    ffts = ffts[:int(f_lim*t)]
    

    plt.plot(freqs, (np.abs(ffts)), "-o")
    plt.plot(freqs, (np.imag(ffts)), "-o")
    plt.plot(freqs, (np.real(ffts)), "-o")
    plt.title("Fourier transform")
    plt.ylabel("Fourier coefficients")
    plt.xlabel("Frequency [Hz]")
    plt.legend(["Power Spectral Density", "Complex component", "Real component"])
    plt.xticks(np.arange(7)*0.04*2 - 0.04)
    plt.savefig("Figures/fft_example.pdf")
    plt.show()

    

def sine_plot2():
    fs = 1000
    t = 4/0.04
    f1 = 0.04
    f2 = 0.1
    f3 = 0.2
    f_lim = 0.5
    n = int(fs*t)
    df = fs/n
    times = np.linspace(0, t, n)
    
    signal = 3*np.cos(2*np.pi*f1*times) + np.cos(2*np.pi*f2*times) + np.cos(2*np.pi*f3*times)
    
    plt.plot(times, signal)
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.title("A combination of cosines")
    plt.savefig("Figures/signal_example_multiple.pdf")
    plt.show()
    
    # freqs = np.linspace(-fs/2, fs/2, n)
    freqs = np.arange(n)*df
    ffts = np.fft.fft(signal)/n*2
    # ffts = np.roll(ffts, int(n/2))

    
    
    freqs = freqs[:int(n/2)]
    ffts = ffts[:int(n/2)]
    
    freqs = freqs[:int(f_lim*t)]
    ffts = ffts[:int(f_lim*t)]
    

    plt.plot(freqs, (np.abs(ffts)), "-o")
    plt.title("PSD")
    plt.ylabel("PSD")
    plt.xlabel("Frequency [Hz]")
    plt.xticks(np.arange(7)*0.04*2 - 0.04)
    plt.savefig("Figures/fft_example_multiple.pdf")
    plt.show()
    
if __name__ == "__main__":
    triple_test()
    