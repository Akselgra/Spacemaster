import numpy as np


def standard_form(f, dec = 2):
    """
    Takes Float
    Returns float on standard form as string
    """
    
    zeros = 0
    string = str(f)
    if f < 0:
        string = string[1:]
    
    for i in range(len(string)-1):
        if string[i+1] == ".":
            break
        zeros += 1
    
    for i in range(len(string[zeros:])):
        if string[i] == ".":
            continue
        if string[i] != "0":
            break
        zeros -= 1
     
    dot_ind = 0
    for i in range(len(string)):
        if string[i] == ".":
            break
        dot_ind += 1
    
    
    if zeros >= 0:
        if 1 <= dot_ind < dec+2:
            temp = string[1:dot_ind] + string[dot_ind+1 : dec+1] + str(int((np.round(float(string[dec+1] + "." + string[dec+2:])))))
        else:
            temp = string[1:dec+1]
            
        bob = string[0] + "." + temp
        if f < 0:
            bob = "-" + bob
    
    else:
        temp = string[np.abs(zeros) + 2 : np.abs(zeros) + dec + 1] + str(int((np.round(float(string[np.abs(zeros) + dec+1] + "." + string[np.abs(zeros) + dec+2:])))))
        bob = string[np.abs(zeros)+1] + "." + temp
        if f < 0:
            bob = "-" + bob
    
    return(bob + " \\cdot 10^{%g}" % zeros)
    

a1 = 0.00024
a2 = 1
a3 = 1
a4 = 1

b1 = 1
b2 = 1
b3 = 1
b4 = 1

f_a1 = standard_form(a1, dec = 1)

print("\\begin{table}[htbp]")
print("\\centering")
print("\\caption{Table of regression coefficients for the linear regressions found in Figure \\ref{fig:multi_std}. The first column has the latitudinal regions, the second column has the regression slopes and the third column has the regression constants.}")
print("\\begin{tabular}{|c|c|c|}")
print("\\hline")
print("Region & Std regression slope [s$^{-1}$] & Std regression constant [s$^{-1}$]\\" + "\\")
print("\\hline")
print("Equatorial & $" + f_a1 +  " \pm 8.1 \cdot 10^{-5}$ & $1.2 \cdot 10^{-1} \pm 7.2 \cdot 10^{-3}$\\" + "\\")
print("\\hline ")
print("Midlatitude & $-8.6 \cdot 10^{-5} \pm 6.5 \cdot 10^{-5}$ & $1.9 \cdot 10^{-1} \pm 5.8 \cdot 10^{-3}$\\" + "\\" )
print("\\hline ")
print("Auroral Oval & $3.5 \cdot 10^{-4} \pm 9.0 \cdot 10^{-5}$ & $1.7 \cdot 10^{-1} \pm 8.0 \cdot 10^{-3}$\\" + "\\")
print("\\hline ")
print("Northern polar cap & $6.7 \cdot 10^{-4} \pm 7.4 \cdot 10^{-5}$ & $1.6 \cdot 10^{-1} \pm 6.6 \cdot 10^{-3}$\\" + "\\")
print("\\hline ")
print("\\end{tabular}")
print("\\label{tab:multi_std}")
print("\\end{table}")
