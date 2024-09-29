import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from math import ceil, sin, radians

data = read_csv("Unknown_Scan.csv")
angles = data.get('&b / °').to_numpy()
scan = data.get('R_0 / 1/s').to_numpy()
angles2 = [sin(radians(ang)) for ang in angles]

alpha_line = 71.08 # n*lambda/pm
beta_line = 63.06  # n=1
angstrom = 100 #pm

peaks_index = find_peaks(scan, prominence=10, height=1, width=(2, 5))[0]
peaks = angles[peaks_index]

plt1, ax1 = plt.subplots()
ax1.plot(angles, scan, c='b')
for peak in peaks:
    ax1.axvline(peak, linestyle='dashed', c='r')
ax1.set_title("Bragg Scan of an Unknown Sample")
ax1.set_xlabel("Angle θ (°)")
ax1.set_ylabel("Counts per second (10 second average)")
ax1.set_yscale("log")
plt1.savefig("Bragg_Scan.png")

def linear_fit(x, m, c):
    return m * x + c

value_pairs = [[], []]

for i, peak in enumerate(peaks):
    n = ceil((i+1)/2)
    if i % 2 == 0:
        value_pairs[0].append(sin(radians(peak)))
        value_pairs[1].append(n*beta_line)
    else:
        value_pairs[0].append(sin(radians(peak)))
        value_pairs[1].append(n*alpha_line)

[fit, cov] = curve_fit(linear_fit, value_pairs[0], value_pairs[1], p0=[565.2, 0])

lattice = fit[0] / angstrom
error = np.sqrt(np.diag(cov))[0] / angstrom


fit_line = [linear_fit(ang, *fit) for ang in angles2]
plt2, ax2 = plt.subplots()
ax2.scatter(value_pairs[0], value_pairs[1], c='r')
ax2.plot(angles2, fit_line, linestyle="dotted", c='b', label=f"Fit line slope: {lattice:.2f} ± {error:.2f} Å")
ax2.set_title("Value pairs as a function of sin θ")
ax2.set_xlabel("sin θ")
ax2.set_ylabel("nλ/pm")
ax2.legend(loc="upper left") 
plt2.savefig("value_Pairs.png")

plt.show()
