import matplotlib.pyplot as plt
import numpy as np
from pandas import read_csv 
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from math import sin, radians, ceil

# Load the data and covert it to a numpy array
data = read_csv('25kV_&_35kV_LiF_Scan.csv')
angles = data.get('&b / °').to_numpy()
scan1 = data.get('R_0 / 1/s').to_numpy()
scan2 = data.get('R_1 / 1/s').to_numpy()

# locate the peaks associated with the alpha and beta lines
peaks_index= find_peaks(scan1, prominence=10, height=1, width=(2, 5))[0]
peaks = angles[peaks_index]

# plot the raw data
_, ax1 = plt.subplots()
for peak in peaks:
    ax1.axvline(x=peak, linestyle='dashed', c='r')
ax1.plot(angles, scan1, c='b', label="35kV Scan")
ax1.plot(angles, scan2, c='g', label="25kV Scan")
ax1.set_title("Bragg reflection at an LiF monocrystal")
ax1.set_xlabel("Angle θ (°)")
ax1.set_ylabel("Average X-Ray Count /s over 10 seconds")
ax1.set_yscale('log')
ax1.legend(loc="upper right")

# Create a subset of the data where 'n' is assumed to be 1
alpha1, beta2 = peaks_index[1], peaks_index[2]
end_index = alpha1 + ceil((beta2 - alpha1)/2)
first_order_y = scan1[0: end_index]
first_order_x = angles[0: end_index]
# grab the tail of with the high energy cut off (reagon that looks mostly linear)
first_order_y_fit_range = scan1[22: 30] # indexes setermined visully
first_order_x_fit_range = angles[22: 30]

# For getting the wavelength of the alpha and beta lines
def linear_fit(x, m, c):
    return m * x + c

# convert from angle to energy
def to_energy(ang, n):
    c   = 299792458 # ms_1              (https://en.wikipedia.org/wiki/Speed_of_light)
    a_0 = 402.7e-12 # pm                (LEYBOLD Physics Leaflets P7.1.2.1)
    h   = 4.135667696e-15 # eV Hz^-1    (https://en.wikipedia.org/wiki/Planck_constant)
    wave_len = (a_0 / n) * sin(radians(ang))
    wave_feq = c / wave_len
    return wave_feq * h

first_order_x = [ to_energy(i, 1) for i in first_order_x]
first_order_x_fit_range = [ to_energy(i, 1) for i in first_order_x_fit_range]

[first_order_fit, first_order_cov] = curve_fit(linear_fit, first_order_x_fit_range, first_order_y_fit_range)
first_order_x_fit = np.arange(20_000, 40_000, 100)
first_order_y_fit = linear_fit(first_order_x_fit, *first_order_fit)
print(first_order_fit)
# add vertical lines for the peaks
_, ax2 = plt.subplots()
for i, peak in enumerate(peaks):
    n = ceil((i+1)/2)
    ax2.scatter(sin(radians(peak)), n, c='r')

# get sine of angle for future processing
peaks = [sin(radians(i)) for i in peaks]


# convert from wavelength to energy
def to_eV(wave_len):
    c   = 299792458 # ms_1              (https://en.wikipedia.org/wiki/Speed_of_light)
    a_0 = 402.7e-12 # pm                (LEYBOLD Physics Leaflets P7.1.2.1)
    h   = 4.135667696e-15 # eV Hz^-1    (https://en.wikipedia.org/wiki/Planck_constant)
    return (c / (a_0 / wave_len) * h)

# Split the values for the alpha and the beta lines
alpha_lines = peaks[1::2]
beta_lines = peaks[0::2]

# Get the wavlength for the alpha and beta lines
[fit_alpha, cov_alpha] = curve_fit(linear_fit, alpha_lines, [1, 2, 3], p0=[17.445, 0])
[fit_beta, cov_beta] = curve_fit(linear_fit, beta_lines, [1, 2, 3], p0=[19.654, 0])



# Plot the fitted vaues
fit_x = np.arange(0, 0.6, 0.05)
fit_y_alpha = linear_fit(fit_x, *fit_alpha)
fit_y_beta = linear_fit(fit_x, *fit_beta)
ax2.plot(fit_x, fit_y_alpha, c='b', label=f"Alpha line: {to_eV(fit_alpha[0])} eV")
ax2.plot(fit_x, fit_y_beta, c='g', label=f"Beta line: {to_eV(fit_beta[0])} eV")
ax2.legend(loc="upper right")

# Plot the high energy cut off
_, ax3 = plt.subplots()
ax3.plot(first_order_x, first_order_y)
ax3.plot(first_order_x_fit, first_order_y_fit)
ax3.plot(first_order_x_fit_range, first_order_y_fit_range)
plt.show()
