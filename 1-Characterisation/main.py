import matplotlib.pyplot as plt
import numpy as np
from pandas import read_csv 
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from math import sin, radians, ceil

#===========================#
#   Function Definitions    #
#===========================#
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

# convert from wavelength to energy
def to_eV(wave_len):
    c   = 299792458 # ms_1              (https://en.wikipedia.org/wiki/Speed_of_light)
    a_0 = 402.7e-12 # pm                (LEYBOLD Physics Leaflets P7.1.2.1)
    h   = 4.135667696e-15 # eV Hz^-1    (https://en.wikipedia.org/wiki/Planck_constant)
    return (c / (a_0 / wave_len) * h)

# Load the data and covert it to a numpy array
data = read_csv('25kV_&_35kV_LiF_Scan.csv')
angles = data.get('&b / °').to_numpy()
scan1 = data.get('R_0 / 1/s').to_numpy()
scan2 = data.get('R_1 / 1/s').to_numpy()


#===============#
#   Raw Data    #
#===============#
# locate the peaks associated with the alpha and beta lines
peaks_index= find_peaks(scan1, prominence=10, height=1, width=(2, 5))[0]
peaks = angles[peaks_index]

# plot the raw data
fig1, ax1 = plt.subplots()
for peak in peaks:
    # add vertical lines for the peaks
    ax1.axvline(x=peak, linestyle='dashed', c='r')
ax1.plot(angles, scan1, c='b', label="35kV Scan")
ax1.plot(angles, scan2, c='g', label="25kV Scan")
ax1.set_title("Bragg reflection at an LiF monocrystal")
ax1.set_xlabel("Angle θ (°)")
ax1.set_ylabel("Counts per second (10 second average)")
ax1.set_yscale('log')
ax1.legend(loc="upper right")
fig1.savefig("Bragg_Scan.png")
#=========================================#
#   Wavelength of alpha and beta lines    #
#=========================================#
# get sine of angle for future processing
sin_peaks = [sin(radians(i)) for i in peaks]

# Split the values for the alpha and the beta lines
alpha_lines = sin_peaks[1::2]
beta_lines = sin_peaks[0::2]

# Get the wavlength for the alpha and beta lines
[fit_alpha, cov_alpha] = curve_fit(linear_fit, alpha_lines, [1, 2, 3], p0=[17.445, 0]) # inital values for m bassed off lab script
[fit_beta, cov_beta] = curve_fit(linear_fit, beta_lines, [1, 2, 3], p0=[19.654, 0])

# compute one standard deviation errors
alpha_error = np.sqrt(np.diag(cov_alpha))
beta_error = np.sqrt(np.diag(cov_beta))

# Plot the fitted vaues
fit_x = np.arange(0, 0.6, 0.05)
fit_y_alpha = linear_fit(fit_x, *fit_alpha)
fit_y_beta = linear_fit(fit_x, *fit_beta)

# Plot the diffraction order against the sin of the angle
fig2, ax2 = plt.subplots()
for i, peak in enumerate(peaks):
    n = ceil((i+1)/2) # what is the diffraction order 
    ax2.scatter(sin(radians(peak)), n, c='r')
ax2.plot(fit_x, fit_y_alpha, c='b', label=f"Alpha line: {to_eV(fit_alpha[0]):,.2f} ± {alpha_error[0]:.2f} eV")
ax2.plot(fit_x, fit_y_beta, c='g', label=f"Beta line: {to_eV(fit_beta[0]):,.2f} ± {beta_error[0]:.2f} eV")
ax2.set_title("Diffraction order as a function of the sin of the angle")
ax2.set_xlabel("Sin(θ)")
ax2.set_ylabel("Diffraction Order 'n'")
ax2.legend(loc="upper right")
fig2.savefig("Alpha_Beta_line.png")
#==========================#
#   igh-energy cut-off     #
#==========================#
# Create a subset of the data where 'n' is assumed to be 1
alpha1, beta2 = peaks_index[1], peaks_index[2]
end_index = alpha1 + ceil((beta2 - alpha1)/2)
s1FO_y = scan1[0: end_index]
s2FO_y = scan2[0: end_index]
s1FO_x = angles[0: end_index]
s2FO_x = angles[0: end_index]

# grab the tail with the high energy cut off (reagon that looks mostly linear on the right)
s1FO_y_fit_range = scan1[20: 30] # indexes determined visully
s2FO_y_fit_range = scan2[40: 50]
s1FO_x_fit_range = angles[20: 30]
s2FO_x_fit_range = angles[40: 50]
# convert to energy
s1FO_x = [ to_energy(i, 1) for i in s1FO_x]
s2FO_x = [ to_energy(i, 1) for i in s2FO_x]
s1FO_x_fit_range = [ to_energy(i, 1) for i in s1FO_x_fit_range]
s2FO_x_fit_range = [ to_energy(i, 1) for i in s2FO_x_fit_range]

[s1FO_fit, cov_s1FO] = curve_fit(linear_fit, s1FO_x_fit_range, s1FO_y_fit_range)
[s2FO_fit, cov_s2FO] = curve_fit(linear_fit, s2FO_x_fit_range, s2FO_y_fit_range)
s1FO_x_fit = np.arange(26_370, 36_370, 10)
s1FO_y_fit = linear_fit(s1FO_x_fit, *s1FO_fit)
s2FO_x_fit = np.arange(20_890, 25_890, 10)
s2FO_y_fit = linear_fit(s2FO_x_fit, *s2FO_fit)

s1FO_HECO = ((0 - s1FO_fit[1])/s1FO_fit[0])
s2FO_HECO = ((0 - s2FO_fit[1])/s2FO_fit[0])
# Plot the high energy cut off
fig3, ax3 = plt.subplots()
ax3.plot(s1FO_x, s1FO_y, label=f"(35k eV) high-energy cut-off: {s1FO_HECO:.0f} eV", c='b')
ax3.plot(s1FO_x_fit, s1FO_y_fit, c='r')
ax3.plot(s2FO_x, s2FO_y, label=f"(25k eV) high-energy cut-off: {s2FO_HECO:.0f} eV", c='g')
ax3.plot(s2FO_x_fit, s2FO_y_fit, c='r')
ax3.set_title("Count rate as a function of X-Ray energy")
ax3.set_xlabel("Energy (eV)")
ax3.set_ylabel("Counts per second (10 second average)")
ax3.legend(loc="upper right")
fig3.savefig("X-Ray_Energy.png")
