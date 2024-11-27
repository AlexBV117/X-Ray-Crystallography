import matplotlib.pyplot as plt
import numpy as np
from pandas import read_csv 
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from math import sin, radians, ceil, log10, floor


#===========================#
#   Function Definitions    #
#===========================#
# For getting the wavelength of the alpha and beta lines
def linear_fit(x, m, c):
    return m * x + c

def round_sig(x, sig=1):
    return round(x, sig-int(floor(log10(abs(x))))-1)

def propagate_error(new_val, old_val, old_err):
    percent = (old_err / old_val) * 100
    new_err = (percent / 100) * new_val
    return new_err

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
    ax1.axvline(x=peak, linestyle='dashed', c='k')
ax1.plot(angles, scan1, marker=".", color=(0.5, 0.5, 0.5), label="35kV Scan", linestyle="solid")
ax1.plot(angles, scan2, marker=".",  color=(0.75, 0.75, 0.75), label="25kV Scan", linestyle="solid")
ax1.set_title("Bragg reflection at an LiF monocrystal")
ax1.set_xlabel("Angle θ (°)")
ax1.set_ylabel("Counts per second (10 second average)")
ax1.set_yscale('log')
ax1.legend(loc="upper right")
fig1.savefig("Bragg_Scan_2.png")
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

a_0 = 402.7e-12
alpha_wavelen = a_0 / fit_alpha[0]
beta_wavelen = a_0 / fit_beta[0]
alpha_wavelen_err = round_sig(propagate_error(alpha_wavelen, fit_alpha[0], alpha_error[0]))
beta_wavelen_err = round_sig(propagate_error(beta_wavelen, fit_beta[0], beta_error[0]))
# Plot the diffraction order against the sin of the angle
fig2, ax2 = plt.subplots()
ax2.plot(fit_x, fit_y_alpha, c='k', linestyle="dashed", label=f"Alpha line: {alpha_wavelen:.3g} ± {alpha_wavelen_err} m")
ax2.plot(fit_x, fit_y_beta, c='k', linestyle="dotted",  label=f"Beta line:  {beta_wavelen:.3g} ± {beta_wavelen_err} m")
for i, peak in enumerate(peaks):
    n = ceil((i+1)/2) # what is the diffraction order 
    ax2.scatter(sin(radians(peak)), n, c='k', marker="D")
ax2.set_title("Diffraction order as a function of the sin of the angle")
ax2.set_xlabel("Sin(θ)")
ax2.set_ylabel("Diffraction Order 'n'")
ax2.legend(loc="upper right")
fig2.savefig("Alpha_Beta_line_2.png")

plt.show()
