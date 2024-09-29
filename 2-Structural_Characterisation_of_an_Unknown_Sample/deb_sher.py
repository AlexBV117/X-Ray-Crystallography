import matplotlib.pyplot as plt
from pandas import read_csv
from scipy.signal import find_peaks
from math import isnan

data1 = read_csv("NaCl_Powder_Scan_1.csv")
angles1 = data1.get('&b / °').to_numpy()
scan1 = data1.get('R_0 / 1/s').to_numpy()

data2 = read_csv("NaCl_Powder_111_Scan_2.csv")
angles2 = data2.get('&b / °').to_numpy()
scan2 = data2.get('R_0 / 1/s').to_numpy()
scan3 = data2.get('R_1 / 1/s').to_numpy()

scan2 = filter(lambda value: isnan(value) == False, scan2)
scan3 = filter(lambda value: isnan(value) == False, scan3)

scan_comb = [*scan2, *scan3]

peaks_index = find_peaks(scan1, prominence=10, height=1, width=(2, 5))[0]
peaks = angles1[peaks_index]

plt1, ax1 = plt.subplots()
ax1.plot(angles1, scan1)
for peak in peaks:
    ax1.axvline(peak, linestyle='dashed', c='r', label=f"200 peak at θ={peak}°")
ax1.set_title("Debye-Scherrer Scan of NaCl Powder")
ax1.set_xlabel("angle θ")
ax1.set_ylabel("Counts per second (10 second average)")
ax1.legend(loc="upper right")
plt1.savefig("DSs_NaCl.png")

print(scan_comb)

plt2, ax2 = plt.subplots()
ax2.plot(angles2, scan_comb)
ax2.set_title("Debye-Scherrer Scan of NaCl Powder for 111 peak")
ax2.set_xlabel("angle θ")
ax2.set_ylabel("Counts per second (10 second average)")
plt2.savefig("DSs_NaCl_11.png")
plt.show()

