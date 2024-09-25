# import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv

data = read_csv("Unknown_Scan.csv")
angles = data.get('&b / °').to_numpy()
scan = data.get('R_0 / 1/s').to_numpy()

plt1, ax1 = plt.subplots()
ax1.plot(angles, scan)

plt.show()

