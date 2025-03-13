from scipy.ndimage import uniform_filter1d
import numpy as np
import pandas as pd

a = np.load("/home/rodrigogavioli/PEEC2025/1stPart/Simulations/AreaChange/Areachangelog3/dflux_inc70.0.npy")
print(a)
time = a[:,0]
flux = a[:,1]


mediamovel = uniform_filter1d(flux, size=50)
flux = flux - mediamovel

a[:,1] = flux
np.save("/home/rodrigogavioli/PEEC2025/Star-Privateer/Dfluxfiltered", a)




# print(uniform_filter1d([2, 8, 0, 4, 1, 9, 9, 0], size=3))