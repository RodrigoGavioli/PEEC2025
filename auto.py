import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.ndimage as scipy


filename =f"/home/rodrigogavioli/PEEC2025/Simulations/Wrot/Alpha0.2/dflux_inc_70.0_8.npy"
filename1 =f"/home/rodrigogavioli/PEEC2025/Simulations/Wrot/Alpha0.2/dflux_inc_70.0_7.npy"
arquivo = np.load(filename)
arquivo1 = np.load(filename)
s = arquivo[:,1]
s1 = arquivo1[:,1]
smooth = scipy.uniform_filter1d(s, size=500)
smooth1 = scipy.uniform_filter1d(s1, size=500)
s = s - smooth
s1 = s1 - smooth1
arquivo[:,1] = s
arquivo1[:,1] = s
s = arquivo[:,1]
t = arquivo[:,0]
s1 = arquivo1[:,1]
t1 = arquivo1[:,0]

plt.plot(t1,s1)



plt.plot(t,s)
plt.show()