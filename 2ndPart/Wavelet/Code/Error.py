import numpy as np
import pickle 
import pandas as pd
import matplotlib.pyplot as plt
import os

df = pd.read_csv("/home/rodrigogavioli/PEEC2025/2ndPart/Wavelet/tables/alpha0.2.csv")

IncertezasPPs = df.iloc[:,2]
IncertezasPcs = df.iloc[:,6]

print(IncertezasPcs)
plt.plot(IncertezasPPs)
plt.plot(IncertezasPcs)
plt.legend(["ePps","ePcs"])
plt.grid()
plt.show()
plt.close()

#desvio padrão
stdPPs = np.std(IncertezasPPs)
stdPcs = np.std(IncertezasPcs)
print(f"StdPPs:{stdPPs}",f"StdPPs:{stdPcs}")


