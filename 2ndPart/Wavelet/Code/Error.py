import numpy as np
import pickle 
import pandas as pd
import matplotlib.pyplot as plt
import os

df = pd.read_csv("/home/rodrigogavioli/PEEC2025/Simulations/table/Alpha-4.9.csv")

IncertezasPPs = df.iloc[:,2]
IncertezasPcs = df.iloc[:,6]



#desvio padrão
stdPPs = np.std(IncertezasPPs)
stdPcs = np.std(IncertezasPcs)
print(f"StdPPs:{stdPPs}",f"StdPcs:{stdPcs}")

plt.plot(IncertezasPPs)
plt.plot(IncertezasPcs)
plt.legend(["ePps","ePcs"])
plt.grid()
plt.show()
plt.close()


