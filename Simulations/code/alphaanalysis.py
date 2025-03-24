import pandas as pd
import numpy as np
import matplotlib.pyplot as plt




arrayAlpha = np.array([])
Alpha = -5.0
arraymeanuPPs = np.array([])
arraymeanuPcs = np.array([])
arraymeanPPPs = np.array([])
arraymeanPPcs = np.array([])

for i in range(-5, 6,1):
    i = float(i)
    arrayAlpha = np.append(arrayAlpha,i )
    df = pd.read_csv(f"/home/rodrigogavioli/PEEC2025/Simulations/table/Alpha{i}.csv")

    IncertezasPPs = df.iloc[:,2]
    IncertezasPcs = df.iloc[:,6]
    PPs = df.iloc[:,1]
    Pcs = df.iloc[:,5]
    meanuPPs = np.mean(IncertezasPPs)
    meanuPcs = np.mean(IncertezasPcs)
    arraymeanuPPs = np.append(arraymeanuPPs,meanuPPs)
    arraymeanuPcs = np.append(arraymeanuPcs,meanuPcs)
    PercentagePPs = IncertezasPPs*100/(PPs)
    PercentagePcs = IncertezasPcs*100/(Pcs)

    arraymeanPPPs = np.append(arraymeanPPPs, np.mean(PercentagePPs))
    arraymeanPPcs = np.append(arraymeanPPcs, np.mean(PercentagePcs))


plt.plot(arrayAlpha,arraymeanuPPs)
plt.plot(arrayAlpha,arraymeanuPcs)
plt.ylabel("Mean Absolute Uncertainty (years)")
plt.xlabel("Alpha")
plt.legend(["meanuPPs","meanuPcs"])
plt.grid()
plt.show()

plt.plot(arrayAlpha, arraymeanPPPs)
plt.plot(arrayAlpha, arraymeanPPcs)
plt.ylabel("Mean Relative Uncertainty")
plt.xlabel("Alpha")
plt.legend(["meanuPPPs","meanuPPcs"])
plt.grid()
plt.show()
