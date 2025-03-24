import numpy as np
import pickle 
import pandas as pd
import matplotlib.pyplot as plt
import os

#
"""""
Alpha = -5
for i in range(Alpha,-Alpha+1)
"""""


df = pd.read_csv("/home/rodrigogavioli/PEEC2025/Simulations/table/Alpha-4.8.csv")

IncertezasPPs = df.iloc[:,2]
IncertezasPcs = df.iloc[:,6]
PPs = df.iloc[:,1]
Pcs = df.iloc[:,5]


#Relativeu
RelativeuPPs = IncertezasPPs*100/(PPs)
RelativeuPcs = IncertezasPcs*100/(Pcs)

#standardDeviation
stdPPs = np.std(IncertezasPPs)
stdPcs = np.std(IncertezasPcs)
RelativeustdPPs = np.std(RelativeuPPs)
RelativeustdPcs = np.std(RelativeuPcs)

#Mean
mPPs = np.array([])
mPcs = np.array([])

for i in range(len(PPs)):
    mPPs = np.append(mPPs, np.mean(PPs[:i]))
    mPcs = np.append(mPcs, np.mean(Pcs[:i]))
meanPPs = np.mean(PPs)
meanPcs = np.mean(Pcs)
print(f"meanPPs:{meanPPs}")
print(f"meanPPs:{meanPcs}")
meanRelativeuPPs = np.array([])
meanRelativeuPcs = np.array([])
meanRelativeuPPs = np.append(meanRelativeuPPs, np.mean(RelativeuPPs))
meanRelativeuPcs = np.append(meanRelativeuPcs, np.mean(RelativeuPcs))

#Variation coefficient
CvPPs = RelativeustdPPs/meanRelativeuPPs *100
CvPcs = RelativeustdPcs/meanRelativeuPcs *100
print(f"this is CvPPs{CvPPs}")
print(f"this is CvPcs{CvPcs}")




print(f"StdPPs:{stdPPs}",f"StdPcs:{stdPcs}")
print(f"RelativeustdPPs:{RelativeustdPPs}",f"RelativeuStdPcs:{RelativeustdPcs}")

#Plotting mean for general PPs and Pcs 
plt.plot(mPPs)
plt.plot(mPcs)
plt.legend(["mediaPPs","mediaPcs"])
plt.xlabel("Number of simulations")
plt.grid()
plt.show()
plt.close()

#Plotting for general uPPs and uPcs
plt.plot(IncertezasPPs)
plt.plot(IncertezasPcs)
plt.legend(["uPps","uPcs"])
plt.xlabel("Number of simulations")
plt.grid()
plt.show()
plt.close()

#Plotting relative  uPPs and uPCs
plt.plot(RelativeuPPs)
plt.plot(RelativeuPcs)
plt.legend(["RelativeuPPs","RelativeuPcs"])
plt.xlabel("Number of simulations")
plt.show()
plt.close()

plt.hist(RelativeuPPs)
plt.hist(RelativeuPcs)
plt.legend(["RelativeuPPs","RelativeuPcs"])
plt.xlabel("Number of simulations")
plt.show()

