import numpy as np 
import star_privateer as sp
import matplotlib.pyplot as plt
import os
#from Saving import *
import pandas as pd
import scipy.ndimage as scipy




filename ="/home/rodrigogavioli/PEEC2025/Simulations/Wrot/Alpha-5.0"
arquivo = np.load(filename)
t = arquivo[:,1]
filter = scipy.uniform_filter1d(t, size=500)
t = t - filter 
arquivo[:,1] = t
t = arquivo[:,1]
s = arquivo[:,0]

folder = "/home/rodrigogavioli/PEEC2025/2ndPart/Wavelet/Graphs"
foldertable = "/home/rodrigogavioli/PEEC2025/2ndPart/Wavelet/tables"
graphname = "teste2"
 

#dt *= 4
t = t[:(len(t) // 4) * 4] #Alterado
t = np.mean(t.reshape(-1, 4), axis=1)
t = t*365.25
s = s[:(len(s) // 4) * 4] #Alterado
s = np.mean (s.reshape (-1,4), axis=1)

(p_wps, p_acf, gwps, wps, acf,
 cs, coi, features, feature_names, _) = sp.analysis_pipeline (t, s, figsize=(8,12),
                                                             wavelet_analysis=True, plot=True,
                                                             xlim=(0,50), normscale='log', ylogscale=True,
                                                             add_periodogram=True)
for i in range(len(features)):
    print(f"{feature_names[i]}: {features[i]}")

saving_graphs(folder, graphname)

dt = (t[1]-t[0])*86400
(periods, wps, gwps,
 coi, scales) = sp.compute_wps(s, dt, normalise=True, mother=None)

(prot_ps, E_prot_ps,
 param_gauss) = sp.compute_prot_err_gaussian_fit (periods, gwps, n_profile=5,
                                                       threshold=0.1)
print("Erro Gwps", E_prot_ps)


fig = sp.plot_wps(t-t[0], periods, wps, gwps, coi=coi,
                     scales=scales, shading='auto', color_coi='darkgrey',
                     ylogscale=True, lw=1, normscale='log',
                     vmin=None, vmax=None, filename=None, dpi=300,
                     figsize=(8,4), ylim=(1, 100), show_contour=False,
                     param_gauss=param_gauss)

#saving_graphs(folder, graphname)



dt = (t[1] - t[0]) * 86400 

(periods, wps, gwps,
 _, scales) = sp.compute_wps(s, dt, normalise=True, mother=None,
                               backend="pywavelets")
(prot_ps, E_prot_ps,
 param_gauss) = sp.compute_prot_err_gaussian_fit (periods, gwps, n_profile=5,
                                                  threshold=0.1)

print("Erro Gwps2", E_prot_ps)



fig = sp.plot_wps(t-t[0], periods, wps, gwps,
                  scales=scales, shading='auto', color_coi='darkgrey',
                  ylogscale=True, lw=1, normscale='log',
                  vmin=None, vmax=None, filename=None, dpi=300,
                  figsize=(8,4), ylim=(1, 100), show_contour=False,
                  param_gauss=param_gauss)


#saving_graphs(folder, graphname)


Pps = features[0]
Pacf = features[1]
Pcs = features[2]

erroPps = features[3]
erroPacf = features[5]
erroPcs = features[7]

df = pd.DataFrame({
    "Simulação": [graphname],
    "P_PS": [Pps],  
    "Incerteza_P_PS": [erroPps],
    "P_ACF": [Pacf],
    "Incerteza_P_ACF": [erroPacf],
    "P_CS": [Pcs],
    "Incerteza_P_CS": [erroPcs]
    })

print(df)

table_path = "/home/rodrigogavioli/PEEC2025/2ndPart/Wavelet/tables/alpha0.2.csv"

df.to_csv(table_path, mode="a", header=not os.path.exists(table_path), index=False)

print("Tabela salva em ", table_path)

#plt.show()