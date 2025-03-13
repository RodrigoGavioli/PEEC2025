import numpy as np 
import star_privateer as sp
import matplotlib.pyplot as plt
import os
filename = sp.get_target_filename (sp.timeseries,"003733735")
t, s, dt = sp.load_resource (filename)

dt *= 4
t = t[:(len(t) // 4) * 4] #Alterado
t = np.mean(t.reshape(-1, 4), axis=1)
#t = t*365.25
s = s[:(len(s) // 4) * 4] #Alterado
s = np.mean (s.reshape (-1,4), axis=1)

(p_wps, p_acf, gwps, wps, acf,
 cs, coi, features, feature_names, _) = sp.analysis_pipeline (t, s, figsize=(8,12),
                                                             wavelet_analysis=True, plot=True,
                                                             xlim=(0,50), normscale='log', ylogscale=True,
                                                             add_periodogram=True)
for i in range(len(features)):
    print(f"{feature_names[i]}: {features[i]}")


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



graph_name = f"teste{len(os.listdir("/home/rodrigogavioli/PEEC2025/Star-Privateer/Testing/Wavelet/Graphs")) + 1}.png"
output_path = os.path.join("/home/rodrigogavioli/PEEC2025/Star-Privateer/Testing/Wavelet/Graphs", graph_name)


#plt.savefig(output_path, dpi=300, bbox_inches='tight')

plt.show()
plt.close("all")


