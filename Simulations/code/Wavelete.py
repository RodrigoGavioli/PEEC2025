import numpy as np 
import star_privateer as sp
import matplotlib.pyplot as plt
import os
#from Saving import *
import pandas as pd
import scipy.ndimage as scipy
import os


Alpha = -4.9

while Alpha <= 5:
    print(Alpha)
    for i in range(20):
        filename =f"/home/rodrigogavioli/PEEC2025/Simulations/Wrot/Alpha{Alpha}/dflux_inc_70.0_{i}.npy"
        arquivo = np.load(filename)
        s = arquivo[:,1]
        smooth = scipy.uniform_filter1d(s, size=500)
        s = s - smooth
        arquivo[:,1] = s
        s = arquivo[:,1]
        t = arquivo[:,0]


        foldergraph = f"/home/rodrigogavioli/PEEC2025/Simulations/graphs/graphAlpha={Alpha}"
        os.makedirs(foldergraph, exist_ok=True)
        foldertable = "/home/rodrigogavioli/PEEC2025/Simulations/table"
        graphname = f"Alpha={Alpha}"
        table_path = f"/home/rodrigogavioli/PEEC2025/Simulations/table/Alpha{Alpha}.csv"



        #dt *= 4
        t = t[:(len(t) // 4) * 4] #Alterado
        t = np.mean(t.reshape(-1, 4), axis=1)
        t = t*365.25
        s = s[:(len(s) // 4) * 4] #Alterado
        s = np.mean (s.reshape (-1,4), axis=1)

        (p_wps, p_acf, gwps, wps, acf,
        cs, coi, features, feature_names, _) = sp.analysis_pipeline (t, s, figsize=(8,12),
                                                                    wavelet_analysis=True, plot=False,
                                                                    xlim=(0,50), normscale='log', ylogscale=True,
                                                                    add_periodogram=True)

        #saving_graphs(folder, graphname)

        dt = (t[1]-t[0])*86400
        (periods, wps, gwps,
        coi, scales) = sp.compute_wps(s, dt, normalise=True, mother=None)

        (prot_ps, E_prot_ps,
        param_gauss) = sp.compute_prot_err_gaussian_fit (periods, gwps, n_profile=5,
                                                            threshold=0.1)
        


        fig = sp.plot_wps(t-t[0], periods, wps, gwps, coi=coi,
                            scales=scales, shading='auto', color_coi='darkgrey',
                            ylogscale=True, lw=1, normscale='log',
                            vmin=None, vmax=None, filename=None, dpi=300,
                            figsize=(8,4), ylim=(1, 100), show_contour=False,
                            param_gauss=param_gauss)
        plt.close(fig)
        #saving_graphs(folder, graphname)



        dt = (t[1] - t[0]) * 86400 

        (periods, wps, gwps,
        _, scales) = sp.compute_wps(s, dt, normalise=True, mother=None,
                                    backend="pywavelets")
        (prot_ps, E_prot_ps,
        param_gauss) = sp.compute_prot_err_gaussian_fit (periods, gwps, n_profile=5,
                                                        threshold=0.1)

        



        fig = sp.plot_wps(t-t[0], periods, wps, gwps,
                        scales=scales, shading='auto', color_coi='darkgrey',
                        ylogscale=True, lw=1, normscale='log',
                        vmin=None, vmax=None, filename=None, dpi=300,
                        figsize=(8,4), ylim=(1, 100), show_contour=False,
                        param_gauss=param_gauss)

        plt.close(fig)
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

       

        df.to_csv(table_path, mode="a", header=not os.path.exists(table_path), index=False)

        print(f"Tabela da simulação {i+1} de 20 salva em ", table_path)
        plt.close()

        #plt.show()
    Alpha += 0.1

    Alpha = round(Alpha,1)
