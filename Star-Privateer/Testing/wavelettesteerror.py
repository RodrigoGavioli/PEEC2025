import numpy as np
import star_privateer as sp
import matplotlib.pyplot as plt

# Carregar o arquivo de dados (presumivelmente já em anos)
filename = sp.get_target_filename(sp.timeseries, "dflux_inc70.0")
t, s, dt = sp.load_resource(filename)

# Garantir que t está em anos
# Se o arquivo já estiver em anos, não precisamos fazer nada aqui.
# Caso contrário, converta t para anos, por exemplo:
# t = t / 365.25  # Caso o arquivo esteja em dias, converta para anos

# Ajustar o cálculo de dt para anos, agora que t está em anos:
dt = (t[1] - t[0])  # Agora em anos

# Ajustar a resolução de t e s (caso necessário, reduzindo a resolução para médias de 4 pontos)
t = t[:(len(t) // 4) * 4]  # Alterado
t = np.mean(t.reshape(-1, 4), axis=1)
t = t / 365.25  # Converter de dias para anos
s = s[:(len(s) // 4) * 4]  # Alterado
s = np.mean(s.reshape(-1, 4), axis=1)

# Executar o pipeline de análise
(p_wps, p_acf, gwps, wps, acf, cs, coi, features, feature_names, _) = sp.analysis_pipeline(
    t, s, figsize=(8, 12), wavelet_analysis=True, plot=True,
    xlim=(0, 50), normscale='log', ylogscale=True, add_periodogram=True
)

# Calcular o WPS e o período de rotação
dt = (t[1] - t[0]) * 86400  # Ajustar para segundos, se necessário (não essencial para anos)
(periods, wps, gwps, coi, scales) = sp.compute_wps(s, dt, normalise=True, mother=None)

# Ajuste do erro na determinação do período de rotação
(prot_ps, E_prot_ps, param_gauss) = sp.compute_prot_err_gaussian_fit(
    periods, gwps, n_profile=5, threshold=0.1
)

# Plotar o espectro de potência de wavelet
fig = sp.plot_wps(t - t[0], periods, wps, gwps, coi=coi,
                  scales=scales, shading='auto', color_coi='darkgrey',
                  ylogscale=True, lw=1, normscale='log',
                  vmin=None, vmax=None, filename=None, dpi=300,
                  figsize=(8, 4), ylim=(1, 100), show_contour=False,
                  param_gauss=param_gauss)

# AQUI: Garantir que o eixo X é rotulado como "Tempo (anos)"
plt.xlabel("Tempo (anos)")  # Eixo X agora está em anos
plt.show()

# Recalcular e plotar o WPS com o backend 'pywavelets'
dt = (t[1] - t[0]) * 86400  # Ajustar para segundos, se necessário
(periods, wps, gwps, _, scales) = sp.compute_wps(s, dt, normalise=True, mother=None, backend="pywavelets")
(prot_ps, E_prot_ps, param_gauss) = sp.compute_prot_err_gaussian_fit(
    periods, gwps, n_profile=5, threshold=0.1
)

# Plotar o WPS novamente
fig = sp.plot_wps(t - t[0], periods, wps, gwps,
                  scales=scales, shading='auto', color_coi='darkgrey',
                  ylogscale=True, lw=1, normscale='log',
                  vmin=None, vmax=None, filename=None, dpi=300,
                  figsize=(8, 4), ylim=(1, 100), show_contour=False,
                  param_gauss=param_gauss)

# AQUI: Garantir novamente que o eixo X é rotulado como "Tempo (anos)"
plt.xlabel("Tempo (anos)")  # Eixo X agora está em anos
plt.show()
