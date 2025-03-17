
from astropy.io import fits
import numpy as np
import scipy.ndimage as scipy
import pandas as pd
import matplotlib.pyplot as pl

dflux = np.load("/home/rodrigogavioli/PEEC2025/1stPart/SimulationsAreaevolWrot/AreaEvolE030/dflux_inc70.0.npy")
time = dflux[:,0]
flux = dflux[:,1]
pl.plot(time, flux)


mediamovel = scipy.uniform_filter1d(flux, size=50)
flux = flux - mediamovel
dflux[:,1] = flux
pl.plot(time, flux)
np.save("/home/rodrigogavioli/PEEC2025/2ndPart/Wavelet/Dfluxfiltered/DfluxfilteredAreaEvolE03", dflux)

pl.xlabel("Tempo (dias)")
pl.ylabel("Fluxo normalizado")
pl.title("Curva de Luz - Comparação Original vs Filtrada")
pl.legend(["Original", "Filtrada"])
pl.grid()
#pl.show()


#--------------------------------------------------------------------------------#
"""
dflux_data = np.load("/home/rodrigogavioli/PEEC2025/2ndPart/Wavelet/Dfluxfiltered/DfluxfilteredAlpha005.npy")

hdu = fits.PrimaryHDU(dflux_data)
hdul = fits.HDUList([hdu])

fits_file = "/home/rodrigogavioli/meu_ambiente/lib/python3.12/site-packages/star_privateer/timeseries/dfluxfilteredAmplitude005.fits"
hdul.writeto(fits_file, overwrite=True)

print(f"Arquivo FITS salvo em: {fits_file}")

with fits.open(fits_file) as hdul:
    dflux_fits = hdul[0].data  


time_fits = dflux_fits[:, 0]
flux_fits = dflux_fits[:, 1]


pl.plot(time, flux, label="Fluxo Filtrado Original", color="blue", alpha=0.5)
pl.plot(time_fits, flux_fits, label="Fluxo do FITS", color="red", linestyle="dashed")
pl.xlabel("Tempo (dias)")
pl.ylabel("Fluxo Normalizado")
pl.title("Curva de Luz - Comparação Original vs FITS")
pl.legend()
pl.grid()
pl.show()
"""""