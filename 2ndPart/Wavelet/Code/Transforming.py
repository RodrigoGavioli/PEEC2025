
from astropy.io import fits
import numpy as np
import scipy.ndimage as scipy
import pandas as pd
import matplotlib.pyplot as pl

dflux = np.load("/home/rodrigogavioli/PEEC2025/1stPart/SimulationsLogANspots/Initial_Data/dflux_inc70.0.npy")
time = dflux[:,0]
flux = dflux[:,1]
pl.plot(time, flux)


mediamovel = scipy.uniform_filter1d(flux, size=50)
flux = flux - mediamovel
dflux[:,1] = flux
pl.plot(time, flux)
np.save("/home/rodrigogavioli/PEEC2025/2ndPart/Wavelet/Dfluxfiltered", dflux)


dflux_data = np.load("/home/rodrigogavioli/PEEC2025/2ndPart/Wavelet/Dfluxfiltered/Dfluxfiltered.npy")

hdu = fits.PrimaryHDU(dflux_data)
hdul = fits.HDUList([hdu])

fits_file = "/home/rodrigogavioli/meu_ambiente/lib/python3.12/site-packages/star_privateer/timeseries/dfluxfilteredInitialData.fits"
hdul.writeto(fits_file, overwrite=True)

print(f"Arquivo FITS salvo em: {fits_file}")

