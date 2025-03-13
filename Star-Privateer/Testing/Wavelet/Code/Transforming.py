
from astropy.io import fits
import numpy as np
import scipy.ndimage as scipy
import pandas as pd
import matplotlib.pyplot as pl

dflux = np.load("/home/rodrigogavioli/PEEC2025/1stPart/Simulations/AreaChange/Areachangelog3/dflux_inc70.0.npy")
time = dflux[:,0]
flux = dflux[:,1]
pl.plot(time, flux)
pl.show()


mediamovel = scipy.uniform_filter1d(flux, size=50)
flux = flux - mediamovel
print("esta é a media movel do scipy:", mediamovel)
dflux[:,1] = flux
pl.plot(time, flux)
pl.show()
np.save("/home/rodrigogavioli/PEEC2025/Star-Privateer/Dfluxfiltered", dflux)


dflux_data = np.load("/home/rodrigogavioli/PEEC2025/Star-Privateer/Dfluxfiltered/Dfluxfiltered.npy")

hdu = fits.PrimaryHDU(dflux_data)
hdul = fits.HDUList([hdu])

fits_file = "/home/rodrigogavioli/meu_ambiente/lib/python3.12/site-packages/star_privateer/timeseries/dfluxfilteredlog3.fits"
hdul.writeto(fits_file, overwrite=True)

print(f"Arquivo FITS salvo em: {fits_file}")

