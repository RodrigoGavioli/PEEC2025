import math
import numpy as np
from functions_and_count import *
from main_spot_simul import *

def dflux(Pcycle,tstep,tmin,overlap_c,nspot,logn_area,muL,sigmaL,D_GW,C_GW,wrot,area_evol,tspan,cad,Amin,Rsun,inc,Io,Cs,limb_D,nbeta,ntheta):
	tstep = tstep / 24. / 365.25  
	beta_step = 360. / (nbeta - 1)  


	vsidef = open('vside.txt', 'w')
	bsidef = open('bside.txt', 'w')
	vsidef.write('# spot data: visible/near side\n')
	vsidef.write('# time, Latitude, Area, Lifetime, Longitude, Longitude step, spot ID number, maximum Area\n')
	bsidef.write('# spot data: both sides - near and far side\n')
	bsidef.write('# time, Latitude, Area, Lifetime, Longitude, Longitude step, spot ID number, maximum Area\n')

	
	bside = simulation(Pcycle, tstep, tmin, overlap_c, nspot, logn_area, muL, sigmaL, D_GW, C_GW, wrot, area_evol, tspan, cad, Amin, vsidef, bsidef)

	vsidef.close()
	bsidef.close()

	
	bside = np.loadtxt('bside.txt')

	
	indx = sorted(range(len(bside[:, 0])), key=lambda k: bside[k, 0])
	data = bside[indx, :]

	time = np.arange(tmin, tmin + tspan, tstep)  

	#dflux = np.zeros((len(time), 2))
	dflux_array = np.zeros((len(time), 2))  # Altered for autocode, so the function was not an array

	for t in range(len(time)):
		ind = np.where(abs(data[:, 0] - time[t]) <= 1e-6)  
		N = len(ind[0]) 

		Area_s = data[ind[0], 2] * 1e-6 * 2. * np.pi * Rsun**2 
		Latitude_s = 90 - data[ind[0], 1]  
		Longitude_s = np.degrees(data[ind[0], 4]) - 90. 

		rs = np.sqrt(Area_s / np.pi) / Rsun
		rd = np.degrees(np.arcsin(rs))

		dfs = 0  
		for n in range(N):
			theta_step = rd[n] / (ntheta - 1)

			mu1 = proj_mu(inc, Latitude_s[n], Longitude_s[n] - rd[n])  
			mu2 = proj_mu(inc, Latitude_s[n], Longitude_s[n] + rd[n])

			if 0. <= mu1 <= 1. or 0. <= mu2 < 1.:
				for i in range(ntheta):
					for j in range(nbeta):
						Area_k = Rsun**2 * np.sin(np.radians(i * theta_step)) * np.radians(theta_step) * np.radians(beta_step)

						theta_st, betha_st = tbstar(np.radians(i * theta_step), np.radians(j * beta_step - 90), np.radians(Latitude_s[n]))

						mu = proj_mu(inc, theta_st, betha_st + Longitude_s[n])

						dfs += dflux_s(Io, mu, Rsun, Area_k, limb_D, Cs)
		#dflux[t, :] = np.array((time[t], dfs))
		dflux_array[t, :] = np.array((time[t], dfs)) #Altered to new variable

	#np.save(f'dflux_inc_{inc}.npy', dflux)
	#np.save(f'dflux_inc_{inc}.npy', dflux_array) # Altered to return the new variable 

	#return dflux
	return dflux_array  # Altered to return the new variable 
