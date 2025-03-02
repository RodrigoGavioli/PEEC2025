import numpy as np
import math
import matplotlib.pyplot as pl

# functions for the spot simulation (S15) and the dflux computation (S17)

#-------------------------------------------------------------------------------#
#                              Spot simulation                                  #
#-------------------------------------------------------------------------------#

def nspots(t,nspot):
    # spot number as a function of time - Eq. 1 in S15
	muN=nspot[0]*(t-nspot[1])**3/(np.exp((t-nspot[1])**2/nspot[2]**2)-nspot[3])
	return muN

def lspots(t,cy,lsp):
    # mean latitude as a function of time - Eq. 3 in S15
    # if there is overlap between consecutive cycles, cy indicates each cycle the spots belongs to:
    # in solar-like active bands: previous cycle - low latitudes; new cyle: mid latitudes
	lm=lsp[0]*np.exp(-((t+cy)-lsp[1])/lsp[2])
	return lm

def siglat(t,latm,mn,slp,Pc):
	# width of the spot formation zone as a function of time and also mean latitude - Eq. 4 in S15
	sig=(slp[0]+slp[1]*(t-mn)/Pc+slp[2]*((t-mn)/Pc)**2)*latm
	return sig

def hems(N):
    # hemisphere: Northern (1) or Southern (-1)
	hm=np.random.randint(0,2,N)
	a=np.where(hm==0)
	hm[a]=-1
	return hm 

def meanLpar(t,Pc,mn,oP,lsp,slp):
	# spot parameters at a given phase of the cycle (t)
	# mL - mean latitude
	# sL - sigma: width of the formation zone
	if t<(mn+oP):                                  # overlap of consecutive cycles
		cyu=np.random.randint(0,2)
		if cyu==1:                                 # current cycle
			mL=lspots(t,0.,lsp)
			sL=siglat(t,mL,mn,slp,Pc)
		else:                                      # previous cycle
			mL=lspots(t,Pc-oP,lsp)
			sL=siglat(t+Pc-oP,mL,mn,slp,Pc)
	else:                                          # current cycle
		mL=lspots(t,0.,lsp)
		sL=siglat(t,mL,mn,slp,Pc)
	return mL,sL

def corGW(dgw,cgw,Area_i):
	# lifetime - modified GW rule - Eqs. 2 and 5 in S15
	if Area_i<85:                                  # correction for small spots
		Lf_c=cgw[0]*np.exp(cgw[1]*Area_i)/365.25
	else:
		Lf_c=Area_i/dgw/365.25
	return Lf_c

def wrott(L_i,wr):
    # angular velocity [radians]
    # Sun's profile - Eq. 16 in S15
	wrot=np.radians(wr[0]+wr[1]*np.sin(np.radians(L_i))**2+wr[2]*np.sin(np.radians(L_i))**4)*365.25
	return wrot

def wrott2(L_i,wr2,weq):
	# angular velocity [radians]
	# simplified profile: Eq. 8 in Santos et al. 2017
	wrot=weq*(1.-wr2*np.sin(np.radians(L_i))**2)
	return wrot

def fsevol(evol,A_k):
    # spot growth or decay - power law (4th paragraph in Sect. 2 of S15)
	if A_k <= 0:
		return 0
	dA=np.exp(evol[0])*A_k**evol[1]
	return dA

def spotevol(Amax,Lfd,aevol,cad,dh):
	# spot evolution [days] - meaning it needs to be "normalized" becayse the power law are define based on daily data
	# but for the general case, not trying to reproduce the solar data, we can simplify it
	# See S15: we account for the full decay, then use the remaining time steps for the growth
	At=[Amax]                                        # list of spot areas over the lifetime of the spot
	Ak=Amax
	for ti in np.arange(1,Lfd,dh*365.25):            # decay
		dA=fsevol(aevol,Ak)/cad                      
		if Ak-dA>=0.:
			At.append(Ak-dA)
			Ak-=dA

	tg=Lfd-len(At)/cad                               # remaining steps to be accounted for the growth 
	Ak=Amax
	for ti in np.arange(0,tg,dh*365.25):             # growth
		dA=fsevol(aevol,Ak)/cad
		At.insert(0,Ak-dA)
		Ak-=dA

	return At	

def nasp_count_n(data,tmin,tstep,tf):
	# computes the number of spots & total spot coverage
	# for the visible side
	dh=tstep/24./365.25                           # coverting to years
	time=np.arange(tmin,tf,dh)                    # time array: from vside we only have the spotted times, we need the unspotted times also
	nasp=np.zeros((len(time),3))                  # array where the data is going to be saved
	k=0
	for t in time:
		ind=np.where(abs(data[:,0]-t)<=1e-6)      # localizing the times # tolerance as t slightly changes depending on the computation
		N=len(ind[0])                             # number of spots present at a given t
		A=sum(data[ind[0],2])                     # total spot area
		nasp[k,:]=np.array((t,N,A))
		k+=1
	return nasp

def nasp_count(vside,tmin,tstep,tf):
	# computes the number of spots & total spot coverage
	# for the visible side
	# while more complicated, this code is faster, because each time step I'm reducing the length of the time array to search from
	dh=tstep/24./365.25                           # coverting to years
	time=np.arange(tmin,tf,dh)                    # time array: from vside we only have the spotted times, we need the unspotted times also
	nasp=np.zeros((len(time),3))                  # array where the data is going to be saved
	indx=sorted(range(len(vside[:,0])), key=lambda k: vside[k,0])
	data=vside[indx,:]                            # ordering the table by t instead of spot ID
	k=0
	kk=0
	for t in time:
		N=0
		A=0
		while data[k,0]<=t+1e-6 and k<len(data[:,0])-1: # localizing the times # tolerance as t slightly changes depending on the computation
			if abs(data[k,0]-t)<=1e-6:
				N+=1                              # number of spots present at a given t 
				A+=data[k,2]                      # total spot area
				k+=1
		#print(t,N)
		nasp[kk,:]=np.array((t,N,A))
		kk+=1
	return nasp

def latid(vside):
	# spot latitudes
	ls=max(vside[:,6])                            # array: length - based on the spot ID number
	indl=[0]
	k=0
	for i in range(len(vside[:,0])):
		if vside[i,6]!=vside[indl[k],6]:          # first time the spot appears in the visible side; only taking a latitude per spot
			indl.append(i)
			k+=1
	latt=vside[indl,:]
	print(ls,' generated groups')
	return latt

#-------------------------------------------------------------------------------#
#                         Flux deficit computation                              #
#-------------------------------------------------------------------------------#

def proj_mu(inc,theta,beta):
	# cosine of the angle between the line of sight and the normal to the surface element
	# Eq 6.4 in the Thesis
	pmu=np.cos(np.radians(inc))*np.cos(np.radians(theta))+\
		np.sin(np.radians(inc))*np.sin(np.radians(theta))*np.cos(np.radians(beta))
	return pmu

def limb_dark(Io,mu,limb_D):
	# limd darkening law for solar-like stars
	# Eq. 6 in S17
	Imu=Io*(1.-limb_D[0]*(1.-mu)+limb_D[1]*(1.-mu)**2)
	return Imu


def dflux_s(Io,mu,Rsun,dS,limb_D,Cs):
	# flux deficit due to spot element
	# Eq. 4 in S17
	if 0<=mu<=1:                                              # visible disc
		fp=limb_dark(Io,mu,limb_D)*mu                         # last part of Eq. 4
		fs=fp*(1.-Cs)*dS/np.pi/Rsun**2                        # Eq. 4
	else:
		fs=0                                                  # if in the far side
	return fs

def tbstar(theta_spot,betha_spot,lat_spot):
	# changing coordinates to the star's referential
	# Eq. 6.7 - 6.9 in the Thesis

	t1=np.cos(theta_spot)*np.cos(lat_spot)                    # first term of line1 in Eq. 6.9
	t2=np.sin(theta_spot)*np.cos(betha_spot)*np.sin(lat_spot) # second term of line2 in Eq. 6.9
	theta_star=math.acos(t1+t2)                               # determining thet_star from Eq. 6.9

	sin_theta_star=np.sqrt(1.-(t1+t2)**2)                     # line2 Eq. 6.9

	betha_star=math.asin(np.sin(theta_spot)*np.sin(betha_spot)/sin_theta_star) # line4 Eq. 6.9

	if lat_spot==0:                                           # at the pole - colatitude==0
		betha_star=betha_spot
		tetha_star=theta_spot

	return np.degrees(theta_star),np.degrees(betha_star)

#-------------------------------------------------------------------------------#
#                                     plots                                     #
#-------------------------------------------------------------------------------#

def plot_simul(tmin,tstep,tspan,Nxlim):
	# ploting the summary figure for the spot simulation

	vside=np.loadtxt('vside.txt')
	
	# counts the total number of spots and total spot area
	nasp=nasp_count(vside,tmin,tstep,tspan)
	# spot latitudes
	latt=latid(vside)

	# plot
	pl.subplots(2,2,figsize=(10,6))
	pl.subplots_adjust(hspace=0.3,wspace=0.3)

	# Number of spots ----------------------------------------------------------
	pl.subplot(2,2,1)
	pl.plot(nasp[:,0],nasp[:,1],'r',linewidth=0.7,alpha=0.8)
	pl.ylim(0,max(nasp[:,1])*1.1)
	pl.xlim(tmin,tmin+tspan)
	pl.xlabel(r'$\rm Time\,\, [years]$')
	pl.ylabel(r'$\rm No.\,\, of\,\, Groups$')

	# Total spot area ----------------------------------------------------------
	pl.subplot(2,2,2)
	pl.plot(nasp[:,0],nasp[:,2],'r',linewidth=0.7,alpha=0.8)
	pl.xlim(tmin,tmin+tspan)
	pl.ylim(0,max(nasp[:,2])*1.1)
	pl.ylabel(r'$\rm Total\,\, group\,\, area\,\, [MSH]$')
	pl.xlabel(r'$\rm Time\,\, [years]$')

	# Latitudes ----------------------------------------------------------------
	pl.subplot(2,2,3)
	pl.plot(latt[:,0],latt[:,1],'.r',ms=4)
	pl.yticks(np.arange(-100,100,25))
	pl.ylim(-max(abs(latt[:,1]))-0.5,max(abs(latt[:,1]))+0.5)
	pl.xlim(tmin,tmin+tspan)
	pl.xlabel(r'$\rm Time\,\, [years]$')
	pl.ylabel(r'$\rm Latitude\,\,[\circ]$')

	# Accumulated area distribution --------------------------------------------
	ax=pl.subplot(2,2,4)
	n,b,p=pl.hist(vside[:,2],bins=np.arange(0,max(vside[:,2])*1.02,max(vside[:,2])*0.02),histtype='stepfilled',color='r',alpha=0.75,edgecolor='none')
	ind=np.where(n>=Nxlim)
	pl.xlim(0,b[ind[0][-1]]*1.05)
	pl.ylim(0,max(n)*1.05)
	pl.xlabel(r'$\rm Group\,\, area\,\, [MSH]$')
	pl.ylabel(r'$\rm Spot\,\, count$')
	pl.savefig('simul_plotNA515.png',bbox_inches='tight',dpi=600)

	return nasp

def plot_LC(time,dflux):
	# ploting the LC
	LC=1-dflux

	# plot
	pl.subplots(figsize=(5,2.5))
	pl.plot(time,LC,'k')
	pl.ylabel(r'$\rm Flux$')
	pl.xlabel(r'$\rm Time\,\,[years]$')
	pl.savefig('light_curveNA515.png',bbox_inches='tight',dpi=600)

	return LC
