import math
import numpy as np
from functions_and_countNA515 import *
from main_spot_simulNA515 import *

#-------------------------------------------------------------------------------#
# main program for the Delta_flux computation                                   #
#-------------------------------------------------------------------------------#

def dflux(Pcycle,tstep,tmin,overlap_c,nspot,logn_area,muL,sigmaL,D_GW,C_GW,wrot,area_evol,tspan,cad,Amin,Rsun,inc,Io,Cs,limb_D,nbeta,ntheta):
	# Output files: vside - visible side; bside - both sides - the output is ordered by spot ID number
	tstep=tstep/24./365.25                                      # converting the time step in years
	beta_step=360./(nbeta-1)                                    # step in longitude - spot referential

	# Output files: vside - visible side; bside - both sides - the output is ordered by spot ID number
	# 1. txt as it is lighter than cvs 2. as it has many columns a header is useful, hence not opting for npy
	vsidef=open('vside.txt','w')
	bsidef=open('bside.txt','w')
	vsidef.write('# spot data: visible/near side\n')
	vsidef.write('# time, Latitude, Area, Lifetime, Longitude, Longitude step, spot ID number, maximum Area\n')
	bsidef.write('# spot data: both sides - near and far side\n')
	bsidef.write('# time, Latitude, Area, Lifetime, Longitude, Longitude step, spot ID number, maximum Area\n')

	# main code for the spot-cycle simulation
	bside=simulation(Pcycle,tstep,tmin,overlap_c,nspot,logn_area,muL,sigmaL,D_GW,C_GW,wrot,area_evol,tspan,cad,Amin,vsidef,bsidef)

	vsidef.close()
	bsidef.close()

	# spot data - both sides: near and far side
	bside=np.loadtxt('bside.txt')

	# sorting according to time, instead of spot ID
	indx=sorted(range(len(bside[:,0])), key=lambda k: bside[k,0])
	data=bside[indx,:]

	k=0
	time=np.arange(tmin,tmin+tspan,tstep)                       # time array: from bside we only have the times with spots, we also need the unspotted times

	# Output files: vside - visible side; bside - both sides - the output is ordered by spot ID number
 
	dflux=np.zeros((len(time),2))                               # array for the output
	for t in range(len(time)):
		ind=np.where(abs(data[:,0]-time[t])<=1e-6)              # considering a tolerance as t slightly changes depending on the computation
		N=len(ind[0])                                           # number of spots present at a given day

		Area_s=data[ind[0],2]*1e-6*2.*np.pi*Rsun**2             # spot areas: converting from muHem to cgs units; Hem area -- 1/2 sphere: 2*pi*r**2
		Latitude_s=90-data[ind[0],1]                            # spot co-latitudes
		Longitude_s=np.degrees(data[ind[0],4])-90.              # spot longitudes: center is 0 (in the spot simulation it starts in 0 at the limb)
		#print(time[t],N)

        # spot radius in degrees -> max latitude in the spot's referential
		rs=np.sqrt(Area_s/np.pi)/Rsun
		rd=np.degrees(np.arcsin(rs))

		dfs=0                                                   # flux deficit due to spots
		for n in range(N):	
			theta_step=rd[n]/(ntheta-1)                         # step in latitude

			# cosine of the angle between the line of sight and the normal to the surface element
			# rough estimate of the spot limits to avoid computations for spots in the farside
			mu1=proj_mu(inc,Latitude_s[n],Longitude_s[n]-rd[n]) # mu1: westmost element
			mu2=proj_mu(inc,Latitude_s[n],Longitude_s[n]+rd[n]) # mu2: eastmost element

			# the element is visible when 0<=mu<=1
			if 0.<=mu1<=1. or 0.<=mu2<1.:
				for i in range(ntheta):                         # for each element
					for j in range(nbeta):
						# area of the element k: beta_step*theta_step [degrees] -- coverting unit
						Area_k=Rsun**2*np.sin(np.radians(i*theta_step))*np.radians(theta_step)*np.radians(beta_step)

						# element coordinates in the star's referencial
						theta_st,betha_st=tbstar(np.radians(i*theta_step),np.radians(j*beta_step-90),np.radians(Latitude_s[n]))

						# cosine of the angle between the line of sight and the normal to the surface element
						mu=proj_mu(inc,theta_st,betha_st+Longitude_s[n])

						# flux deficit due to each element
						dfs+=dflux_s(Io,mu,Rsun,Area_k,limb_D,Cs)

		dflux[t,:]=np.array((time[t],dfs))

	# here I chose npy, only two columns and it is lighter than txt or cvs
	np.save('dflux_inc'+str(inc)+'.npy',dflux)
	
	return dflux