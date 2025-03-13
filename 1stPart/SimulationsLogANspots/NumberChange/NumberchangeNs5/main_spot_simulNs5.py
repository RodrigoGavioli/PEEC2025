import numpy as np
from functions_and_countNs5 import *

#-------------------------------------------------------------------------------#
# main program for the spot generation                                          #
#-------------------------------------------------------------------------------#

def simulation(Pcycle,tstep,tmin,overlap_c,nspot,logn_area,muL,sigmaL,D_GW,C_GW,wrot,area_evol,tf,cad,Amin,vsidef,bsidef):
	t=tmin
	vside=[]                                                               # visible side
	bside=[]                                                               # both sides
	sid=0                                                                  # spot identification number
	while t<tf:  
		#print(t, tf)
		# At each time t, X new spots are generated, for each we determine its parameters over its lifetime
		mL,sL=meanLpar(t%Pcycle,Pcycle,tmin,overlap_c,muL,sigmaL)          # mean latitude and width of the formatiopn zone at given time t
		mN=nspots(t%Pcycle,nspot)                                          # mean number of spots to be GENERATED at time t (different from total number of spots at time t)
		N=np.random.poisson(mN)                                            # randomly generated spots - Poisson distributon

		if N>0:
			# if spots are generated - random parameters with mean values that depend on t (or change it to be fixed values)
			Longitude=np.random.uniform(0,2*np.pi,N)                     # longitude - uniform distribution:  0-pi visible; pi-2pi far-side
			Latitude=np.random.normal(mL,sL,N)*hems(N)                     # latitude - Gaussian distribution + random hemisphere
			Area=np.random.lognormal(logn_area[0],logn_area[1],N)          # area - lognormal distribution
			sid=sid+np.arange(1,N+1)                                       # spot ID number 

			for i in range(N):
				Life=corGW(D_GW,C_GW,Area[i])                              # lifetime: corrected GW rule
				dw=wrott2(Latitude[i],wrot[0],2*np.pi/wrot[1]*365.25)*tstep   # step in longitude derived from the angular velocity -- Peq being coverted to angular velocity and in years
				Att=spotevol(Area[i],int(Life*365.25)+1,area_evol,cad,tstep)  # spot evolution over its lifetime

				tt=t
				ii=0                                                       # index for the varying spot area Att
				while tt<=t+Life<=tf:                                      # over the spot lifetime
					if Longitude[i]<np.pi and Att[ii]>Amin:                # if in the visible side and A>Amin (visibility threshold)
						vsidef.write(str(tt)+'\t'+str(Latitude[i])+'\t'+str(Att[ii])+'\t'+str(Life)+'\t'+str(Longitude[i])+'\t'+str(dw)+'\t'+str(sid[i])+'\t'+str(Area[i])+'\n')
					if Att[ii]>Amin:                                       # both sides
						bsidef.write(str(tt)+'\t'+str(Latitude[i])+'\t'+str(Att[ii])+'\t'+str(Life)+'\t'+str(Longitude[i])+'\t'+str(dw)+'\t'+str(sid[i])+'\t'+str(Area[i])+'\n')
					tt+=tstep 
					ii+=1
					Longitude[i]+=dw
					if Longitude[i]>2.*np.pi:
						Longitude[i]-=2.*np.pi                             # start in 0 again
			sid=sid[N-1]
		t+=tstep
	return N