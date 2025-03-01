#--------------------------------------------------------------------------------#
# A.R.G. Santos                                                                  #
#                                                                                #
# Synthetic light curves of spotted stars, combining the tools for spot          #
# simulations and for the flux deficit due to spots.                             #
#                                                                                #
# Empirical tool to simulate a spot cycle. The input parameters are described in #
# Santos et al. 2015,A&A,580,A62 (hereafter, S15).                               #
#                                                                                #
# The functional forms for the diferent cycle peroperties are based on the       #
# properties of the solar cycle.                                                 #
# The parameters have been changed to more general cases. If the parameters of   #
# SC23 or SC22 are needed, ask A.R.G. Santos.                                    #
#                                                                                #
# Tool for the flux decrease computation is based on Santos et al. 2017,         #
# A&A,599,A1 (hereafter S17) and the PhD Thesis.                                 #
#--------------------------------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as pl
from main_dfluxNs1 import *
from functions_and_countNs1 import *

#--------------------------------parameters--------------------------------------#
plot_vis=1                                # plot visualization: if 1, the summary figure is produced - visible side
Nxlim=100.                                # minimum count to set the xlim of the area distribution

Pcycle=11.                                # cycle period [years]
tstep=4.                                  # time step [hours]
tmin=0                                    # initial time [years]
tspan=tmin+Pcycle                         # final time - full cycle or less [years]
overlap_c=Pcycle*0.09                     # overlap of consecutive cycles [years]
Amin=1.                                   # minimum spot area [muHem]

# defining a normalizing constant to adjust the parameters of the solar functional forms
# this is to be consistent with the sunspot cycle
cad=1./tstep*24.

# SPOT NUMBER: assymetric "Gaussian" (Eq. 1 in S15)
# in order: amplitude; starting time - earlier than tmin; scale related to Pcycle; assymetry
nspot=[0.1/10/cad,tmin-1.,Pcycle*0.4,-0.3]          

# MAXIMUM AREA: lognormal distribution
# in order: mu; sigma
logn_area=[4.,1.]

# SPOT EVOLUTION - spot varying size: power law (4th paragraph in Sect. 2 of S15)
# in order: amplitude; exponent
# for simplification: equal growth and decay rates
area_evol=[0.17,0.47]

# MEAN LATITUDE: exponential (Eq. 3 in S15)
# in order: imitial latitude; initial time; scale related to Pcycle
muL=[35,tmin,Pcycle*0.7]

# WIDTH - LATITUDE: 2nd order polynomial (Eq. 4 in S15)
# in order: coefficientsof the polynomial
sigmaL=[0.1,1.2,-1]

# LIFETIME: GW rule with a correction for the small spots (Eqs. 2 and 5 in S15)
D_GW=10.                                                      # GW constant
C_GW=[5.,6.3e-3]                                              # correction to the GW rule
# ANGULAR VELOCITY: rotation profile
# For the Sun, consider a more accurate relation: Eq. 16 in S15. But it is being generalized
# to the Eq. 8 in Santos et al. 2017.
# in order: alpha (shear); P_eq (rotation period at the equator in days) 
wrot=[0.2,25.] 

# parameters for the flux computation
Rsun=6.955e10                             # solar radius [cgs]
inc=70.                                   # inclination angle
Io=1.									  # intesity at the center of the disc
Cs=0.67                                   # spot intensity contrast
limb_D=[0.5287,0.2175]                    # limb darkening -- Eq. 6 in S17 -- solar-like values

# slicing the spot in surface elements
# in the spot's referential -- latitude varying from 0-R_spot -- longitude varying from 0-360
# Fig. 6.2 in the Thesis and nearby equations
nbeta=30                                  # number of elements in longitude 
ntheta=30                                 # number of elements in latitude
# This choise of elements leads to an overestimating of the flux deficit of about 5% in relation
# to the reference (Notebook around 14/04/2022), which I considered to be reasonable. Note that
# while comparing artificial data, this parameters MUST be unchanged as they change the relative
# flux decrease.

#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

# main code for the light curve generation
dflux=dflux(Pcycle,tstep,tmin,overlap_c,nspot,logn_area,muL,sigmaL,D_GW,C_GW,wrot,area_evol,\
	tspan,cad,Amin,Rsun,inc,Io,Cs,limb_D,nbeta,ntheta)

#-------------------------------------------------------------------------------#
# plot visualization
if plot_vis==1:
	nasp=plot_simul(tmin,tstep,tspan,Nxlim)

	Lc=plot_LC(dflux[:,0],dflux[:,1])

