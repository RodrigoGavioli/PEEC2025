import os
import time
import numpy as np
import matplotlib.pyplot as plt
import tqdm
from main_dfluxAutoSimul import dflux

Nruns = 20
Alpha = -5



while Alpha <= 5:
   for i in tqdm.tqdm(range(Nruns), desc=f"Simulando Alpha={Alpha:.1f}", unit="simulação"): #discovered in starprivateer tqdm and surprisingly works
        timestamp = time.strftime("%d%m%Y-%H%M%S")
        print(timestamp) 
        file = f"/home/rodrigogavioli/PEEC2025/Simulations/Wrot/Alpha{Alpha:.1f}"
        os.makedirs(file, exist_ok=True)
        
        print("Executing simulation:", i+1, "Of", Nruns)
        
        dflux_data = dflux( #all copied from synthetic_LC 
            Pcycle=11.,
            tstep=4.,
            tmin=0,
            overlap_c=11.*0.09,
            nspot=[0.26/10/(1./4*24), 0-1., 11.*0.4, -0.3],
            logn_area=[4., 1.],
            muL=[35, 0, 11.*0.7],
            sigmaL=[0.1, 1.2, -1],
            D_GW=10.,
            C_GW=[5., 6.3e-3],
            wrot=[Alpha, 25.],
            area_evol=[0.17, 0.47],
            tspan=11.,
            cad=1./4*24,
            Amin=1.,
            Rsun=6.955e10,
            inc=70.,
            Io=1.,
            Cs=0.67,
            limb_D=[0.5287, 0.2175],
            nbeta=30,
            ntheta=30
        ) 

        
        np.save(os.path.join(file,f"dflux_inc_70.0_{i}.npy"), dflux_data) #saving

        if os.path.exists("vside.txt"):
            os.rename("vside.txt", os.path.join(file, f"vside_{i}.txt"))
        if os.path.exists("bside.txt"):
            os.rename("bside.txt", os.path.join(file, f"bside_{i}.txt"))

   
   Alpha += 0.1 #start alter the alpha
   Alpha = round(Alpha, 1) #round for naming directorys effects


 