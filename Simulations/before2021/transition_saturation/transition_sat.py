import numpy as np
import scipy.constants as sc

#####-----#####-----#####-----#####
# The laser pulse duration is < 10 ns, the radiative lifetime of OH
# for the A2Π(v = 1) − X2Π(v = 1) transision is 0.717 ns

###--- Laser ---
## 10. Sep. 2018 we had 4 mW at 10 Hz = 0.4 mW / pulse = 400 uJ
## 19. Nov. 2018 we had 200 uJ

###--- Dominiks folder pinhole measurment
#  0.4*1/e^2 = 0.5 mm (radius)
#  
E_laser   = 200*1e-6 #laser pulse energy, microJoule


wavelength = 282e-9 
e_ph = sc.h*sc.c/wavelength
Nbr_Ph = E_laser/e_ph

dia_laser = 0.001 #laser diameter in meters
A_laser   = dia_laser**2*np.pi/4




###--- Molecule ---
### Absorption cross section from 10.1051/analusis:1999270328
sig_oh = 2.8*1e-16*(1e-2)**2



nbr_pho_per_mol = Nbr_Ph/(A_laser/sig_oh)

print('Nbr of photons per mol: '+str(np.round(nbr_pho_per_mol,2)))
