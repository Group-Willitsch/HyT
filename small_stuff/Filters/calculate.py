import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import interp1d

plt.close('all')
def bug():
    plt.show()
    sys.exit()

#filters      = ['thorlabs_fguv5-uv','thorlabs_fguv11-uv','thorlabs_fgb37']
filter_names  = ['thorlabs_fgb37',
				'thorlabs_fgb37',
				'thorlabs_fguv5-uv'
				]


###--- Fitting parameters
check_fit = -1    ####-1 => false
inerpolation = 'cubic'  ##'cubic'
wavelength_range = np.linspace(200,1000,1001)

#nh_1 = np.array((305,338))  ####

pump_wl  = 305
probe_wl = 338

plt.figure()
for i in range(len(filter_names)):
    print(filter_names[i])
    try:
        parameters   = np.genfromtxt(str(filter_names[i])+str('.txt'))

    except:
        print('No such filter in list u fool....')
        bug()

    if i == check_fit:
        plt.plot(parameters[:,0],parameters[:,1]/100.,label='Raw')

    ###--- interpolate
    if i == 0:
        f2 = interp1d(parameters[:,0],parameters[:,1]/100., kind=str(inerpolation))
        first = f2(wavelength_range)
        filters = first
        pump  = f2(pump_wl)
        probe = f2(probe_wl)
    else:
        f2_new = interp1d(parameters[:,0],parameters[:,1]/100., kind=str(inerpolation))
        f2 = np.append(f2,f2_new)
        second = f2_new(wavelength_range)

        if i == 1:
            filters = np.vstack((first,second))
        if i >1:
            filters = np.vstack((filters,second))
        pump  = np.vstack((pump,f2_new(pump_wl)))
        probe = np.vstack((probe,f2_new(probe_wl)))


print(np.shape(filters))
try:
    plt.plot(wavelength_range,filters)
except:
    plt.plot(wavelength_range,np.prod(filters,axis=0))

plt.plot((pump_wl,pump_wl),(1e-9,1),label='Pump WL')
plt.plot((probe_wl,probe_wl),(1e-9,1),label='LIF')

plt.ylabel('Transmisson (%)')
plt.xlabel('Wavelength (nm)')
plt.yscale('Log')
plt.legend()
plt.grid()


print('Transmission pump:  '+str(np.prod(pump)))
print('Transmission probe:  '+str(np.prod(probe)))
print('Ratio Probe/Pump:   '+str(np.prod(probe)/np.prod(pump)))




plt.show()
