import matplotlib.pyplot as plt
import numpy as np


from main  import source
from scipy import integrate


plt.close('all')




def eulerint_free_flight(pos,vel,dt):
    pos += vel * dt
    return(np.append(pos,vel,axis=0))

source = source() 
particles = source.gauss() 

pos  = particles[:3,:]  ### Position X lab-frame
vel  = particles[3:,:]  ### Position Y lab-frame

print(pos)
dt = 1e-6
bla = eulerint_free_flight(pos,vel,dt)
print(pos)


def derivate(x,t):
    return(x)


X_odeint = integrate.odeint(derivate, pos, dt)



#
#
#
#vx0, vy0 = 1., 1.
#dt = 0.02 
#X0 = np.array([0., 0., vx0, vy0])
#tmax = 0.2
#nt = int(tmax/dt) 
#ti = np.linspace(0., nt * dt, nt)
#
#
#def Euler(func, X0, t):
#  dt = t[1] - t[0]
#  nt = len(t)
#  X  = np.zeros([nt, len(X0)])
#  X[0] = X0
#  for i in range(nt-1):
#    X[i+1] = X[i] + func(X[i], t[i]) * dt
#  return X
#
#
#
#X_euler = Euler(derivate, X0, ti)
#

#
#def derivate(X, t):
#  return np.array([X[2], X[3], 0., -g])
#
#
#X_odeint = integrate.odeint(derivate, X0, ti)
#
#
#
#plt.figure()
#plt.plot(X_odeint[:,0],X_odeint[:,1], "mv", label = "ODEint")
#plt.plot(X_odeint[:,0],X_odeint[:,1], "mv", label = "ODEint")
#
#
#
#plt.grid()
#plt.xlabel("x")
#plt.ylabel("y")
#plt.legend()
#plt.show()
