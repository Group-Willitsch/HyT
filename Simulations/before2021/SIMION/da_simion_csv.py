import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.axes3d as p3
import sys


plt.close('all')
def bug():
    plt.show()
    sys.exit()

data = np.genfromtxt('/Users/thomas/ownCloud/Lab/Data/2020-02-05/data-04-02-2020.csv',delimiter=',')


time_steps = np.unique(data[:,0])

bla = np.empty(len(time_steps-1))
for i in range(len(time_steps)-1):
    bla[i] = time_steps[i+1]-time_steps[i]

###--- Remove lost ions
ion_nbr = np.unique(data[:,1])
counts  = np.zeros(len(ion_nbr))
for i in range(len(ion_nbr)):
    counts[i] = list(data[:,1].flatten()).count(ion_nbr[i])

remove_ions = ion_nbr[counts<np.max(counts)]
for i in range(len(remove_ions)):
    data = data[data[:,1]!=remove_ions[i]]
ion_nbr = np.unique(data[:,1])


####--- Sort for movie
bins = 20
two_d_movie_dat_1 = np.zeros((len(time_steps),bins,bins))
two_d_movie_dat_2 = np.zeros((len(time_steps),bins*2,bins*2))

range_1 = ((2.5,3.5),(2.5,3.5))   ###--- Range for the first histogram
range_2 = ((2.5,3.5),(31.,32.))   ###--- Range for the first histogram

moviedat = np.zeros((len(time_steps),len(ion_nbr),3))
for i in range(len(time_steps)-1):
    moviedat[i,:,0] = data[i*len(ion_nbr):(i+1)*len(ion_nbr),2]
    moviedat[i,:,1] = data[i*len(ion_nbr):(i+1)*len(ion_nbr),3]
    moviedat[i,:,2] = data[i*len(ion_nbr):(i+1)*len(ion_nbr),4]
    two_d_movie_dat_1[i,],x,z = np.histogram2d(moviedat[i,:,0],moviedat[i,:,1],range=range_1,bins=bins)
    two_d_movie_dat_2[i,],x,z = np.histogram2d(moviedat[i,:,1],moviedat[i,:,2],range=range_2,bins=bins*2)


###  100*1000*3   --- 100 Zeitschritte
print('shape of moviedat', moviedat.shape)
gs = plt.GridSpec(1,3)
fig = plt.figure(figsize=(18,5))
ax = fig.add_subplot(gs[0,0], projection='3d')

# Setting the axes properties
ax.set_xlim3d([2.0, 4.0])
ax.set_xlabel('X')
ax.set_ylim3d([2.0, 4.0])
ax.set_ylabel('Y')
ax.set_zlim3d([30.5, 32.5])
ax.set_zlabel('Z')

ax.set_title('3D Test')

ploet = ax.scatter(moviedat[0,:,0], moviedat[0,:,1], moviedat[0,:,2], marker = 'o', alpha = 0.05)
#plt.title('timestamp: %.3f'%(time_steps[0]/1000) + 'ms')

ax2 = fig.add_subplot(gs[0,1])
plot_hist = ax2.imshow(two_d_movie_dat_1[0,:,:],interpolation='none',cmap='hot',extent=(range_1[1][0],range_1[1][1],range_1[0][0],range_1[0][1]))
plt.title('timestamp: %.3f'%(time_steps[0]/1000) + 'ms')
plt.xlabel('X (mm)')
plt.ylabel('Y (mm)')

ax3 = fig.add_subplot(gs[0,2])
plot_hist_fu = ax3.imshow(two_d_movie_dat_2[0,:,:],interpolation='none',cmap='hot',extent=(range_2[1][0],range_2[1][1],range_2[0][0],range_2[0][1]),vmin=0,vmax=15)
plt.title('timestamp: %.3f'%(time_steps[0]/1000) + 'ms')
plt.xlabel('Z (mm)')
plt.ylabel('Y (mm)')

fig.set_tight_layout(True)
#~ lines = [ax.plot(dat[:,0], dat[:,1], dat[:,2])[0] for dat in moviedat]

#~ plt.show()
def update_graph(num):
	data_1 = two_d_movie_dat_1[num,:,:]
	data_2 = two_d_movie_dat_2[num,:,:]

	plt.title('timestamp: %.3f'%(time_steps[num]/1000) + 'ms')
	ploet._offsets3d = (moviedat[num,:,0], moviedat[num,:,1], moviedat[num,:,2])
	plot_hist.set_data(data_1)
	plot_hist_fu.set_data(data_2)

    #return


anim = FuncAnimation(fig, update_graph, frames=np.arange(0, moviedat.shape[0]),interval=5, blit=False)
#anim = FuncAnimation(fig, update_graph, frames=np.arange(0, 20),interval=5, blit=False)

#anim.save('test.gif', writer='imagemagick', fps=20)


plt.figure(5)
gs = plt.GridSpec(1,3)
fig_s = plt.figure(figsize=(18,5))
ax_s = fig_s.add_subplot(gs[0,0])  ##axis static
ax_s = plt.hist(moviedat[-2,:,0],bins=50,histtype='step')

ax_s = fig_s.add_subplot(gs[0,1])  ##axis static
ax_s = plt.hist(moviedat[-2,:,1],bins=50,histtype='step')

ax_s = fig_s.add_subplot(gs[0,2])  ##axis static
ax_s = plt.hist(moviedat[-2,:,2],bins=50,histtype='step')




plt.show()
