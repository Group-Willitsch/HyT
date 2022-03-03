'''
	plotsequence
'''
import numpy as np
import matplotlib.pyplot as plt


def plotsequence(output,checkboxvals):
	'''
		takes output from pulses and plots the pulsesequences
	'''
	xdata=[]
	ydata=[]

	for i in range(len(output[0])):
		y=[]
		xdata.append(float(output[1][i]))
		for j in range(len(output[0][i])):
			y.append(float(output[0][i][j])+2*j)
		ydata.append(y)
		
		
	ydatanp=np.array([[]])
	for y in ydata:
		ydatanp=np.append(ydatanp,y)
		
	xdatanp=np.array(output[1]).astype(float)
	xdatanp=np.cumsum(xdatanp)
	xdatanp=np.append(xdatanp,[0])
	xdatanp=np.sort(xdatanp)
	xdatanp=np.delete(xdatanp,xdatanp.shape[0]-1)
	ydatanp=np.reshape(ydatanp,(len(ydata),len(ydata[0]))).T
	
	
	steps=10000
	stepsize=1e-7
	dt=1e-7
	xcont=[]
	ycont=[]
	
	newlabels=[]
	no=0
	for channel in checkboxvals:
		y=[]
		x=[]
		time=0
		event=0
		if channel:
			newlabels.append(no+1)
			
				
			while True:
				time=xdatanp[event]
				x.append(time)
				y.append(ydatanp[no][event])
					
				event+=1
				try:
					time=xdatanp[event] - dt
					x.append(time)
					y.append(ydatanp[no][event-1])
				
					time=xdatanp[event]
				except IndexError:
					#~ print('indexerr')
					break
					
			x.append(xdatanp[-1])
			y.append(ydatanp[no][-1])

				
				
			ycont.append(y)
			xcont.append(x)
		no+=1
	
	return xcont, ycont, newlabels
	
	

