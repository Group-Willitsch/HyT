import matplotlib.pyplot as plt
from Phidget22.PhidgetException import *
from Phidget22.Phidget import *
from Phidget22.Devices.DigitalInput import *
import traceback
import time

#Declare any event handlers here. These will be called every time the associated event occurs.


def onDigitalInput0_StateChange(self, state):
	time_elapsed.append(time.clock()-t0)
	print(str(time.clock()-t0) +": "+ str(state))

def onDigitalInput0_Attach(self):
	print("Attach!")

def onDigitalInput0_Detach(self):
	print("Detach!")

def onDigitalInput0_Error(self, code, description):
	print("Code: " + ErrorEventCode.getName(code))
	print("Description: " + str(description))
	print("----------")

def main():

	try:
		#Create your Phidget channels
		digitalInput0 = DigitalInput()

		#Set addressing parameters to specify which channel to open (if any)

		#Assign any event handlers you need before calling open so that no events are missed.
		digitalInput0.setOnStateChangeHandler(onDigitalInput0_StateChange)
		digitalInput0.setOnAttachHandler(onDigitalInput0_Attach)
		digitalInput0.setOnDetachHandler(onDigitalInput0_Detach)
		digitalInput0.setOnErrorHandler(onDigitalInput0_Error)

		#Open your Phidgets and wait for attachment
		digitalInput0.openWaitForAttachment(5000)

		#Do stuff with your Phidgets here or in your event handlers.

		try:
			input("Press Enter to Stop\n")
		except (Exception, KeyboardInterrupt):
			pass

		#Close your Phidgets once the program is done.

		digitalInput0.close()
		# print(time_elapsed)
		for t in time_elapsed:
			f.write(str(t)+'\n')

		time_inteval = []
		for i in range(len(time_elapsed)-1):
			time_inteval.append(time_elapsed[i+1]-time_elapsed[i])
			print(time_elapsed[i+1]-time_elapsed[i])

		plt.figure()
		plt.hist(time_inteval,bins=1000)
		plt.show()

	except PhidgetException as ex:
		#We will catch Phidget Exceptions here, and print the error informaiton.
		traceback.print_exc()
		print("")
		print("PhidgetException " + str(ex.code) + " (" + ex.description + "): " + ex.details)


if __name__ == "__main__":
	f = open('test.txt', 'w')
	time_elapsed = []
	t0 = time.clock()
	main()

