

'''
Stupid example to learn how to use 2 threads at the same time.
One thread (gui, the main one) will generate one random number every second.
The second thread will check continuosly if u press a print button. If pressed, it will print the number on a second scree.
Author: PV, October 2020
'''


from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QThread, pyqtSignal
from TwoThrExample import Ui_MainWindow
import pyqtgraph as pg
import sys
import time
import numpy as np
#from numba import jit #try numba via decorator @jit(nopython=True) for improved speed with numpy functions.



class GuiThread(QtWidgets.QMainWindow, Ui_MainWindow):
#    get_rnd_num = pyqtSignal(object)

    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)    #mumbo jumbo to inizialize the class
        Ui_MainWindow.__init__(self)            #same
        self.setupUi(self)                      #same

        #variables
        self.current_rnd_num = 0.1

        #worker thread
        self.worker = Worker()  #init
        self.worker.start()     #mumbo jumbo^2: start always calls the run() method, which will call the exec()

        #screen init (start always with current random number)
        self.lcdNumber_2.display(self.current_rnd_num)
        #Graph init
        self.graphWidget

        #connections
        self.worker.send_rnd_num.connect(self.update_display_num)
        self.pushButton.clicked.connect(self.push_button_has_been_clicked)
        self.killer_button.clicked.connect(self.killer_button_has_been_clicked)
        #end of constructor

    def update_display_num(self, new_rnd_num_guiThread):
        #self.current_rnd_num = self.worker.rnd_num_generator
        self.current_rnd_num = new_rnd_num_guiThread
        self.lcdNumber_2.display(self.current_rnd_num)
        return

    def push_button_has_been_clicked(self):
        print("push button clicked")
        print(self.current_rnd_num)
        self.lcdNumber.display(self.current_rnd_num)
        return

    def killer_button_has_been_clicked(self):
        print("Killer button pressed\nQuitting the program...\nCiao frate'\r\n")
        sys.exit() #may be dirty (leave threads occupied??).
        #The exit code is 0 tough, u can see it on win via echo  %ErrorLevel%.  





class Worker(QThread):
    '''
    This is the thread that makes calculations:
    1) generates one rnd number every second
    2)
    '''
    send_rnd_num = pyqtSignal(object) #to send the generated number to the GuiThread
    #NOTE: pyqtSignal HAS TO BE OUTSIDE THE CONSTRUCTOR, FOR SOME UNKNOWN REASON. It's related to difference between bound/unbound signals.

    def __init__(self):
        QThread.__init__(self) #mumbo jumbo to call the init of QThread
        self.new_rnd_num = 0.
        self.waiting_time_rnd = .5 #waiting time to be passed to time.sleep, in seconds

    #QThreads begin executing in run(). By default, run() starts the event loop by calling exec() and runs a Qt event loop inside the thread.
    #REF: https://www.riverbankcomputing.com/static/Docs/PyQt5/api/qtcore/qthread.html?highlight=qthread
    def run(self):
        while True:
            time.sleep(self.waiting_time_rnd)
            self.new_rnd_num = np.random.normal()
#            print(self.new_rnd_num)
            self.send_rnd_num.emit(self.new_rnd_num)


















def main():
    app = QtWidgets.QApplication(sys.argv)
    window = GuiThread()
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
