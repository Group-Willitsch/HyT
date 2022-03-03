import numpy as np
import sys
from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import QThread, pyqtSignal, Qt

from Mainwindow import Ui_MainWindow
import locale
import time
from xmlrpc.client import ServerProxy
import numpy as np

class GuiThread(QtWidgets.QMainWindow, Ui_MainWindow):
    ###pyqtSignals are used to cummunicate between threads, here GuiThread and testprint
    stopper_sig = pyqtSignal()         ###Stopping command, no object (or data) required
    set_WL_sig  = pyqtSignal(object)   ###Send the wavelength (WL, aka here just a number) from one thread to another

    def __init__(self):
        super(self.__class__, self).__init__()       ###  Some standard lines to set up the gui
        self.setupUi(self)

        ###Variables
        self.show_nbr = 0                             ### Variable to display the number
        self.lcdNumber.display(str(self.show_nbr))    ### Connect variable to display
        self.set_WL = locale.atof(self.Set_wavelength.text())   ###Read the value. Had a problem with the "," in german language, that is why locale.atof

        ###Buttons actions
        self.start_button.clicked.connect(self.StartPIDen)    ####Button names are defined in mainwindow.py
        self.stop_button.clicked.connect(self.EndPIDen)       ####Button names are defined in mainwindow.py

        ###Threads
        self.printer = testprint()                            ####This is the second thread

        ###Connect threads with methods (fcts)
        self.printer.wl_printer.connect(self.update_display_wl)  ###Connect from "emitter -> reciever"
        self.stopper_sig.connect(self.printer.stop_jump)         ###Connect from "emitter -> reciever"
        self.Set_wavelength.valueChanged.connect(self.update_wl)


    def StartPIDen(self):                                      ###Start executing
        self.printer.start()                                   #### .start() executes ALWAYS THE "run():" method !!!!!!!!!!!!!!!
    def EndPIDen(self):
        self.stopper_sig.emit()                                ###Emit the stopping signa

    def update_display_wl(self,current_wl):                    ###Update the display of the ui
        self.show_nbr = current_wl
        self.lcdNumber.display(str(self.show_nbr))

    def update_wl(self):
        self.set_WL = locale.atof(self.Set_wavelength.text())
        print(str(self.set_WL))
        return


class testprint(QThread):                                       ###QThread is used to initilize this as a Thread
    wl_printer = pyqtSignal(object)                             ###Used to send the current wl to the GUI Thread

    def __init__(self):
        QThread.__init__(self)
        self.printing = True
        self.counter = 0

    def run(self):                                              #####--- Put here the method to be executed in this new thread
        self.printing = True
        while self.printing == True:
            print(str(self.counter))
            self.wl_printer.emit(self.counter)
            self.counter+=1
            time.sleep(1.0)

    def stop_jump(self):
        self.printing = False
        self.counter = 0



def main():
    app = QtWidgets.QApplication(sys.argv)
    window = GuiThread()
    window.show()
    app.exec()


if __name__ == '__main__':
    main()
