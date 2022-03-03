#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QThread, pyqtSignal
import gui
import sys

class GuiThread(QtWidgets.QMainWindow, gui.Ui_MainWindow):
    workerfunc_gui = pyqtSignal(object)  ###This is a signal used to cummunicate between threads

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)

        self.button.clicked.connect(self.OnClick)
        self.start_working()

    def start_working(self):
        self.workT = WorkThread()

        self.workerfunc_gui.connect(self.workT.workerfunc_work)
        self.workT.Log_work.connect(self.Log)

        self.workT.start()
        return

    def Log(self,msg):
        self.textbox.append(msg)
        return

    def OnClick(self):
        self.somefunction()
        return

    def somefunction(self):
        self.workerfunc_gui.emit('lala')
        return



class WorkThread(QThread):
    Log_work = pyqtSignal(object)
    

    def __init__(self):
        QThread.__init__(self)

    def __del__(self):
        self.wait()

    def run(self):
		#~ start camera
        return

    def workerfunc_work(self,msg):
        self.Log_work.emit(msg)
        return

def main():
    app = QtWidgets.QApplication(sys.argv)
    form = GuiThread()
    form.show()
    app.exec_()

if __name__ == '__main__':
    main()
