from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        self.title = 'PyQt5 textbox - pythonspot.com'
        self.left = 10
        self.top = 10
        self.width = 400
        self.height = 140
        #~ self.initUI()

        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
 
        # Create textbox
        self.textbox = QtWidgets.QTextBrowser(self)
        self.textbox.move(20, 20)
        self.textbox.resize(280,40)
 
        # Create a button in the window
        self.button = QtWidgets.QPushButton('Show text', self)
        self.button.move(20,80)
 
        # connect button to function on_click
