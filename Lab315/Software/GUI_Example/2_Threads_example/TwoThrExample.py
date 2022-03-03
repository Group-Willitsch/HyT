# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'TwoThrExample.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(981, 821)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(80, 130, 641, 221))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.lcdNumber_2 = QtWidgets.QLCDNumber(self.verticalLayoutWidget)
        self.lcdNumber_2.setObjectName("lcdNumber_2")
        self.horizontalLayout.addWidget(self.lcdNumber_2)
        self.explanation = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.explanation.setObjectName("explanation")
        self.horizontalLayout.addWidget(self.explanation)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.pushButton = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.pushButton.setObjectName("pushButton")
        self.horizontalLayout_2.addWidget(self.pushButton)
        self.lcdNumber = QtWidgets.QLCDNumber(self.verticalLayoutWidget)
        self.lcdNumber.setObjectName("lcdNumber")
        self.horizontalLayout_2.addWidget(self.lcdNumber)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.killer_button = QtWidgets.QPushButton(self.centralwidget)
        self.killer_button.setGeometry(QtCore.QRect(780, 130, 161, 32))
        self.killer_button.setObjectName("killer_button")
        self.horizontalLayoutWidget_3 = QtWidgets.QWidget(self.centralwidget)
        self.horizontalLayoutWidget_3.setGeometry(QtCore.QRect(130, 410, 791, 321))
        self.horizontalLayoutWidget_3.setObjectName("horizontalLayoutWidget_3")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_3)
        self.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.graphWidget = PlotWidget(self.horizontalLayoutWidget_3)
        self.graphWidget.setObjectName("graphWidget")
        self.horizontalLayout_3.addWidget(self.graphWidget)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.clear_plot_button = QtWidgets.QPushButton(self.horizontalLayoutWidget_3)
        self.clear_plot_button.setObjectName("clear_plot_button")
        self.verticalLayout_2.addWidget(self.clear_plot_button)
        self.checkBox_logscaley = QtWidgets.QCheckBox(self.horizontalLayoutWidget_3)
        self.checkBox_logscaley.setObjectName("checkBox_logscaley")
        self.verticalLayout_2.addWidget(self.checkBox_logscaley)
        self.horizontalLayout_3.addLayout(self.verticalLayout_2)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.horizontalLayout_3.addLayout(self.verticalLayout_3)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 981, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.explanation.setText(_translate("MainWindow", "One random number generated every second."))
        self.pushButton.setText(_translate("MainWindow", "Pass this number to the screen on the right ->"))
        self.killer_button.setText(_translate("MainWindow", "Kill this program"))
        self.clear_plot_button.setText(_translate("MainWindow", "clear plot"))
        self.checkBox_logscaley.setText(_translate("MainWindow", "Enable logscale y"))

from pyqtgraph import PlotWidget
