# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'under_pressure.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_my_pressures(object):
    def setupUi(self, my_pressures):
        my_pressures.setObjectName("my_pressures")
        my_pressures.resize(924, 735)
        my_pressures.setIconSize(QtCore.QSize(100, 100))
        self.centralWidget = QtWidgets.QWidget(my_pressures)
        self.centralWidget.setObjectName("centralWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralWidget)
        self.horizontalLayout.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout.setSpacing(6)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setSpacing(6)
        self.verticalLayout.setObjectName("verticalLayout")
        self.frame_2 = QtWidgets.QFrame(self.centralWidget)
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_2.setObjectName("frame_2")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.frame_2)
        self.horizontalLayout_2.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout_2.setSpacing(6)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout()
        self.verticalLayout_4.setSpacing(6)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.frame_9 = QtWidgets.QFrame(self.frame_2)
        self.frame_9.setMinimumSize(QtCore.QSize(200, 0))
        self.frame_9.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_9.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_9.setObjectName("frame_9")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.frame_9)
        self.verticalLayout_6.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_6.setSpacing(6)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.frame_15 = QtWidgets.QFrame(self.frame_9)
        self.frame_15.setMinimumSize(QtCore.QSize(250, 0))
        self.frame_15.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_15.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_15.setObjectName("frame_15")
        self.label_16 = QtWidgets.QLabel(self.frame_15)
        self.label_16.setGeometry(QtCore.QRect(10, 10, 81, 16))
        self.label_16.setObjectName("label_16")
        self.COMbox = QtWidgets.QComboBox(self.frame_15)
        self.COMbox.setGeometry(QtCore.QRect(10, 30, 69, 22))
        self.COMbox.setObjectName("COMbox")
        self.ClearButton = QtWidgets.QPushButton(self.frame_15)
        self.ClearButton.setGeometry(QtCore.QRect(10, 100, 75, 23))
        self.ClearButton.setObjectName("ClearButton")
        self.label_17 = QtWidgets.QLabel(self.frame_15)
        self.label_17.setGeometry(QtCore.QRect(10, 130, 47, 13))
        self.label_17.setObjectName("label_17")
        self.PlotBox = QtWidgets.QComboBox(self.frame_15)
        self.PlotBox.setGeometry(QtCore.QRect(10, 150, 69, 22))
        self.PlotBox.setObjectName("PlotBox")
        self.PlotBox.addItem("")
        self.PlotBox.addItem("")
        self.label_18 = QtWidgets.QLabel(self.frame_15)
        self.label_18.setGeometry(QtCore.QRect(10, 180, 161, 16))
        self.label_18.setObjectName("label_18")
        self.twindowTB = QtWidgets.QTextBrowser(self.frame_15)
        self.twindowTB.setGeometry(QtCore.QRect(10, 200, 81, 31))
        self.twindowTB.setObjectName("twindowTB")
        self.tstepTB = QtWidgets.QTextEdit(self.frame_15)
        self.tstepTB.setGeometry(QtCore.QRect(10, 250, 71, 31))
        self.tstepTB.setObjectName("tstepTB")
        self.label_19 = QtWidgets.QLabel(self.frame_15)
        self.label_19.setGeometry(QtCore.QRect(10, 230, 101, 16))
        self.label_19.setObjectName("label_19")
        self.LogScaleCB = QtWidgets.QCheckBox(self.frame_15)
        self.LogScaleCB.setGeometry(QtCore.QRect(10, 300, 70, 17))
        self.LogScaleCB.setObjectName("LogScaleCB")
        self.COMButton = QtWidgets.QPushButton(self.frame_15)
        self.COMButton.setGeometry(QtCore.QRect(90, 60, 75, 23))
        self.COMButton.setObjectName("COMButton")
        self.ConCOMButton = QtWidgets.QPushButton(self.frame_15)
        self.ConCOMButton.setGeometry(QtCore.QRect(10, 60, 75, 23))
        self.ConCOMButton.setObjectName("ConCOMButton")
        self.label_20 = QtWidgets.QLabel(self.frame_15)
        self.label_20.setGeometry(QtCore.QRect(90, 10, 91, 16))
        self.label_20.setObjectName("label_20")
        self.COMbox2 = QtWidgets.QComboBox(self.frame_15)
        self.COMbox2.setGeometry(QtCore.QRect(90, 30, 69, 22))
        self.COMbox2.setObjectName("COMbox2")
        self.pic_frame = QtWidgets.QFrame(self.frame_15)
        self.pic_frame.setGeometry(QtCore.QRect(0, 320, 241, 221))
        self.pic_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.pic_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.pic_frame.setObjectName("pic_frame")
        self.verticalLayout_13 = QtWidgets.QVBoxLayout(self.pic_frame)
        self.verticalLayout_13.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_13.setSpacing(6)
        self.verticalLayout_13.setObjectName("verticalLayout_13")
        self.pic_frame_layout = QtWidgets.QVBoxLayout()
        self.pic_frame_layout.setSpacing(6)
        self.pic_frame_layout.setObjectName("pic_frame_layout")
        self.verticalLayout_13.addLayout(self.pic_frame_layout)
        self.verticalLayout_6.addWidget(self.frame_15)
        self.verticalLayout_4.addWidget(self.frame_9)
        self.horizontalLayout_2.addLayout(self.verticalLayout_4)
        self.verticalLayout_5 = QtWidgets.QVBoxLayout()
        self.verticalLayout_5.setSpacing(6)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.frame_8 = QtWidgets.QFrame(self.frame_2)
        self.frame_8.setMinimumSize(QtCore.QSize(200, 0))
        self.frame_8.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_8.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_8.setObjectName("frame_8")
        self.Ch1TB = QtWidgets.QTextBrowser(self.frame_8)
        self.Ch1TB.setGeometry(QtCore.QRect(10, 30, 71, 31))
        self.Ch1TB.setObjectName("Ch1TB")
        self.label = QtWidgets.QLabel(self.frame_8)
        self.label.setGeometry(QtCore.QRect(10, 10, 71, 16))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.frame_8)
        self.label_2.setGeometry(QtCore.QRect(90, 40, 47, 13))
        self.label_2.setObjectName("label_2")
        self.verticalLayout_5.addWidget(self.frame_8)
        self.frame_6 = QtWidgets.QFrame(self.frame_2)
        self.frame_6.setMinimumSize(QtCore.QSize(200, 0))
        self.frame_6.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_6.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_6.setObjectName("frame_6")
        self.label_3 = QtWidgets.QLabel(self.frame_6)
        self.label_3.setGeometry(QtCore.QRect(90, 40, 47, 13))
        self.label_3.setObjectName("label_3")
        self.Ch2TB = QtWidgets.QTextBrowser(self.frame_6)
        self.Ch2TB.setGeometry(QtCore.QRect(10, 30, 71, 31))
        self.Ch2TB.setObjectName("Ch2TB")
        self.label_4 = QtWidgets.QLabel(self.frame_6)
        self.label_4.setGeometry(QtCore.QRect(10, 10, 71, 16))
        self.label_4.setObjectName("label_4")
        self.verticalLayout_5.addWidget(self.frame_6)
        self.frame_11 = QtWidgets.QFrame(self.frame_2)
        self.frame_11.setMinimumSize(QtCore.QSize(200, 0))
        self.frame_11.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_11.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_11.setObjectName("frame_11")
        self.label_5 = QtWidgets.QLabel(self.frame_11)
        self.label_5.setGeometry(QtCore.QRect(90, 40, 47, 13))
        self.label_5.setObjectName("label_5")
        self.Ch3TB = QtWidgets.QTextBrowser(self.frame_11)
        self.Ch3TB.setGeometry(QtCore.QRect(10, 30, 71, 31))
        self.Ch3TB.setObjectName("Ch3TB")
        self.label_6 = QtWidgets.QLabel(self.frame_11)
        self.label_6.setGeometry(QtCore.QRect(10, 10, 71, 16))
        self.label_6.setObjectName("label_6")
        self.verticalLayout_5.addWidget(self.frame_11)
        self.frame_10 = QtWidgets.QFrame(self.frame_2)
        self.frame_10.setMinimumSize(QtCore.QSize(200, 0))
        self.frame_10.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_10.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_10.setObjectName("frame_10")
        self.label_7 = QtWidgets.QLabel(self.frame_10)
        self.label_7.setGeometry(QtCore.QRect(90, 40, 47, 13))
        self.label_7.setObjectName("label_7")
        self.Ch4TB = QtWidgets.QTextBrowser(self.frame_10)
        self.Ch4TB.setGeometry(QtCore.QRect(10, 30, 71, 31))
        self.Ch4TB.setObjectName("Ch4TB")
        self.label_8 = QtWidgets.QLabel(self.frame_10)
        self.label_8.setGeometry(QtCore.QRect(10, 10, 91, 16))
        self.label_8.setObjectName("label_8")
        self.verticalLayout_5.addWidget(self.frame_10)
        self.frame_7 = QtWidgets.QFrame(self.frame_2)
        self.frame_7.setMinimumSize(QtCore.QSize(200, 0))
        self.frame_7.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_7.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_7.setObjectName("frame_7")
        self.label_9 = QtWidgets.QLabel(self.frame_7)
        self.label_9.setGeometry(QtCore.QRect(90, 40, 47, 13))
        self.label_9.setObjectName("label_9")
        self.Ch5TB = QtWidgets.QTextBrowser(self.frame_7)
        self.Ch5TB.setGeometry(QtCore.QRect(10, 30, 71, 31))
        self.Ch5TB.setObjectName("Ch5TB")
        self.label_10 = QtWidgets.QLabel(self.frame_7)
        self.label_10.setGeometry(QtCore.QRect(10, 10, 91, 16))
        self.label_10.setObjectName("label_10")
        self.verticalLayout_5.addWidget(self.frame_7)
        self.frame_18 = QtWidgets.QFrame(self.frame_2)
        self.frame_18.setMinimumSize(QtCore.QSize(200, 0))
        self.frame_18.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_18.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_18.setObjectName("frame_18")
        self.Ch6TB = QtWidgets.QTextBrowser(self.frame_18)
        self.Ch6TB.setGeometry(QtCore.QRect(10, 30, 71, 31))
        self.Ch6TB.setObjectName("Ch6TB")
        self.label_21 = QtWidgets.QLabel(self.frame_18)
        self.label_21.setGeometry(QtCore.QRect(90, 40, 47, 13))
        self.label_21.setObjectName("label_21")
        self.label_22 = QtWidgets.QLabel(self.frame_18)
        self.label_22.setGeometry(QtCore.QRect(10, 10, 131, 16))
        self.label_22.setObjectName("label_22")
        self.verticalLayout_5.addWidget(self.frame_18)
        self.frame_14 = QtWidgets.QFrame(self.frame_2)
        self.frame_14.setMinimumSize(QtCore.QSize(200, 0))
        self.frame_14.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_14.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_14.setObjectName("frame_14")
        self.Ch7TB = QtWidgets.QTextBrowser(self.frame_14)
        self.Ch7TB.setGeometry(QtCore.QRect(10, 30, 71, 31))
        self.Ch7TB.setObjectName("Ch7TB")
        self.label_11 = QtWidgets.QLabel(self.frame_14)
        self.label_11.setGeometry(QtCore.QRect(90, 40, 47, 13))
        self.label_11.setObjectName("label_11")
        self.label_12 = QtWidgets.QLabel(self.frame_14)
        self.label_12.setGeometry(QtCore.QRect(10, 10, 151, 16))
        self.label_12.setObjectName("label_12")
        self.verticalLayout_5.addWidget(self.frame_14)
        self.horizontalLayout_2.addLayout(self.verticalLayout_5)
        self.verticalLayout.addWidget(self.frame_2)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setSpacing(6)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.frame = QtWidgets.QFrame(self.centralWidget)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame.setObjectName("frame")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.frame)
        self.verticalLayout_3.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_3.setSpacing(6)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.frame_4 = QtWidgets.QFrame(self.frame)
        self.frame_4.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_4.setObjectName("frame_4")
        self.verticalLayout_8 = QtWidgets.QVBoxLayout(self.frame_4)
        self.verticalLayout_8.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_8.setSpacing(6)
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.PlotLayout1 = QtWidgets.QVBoxLayout()
        self.PlotLayout1.setSpacing(6)
        self.PlotLayout1.setObjectName("PlotLayout1")
        self.verticalLayout_8.addLayout(self.PlotLayout1)
        self.verticalLayout_3.addWidget(self.frame_4)
        self.frame_5 = QtWidgets.QFrame(self.frame)
        self.frame_5.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_5.setObjectName("frame_5")
        self.verticalLayout_10 = QtWidgets.QVBoxLayout(self.frame_5)
        self.verticalLayout_10.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_10.setSpacing(6)
        self.verticalLayout_10.setObjectName("verticalLayout_10")
        self.PlotLayout2 = QtWidgets.QVBoxLayout()
        self.PlotLayout2.setSpacing(6)
        self.PlotLayout2.setObjectName("PlotLayout2")
        self.verticalLayout_10.addLayout(self.PlotLayout2)
        self.verticalLayout_3.addWidget(self.frame_5)
        self.frame_12 = QtWidgets.QFrame(self.frame)
        self.frame_12.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_12.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_12.setObjectName("frame_12")
        self.verticalLayout_12 = QtWidgets.QVBoxLayout(self.frame_12)
        self.verticalLayout_12.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_12.setSpacing(6)
        self.verticalLayout_12.setObjectName("verticalLayout_12")
        self.PlotLayout3 = QtWidgets.QVBoxLayout()
        self.PlotLayout3.setSpacing(6)
        self.PlotLayout3.setObjectName("PlotLayout3")
        self.verticalLayout_12.addLayout(self.PlotLayout3)
        self.verticalLayout_3.addWidget(self.frame_12)
        self.frame_13 = QtWidgets.QFrame(self.frame)
        self.frame_13.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_13.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_13.setObjectName("frame_13")
        self.verticalLayout_14 = QtWidgets.QVBoxLayout(self.frame_13)
        self.verticalLayout_14.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_14.setSpacing(6)
        self.verticalLayout_14.setObjectName("verticalLayout_14")
        self.PlotLayout4 = QtWidgets.QVBoxLayout()
        self.PlotLayout4.setSpacing(6)
        self.PlotLayout4.setObjectName("PlotLayout4")
        self.verticalLayout_14.addLayout(self.PlotLayout4)
        self.verticalLayout_3.addWidget(self.frame_13)
        self.frame_17 = QtWidgets.QFrame(self.frame)
        self.frame_17.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_17.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_17.setObjectName("frame_17")
        self.verticalLayout_9 = QtWidgets.QVBoxLayout(self.frame_17)
        self.verticalLayout_9.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_9.setSpacing(6)
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.PlotLayout5 = QtWidgets.QVBoxLayout()
        self.PlotLayout5.setSpacing(6)
        self.PlotLayout5.setObjectName("PlotLayout5")
        self.verticalLayout_9.addLayout(self.PlotLayout5)
        self.verticalLayout_3.addWidget(self.frame_17)
        self.frame_3 = QtWidgets.QFrame(self.frame)
        self.frame_3.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_3.setObjectName("frame_3")
        self.verticalLayout_16 = QtWidgets.QVBoxLayout(self.frame_3)
        self.verticalLayout_16.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_16.setSpacing(6)
        self.verticalLayout_16.setObjectName("verticalLayout_16")
        self.PlotLayout6 = QtWidgets.QVBoxLayout()
        self.PlotLayout6.setSpacing(6)
        self.PlotLayout6.setObjectName("PlotLayout6")
        self.verticalLayout_16.addLayout(self.PlotLayout6)
        self.verticalLayout_3.addWidget(self.frame_3)
        self.frame_19 = QtWidgets.QFrame(self.frame)
        self.frame_19.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_19.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame_19.setObjectName("frame_19")
        self.verticalLayout_11 = QtWidgets.QVBoxLayout(self.frame_19)
        self.verticalLayout_11.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_11.setSpacing(6)
        self.verticalLayout_11.setObjectName("verticalLayout_11")
        self.PlotLayout7 = QtWidgets.QVBoxLayout()
        self.PlotLayout7.setSpacing(6)
        self.PlotLayout7.setObjectName("PlotLayout7")
        self.verticalLayout_11.addLayout(self.PlotLayout7)
        self.verticalLayout_3.addWidget(self.frame_19)
        self.verticalLayout_2.addWidget(self.frame)
        self.horizontalLayout.addLayout(self.verticalLayout_2)
        my_pressures.setCentralWidget(self.centralWidget)
        self.menuBar = QtWidgets.QMenuBar(my_pressures)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 924, 21))
        self.menuBar.setObjectName("menuBar")
        self.menuFile = QtWidgets.QMenu(self.menuBar)
        self.menuFile.setObjectName("menuFile")
        self.menuSave = QtWidgets.QMenu(self.menuFile)
        self.menuSave.setObjectName("menuSave")
        my_pressures.setMenuBar(self.menuBar)
        self.mainToolBar = QtWidgets.QToolBar(my_pressures)
        self.mainToolBar.setObjectName("mainToolBar")
        my_pressures.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        self.statusBar = QtWidgets.QStatusBar(my_pressures)
        self.statusBar.setObjectName("statusBar")
        my_pressures.setStatusBar(self.statusBar)
        self.actionYou_wish = QtWidgets.QAction(my_pressures)
        self.actionYou_wish.setCheckable(True)
        self.actionYou_wish.setObjectName("actionYou_wish")
        self.menuSave.addAction(self.actionYou_wish)
        self.menuFile.addAction(self.menuSave.menuAction())
        self.menuBar.addAction(self.menuFile.menuAction())

        self.retranslateUi(my_pressures)
        QtCore.QMetaObject.connectSlotsByName(my_pressures)

    def retranslateUi(self, my_pressures):
        _translate = QtCore.QCoreApplication.translate
        my_pressures.setWindowTitle(_translate("my_pressures", "My Pressures"))
        self.label_16.setText(_translate("my_pressures", "Port Pressure"))
        self.ClearButton.setText(_translate("my_pressures", "Clear graphs"))
        self.label_17.setText(_translate("my_pressures", "Plot mode"))
        self.PlotBox.setItemText(0, _translate("my_pressures", "window"))
        self.PlotBox.setItemText(1, _translate("my_pressures", "all"))
        self.label_18.setText(_translate("my_pressures", "time window size (dd:hh:mm:ss)"))
        self.label_19.setText(_translate("my_pressures", "time step (mm:ss)"))
        self.LogScaleCB.setText(_translate("my_pressures", "LogScale"))
        self.COMButton.setText(_translate("my_pressures", "Check COM"))
        self.ConCOMButton.setText(_translate("my_pressures", "connect"))
        self.label_20.setText(_translate("my_pressures", "Port Temperature"))
        self.label.setText(_translate("my_pressures", "Ch 1 (Source)"))
        self.label_2.setText(_translate("my_pressures", "mbar"))
        self.label_3.setText(_translate("my_pressures", "mbar"))
        self.label_4.setText(_translate("my_pressures", "Ch 2 (Dec)"))
        self.label_5.setText(_translate("my_pressures", "mbar"))
        self.label_6.setText(_translate("my_pressures", "Ch 3 (Trap)"))
        self.label_7.setText(_translate("my_pressures", "mbar"))
        self.label_8.setText(_translate("my_pressures", "Ch 4 (Foreline 1)"))
        self.label_9.setText(_translate("my_pressures", "mbar"))
        self.label_10.setText(_translate("my_pressures", "Ch 5 (Foreline 2)"))
        self.label_21.setText(_translate("my_pressures", "mbar"))
        self.label_22.setText(_translate("my_pressures", "Ch 6 (Nothing so far...)"))
        self.label_11.setText(_translate("my_pressures", "degC"))
        self.label_12.setText(_translate("my_pressures", "Ch 7 (Chamber Temperature)"))
        self.menuFile.setTitle(_translate("my_pressures", "File"))
        self.menuSave.setTitle(_translate("my_pressures", "save"))
        self.actionYou_wish.setText(_translate("my_pressures", "you wish..."))
