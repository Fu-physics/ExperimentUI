# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\test1.ui'
#
# Created by: PyQt5 UI code generator 5.11.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1008, 809)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout()
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.tab)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.splitter_2 = QtWidgets.QSplitter(self.tab)
        self.splitter_2.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_2.setObjectName("splitter_2")
        self.verticalLayoutWidget = QtWidgets.QWidget(self.splitter_2)
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox = QtWidgets.QGroupBox(self.verticalLayoutWidget)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.pushButton_ConnectCamera = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_ConnectCamera.setObjectName("pushButton_ConnectCamera")
        self.verticalLayout_3.addWidget(self.pushButton_ConnectCamera)
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.spinBox_ExpTime = QtWidgets.QSpinBox(self.groupBox)
        self.spinBox_ExpTime.setMaximum(1000000)
        self.spinBox_ExpTime.setProperty("value", 10)
        self.spinBox_ExpTime.setObjectName("spinBox_ExpTime")
        self.gridLayout_2.addWidget(self.spinBox_ExpTime, 0, 2, 1, 2)
        self.label_10 = QtWidgets.QLabel(self.groupBox)
        self.label_10.setObjectName("label_10")
        self.gridLayout_2.addWidget(self.label_10, 0, 0, 1, 2)
        self.spinBox_Width = QtWidgets.QSpinBox(self.groupBox)
        self.spinBox_Width.setMaximum(1280)
        self.spinBox_Width.setProperty("value", 1280)
        self.spinBox_Width.setObjectName("spinBox_Width")
        self.gridLayout_2.addWidget(self.spinBox_Width, 1, 1, 1, 1)
        self.spinBox_Hight = QtWidgets.QSpinBox(self.groupBox)
        self.spinBox_Hight.setMaximum(1024)
        self.spinBox_Hight.setProperty("value", 1024)
        self.spinBox_Hight.setObjectName("spinBox_Hight")
        self.gridLayout_2.addWidget(self.spinBox_Hight, 2, 1, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.groupBox)
        self.label_7.setObjectName("label_7")
        self.gridLayout_2.addWidget(self.label_7, 1, 2, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.groupBox)
        self.label_6.setObjectName("label_6")
        self.gridLayout_2.addWidget(self.label_6, 1, 0, 1, 1)
        self.spinBox_Yshift = QtWidgets.QSpinBox(self.groupBox)
        self.spinBox_Yshift.setMaximum(1024)
        self.spinBox_Yshift.setObjectName("spinBox_Yshift")
        self.gridLayout_2.addWidget(self.spinBox_Yshift, 2, 3, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.groupBox)
        self.label_8.setObjectName("label_8")
        self.gridLayout_2.addWidget(self.label_8, 2, 2, 1, 1)
        self.spinBox_XShift = QtWidgets.QSpinBox(self.groupBox)
        self.spinBox_XShift.setMaximum(1280)
        self.spinBox_XShift.setObjectName("spinBox_XShift")
        self.gridLayout_2.addWidget(self.spinBox_XShift, 1, 3, 1, 1)
        self.label_9 = QtWidgets.QLabel(self.groupBox)
        self.label_9.setObjectName("label_9")
        self.gridLayout_2.addWidget(self.label_9, 2, 0, 1, 1)
        self.verticalLayout_3.addLayout(self.gridLayout_2)
        self.pushButton_SetCamera = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_SetCamera.setObjectName("pushButton_SetCamera")
        self.verticalLayout_3.addWidget(self.pushButton_SetCamera)
        self.label_3 = QtWidgets.QLabel(self.groupBox)
        self.label_3.setObjectName("label_3")
        self.verticalLayout_3.addWidget(self.label_3)
        self.textBrowser_SetMeasureInf = QtWidgets.QTextBrowser(self.groupBox)
        self.textBrowser_SetMeasureInf.setObjectName("textBrowser_SetMeasureInf")
        self.verticalLayout_3.addWidget(self.textBrowser_SetMeasureInf)
        self.gridLayout_4 = QtWidgets.QGridLayout()
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.spinBox_FrameNum = QtWidgets.QSpinBox(self.groupBox)
        self.spinBox_FrameNum.setMaximum(500000)
        self.spinBox_FrameNum.setProperty("value", 5000)
        self.spinBox_FrameNum.setObjectName("spinBox_FrameNum")
        self.gridLayout_4.addWidget(self.spinBox_FrameNum, 0, 1, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.label_2.setObjectName("label_2")
        self.gridLayout_4.addWidget(self.label_2, 1, 0, 1, 1)
        self.label = QtWidgets.QLabel(self.groupBox)
        self.label.setObjectName("label")
        self.gridLayout_4.addWidget(self.label, 0, 0, 1, 1)
        self.pushButton_RunMultiFrames = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_RunMultiFrames.setObjectName("pushButton_RunMultiFrames")
        self.gridLayout_4.addWidget(self.pushButton_RunMultiFrames, 2, 0, 1, 1)
        self.spinBox_SegmentNum = QtWidgets.QSpinBox(self.groupBox)
        self.spinBox_SegmentNum.setMaximum(100000)
        self.spinBox_SegmentNum.setProperty("value", 1000)
        self.spinBox_SegmentNum.setObjectName("spinBox_SegmentNum")
        self.gridLayout_4.addWidget(self.spinBox_SegmentNum, 1, 1, 1, 1)
        self.pushButton_Save2Nc = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_Save2Nc.setObjectName("pushButton_Save2Nc")
        self.gridLayout_4.addWidget(self.pushButton_Save2Nc, 3, 0, 1, 2)
        self.verticalLayout_3.addLayout(self.gridLayout_4)
        self.verticalLayout.addWidget(self.groupBox)
        self.verticalLayoutWidget_2 = QtWidgets.QWidget(self.splitter_2)
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 1, 1, 1, 1)
        self.pushButton_StopImaging = QtWidgets.QPushButton(self.verticalLayoutWidget_2)
        font = QtGui.QFont()
        font.setBold(True)
        font.setUnderline(True)
        font.setWeight(75)
        font.setStrikeOut(False)
        self.pushButton_StopImaging.setFont(font)
        self.pushButton_StopImaging.setStyleSheet("background-color: #A3C1DA; color: red")
        self.pushButton_StopImaging.setObjectName("pushButton_StopImaging")
        self.gridLayout.addWidget(self.pushButton_StopImaging, 1, 2, 1, 1)
        self.pushButton_StartImaging = QtWidgets.QPushButton(self.verticalLayoutWidget_2)
        font = QtGui.QFont()
        font.setBold(True)
        font.setUnderline(True)
        font.setWeight(75)
        self.pushButton_StartImaging.setFont(font)
        self.pushButton_StartImaging.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))
        self.pushButton_StartImaging.setAutoFillBackground(False)
        self.pushButton_StartImaging.setStyleSheet("background-color: #A3C1DA; color: green")
        self.pushButton_StartImaging.setObjectName("pushButton_StartImaging")
        self.gridLayout.addWidget(self.pushButton_StartImaging, 1, 0, 1, 1)
        self.verticalLayout_ImageShow = QtWidgets.QVBoxLayout()
        self.verticalLayout_ImageShow.setObjectName("verticalLayout_ImageShow")
        self.gridLayout.addLayout(self.verticalLayout_ImageShow, 0, 0, 1, 3)
        self.verticalLayout_2.addLayout(self.gridLayout)
        self.horizontalLayout.addWidget(self.splitter_2)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.tab_2)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.splitter = QtWidgets.QSplitter(self.tab_2)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.verticalLayoutWidget_3 = QtWidgets.QWidget(self.splitter)
        self.verticalLayoutWidget_3.setObjectName("verticalLayoutWidget_3")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_3)
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.pushButton_OpenFile = QtWidgets.QPushButton(self.verticalLayoutWidget_3)
        self.pushButton_OpenFile.setObjectName("pushButton_OpenFile")
        self.verticalLayout_4.addWidget(self.pushButton_OpenFile)
        self.verticalLayout_AverageImage = QtWidgets.QVBoxLayout()
        self.verticalLayout_AverageImage.setObjectName("verticalLayout_AverageImage")
        self.verticalLayout_4.addLayout(self.verticalLayout_AverageImage)
        self.pushButton_AverageImage = QtWidgets.QPushButton(self.verticalLayoutWidget_3)
        self.pushButton_AverageImage.setObjectName("pushButton_AverageImage")
        self.verticalLayout_4.addWidget(self.pushButton_AverageImage)
        self.label_4 = QtWidgets.QLabel(self.verticalLayoutWidget_3)
        self.label_4.setObjectName("label_4")
        self.verticalLayout_4.addWidget(self.label_4)
        self.textBrowser = QtWidgets.QTextBrowser(self.verticalLayoutWidget_3)
        self.textBrowser.setObjectName("textBrowser")
        self.verticalLayout_4.addWidget(self.textBrowser)
        self.verticalLayoutWidget_4 = QtWidgets.QWidget(self.splitter)
        self.verticalLayoutWidget_4.setObjectName("verticalLayoutWidget_4")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.verticalLayoutWidget_4)
        self.gridLayout_3.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.widget_2 = QtWidgets.QWidget(self.verticalLayoutWidget_4)
        self.widget_2.setObjectName("widget_2")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.widget_2)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.verticalLayout_8 = QtWidgets.QVBoxLayout()
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.widget_4 = QtWidgets.QWidget(self.widget_2)
        self.widget_4.setObjectName("widget_4")
        self.verticalLayoutWidget_7 = QtWidgets.QWidget(self.widget_4)
        self.verticalLayoutWidget_7.setGeometry(QtCore.QRect(10, 10, 421, 291))
        self.verticalLayoutWidget_7.setObjectName("verticalLayoutWidget_7")
        self.verticalLayout_Distribution = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_7)
        self.verticalLayout_Distribution.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_Distribution.setObjectName("verticalLayout_Distribution")
        self.verticalLayout_8.addWidget(self.widget_4)
        self.pushButton_Distribution = QtWidgets.QPushButton(self.widget_2)
        self.pushButton_Distribution.setObjectName("pushButton_Distribution")
        self.verticalLayout_8.addWidget(self.pushButton_Distribution)
        self.verticalLayout_6.addLayout(self.verticalLayout_8)
        self.verticalLayout_16 = QtWidgets.QVBoxLayout()
        self.verticalLayout_16.setObjectName("verticalLayout_16")
        self.widget_5 = QtWidgets.QWidget(self.widget_2)
        self.widget_5.setObjectName("widget_5")
        self.verticalLayoutWidget_8 = QtWidgets.QWidget(self.widget_5)
        self.verticalLayoutWidget_8.setGeometry(QtCore.QRect(10, 10, 421, 291))
        self.verticalLayoutWidget_8.setObjectName("verticalLayoutWidget_8")
        self.verticalLayout_CoherentLength = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_8)
        self.verticalLayout_CoherentLength.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_CoherentLength.setObjectName("verticalLayout_CoherentLength")
        self.verticalLayout_16.addWidget(self.widget_5)
        self.pushButton_CoherentLength = QtWidgets.QPushButton(self.widget_2)
        self.pushButton_CoherentLength.setObjectName("pushButton_CoherentLength")
        self.verticalLayout_16.addWidget(self.pushButton_CoherentLength)
        self.verticalLayout_6.addLayout(self.verticalLayout_16)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem1)
        self.label_5 = QtWidgets.QLabel(self.widget_2)
        self.label_5.setObjectName("label_5")
        self.horizontalLayout_4.addWidget(self.label_5)
        self.spinBox_PixelX = QtWidgets.QSpinBox(self.widget_2)
        self.spinBox_PixelX.setMaximum(1280)
        self.spinBox_PixelX.setObjectName("spinBox_PixelX")
        self.horizontalLayout_4.addWidget(self.spinBox_PixelX)
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem2)
        self.label_11 = QtWidgets.QLabel(self.widget_2)
        self.label_11.setObjectName("label_11")
        self.horizontalLayout_4.addWidget(self.label_11)
        self.spinBox_PixelY = QtWidgets.QSpinBox(self.widget_2)
        self.spinBox_PixelY.setMaximum(1024)
        self.spinBox_PixelY.setObjectName("spinBox_PixelY")
        self.horizontalLayout_4.addWidget(self.spinBox_PixelY)
        self.verticalLayout_6.addLayout(self.horizontalLayout_4)
        self.gridLayout_3.addWidget(self.widget_2, 1, 0, 1, 1)
        self.horizontalLayout_2.addWidget(self.splitter)
        self.tabWidget.addTab(self.tab_2, "")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.verticalLayout_12 = QtWidgets.QVBoxLayout(self.tab_3)
        self.verticalLayout_12.setObjectName("verticalLayout_12")
        self.verticalLayout_9 = QtWidgets.QVBoxLayout()
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.widget = QtWidgets.QWidget(self.tab_3)
        self.widget.setObjectName("widget")
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout(self.widget)
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.splitter_3 = QtWidgets.QSplitter(self.widget)
        self.splitter_3.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_3.setObjectName("splitter_3")
        self.verticalLayoutWidget_5 = QtWidgets.QWidget(self.splitter_3)
        self.verticalLayoutWidget_5.setObjectName("verticalLayoutWidget_5")
        self.verticalLayout_10 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_5)
        self.verticalLayout_10.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_10.setObjectName("verticalLayout_10")
        self.pushButton_OpenBackground = QtWidgets.QPushButton(self.verticalLayoutWidget_5)
        self.pushButton_OpenBackground.setObjectName("pushButton_OpenBackground")
        self.verticalLayout_10.addWidget(self.pushButton_OpenBackground)
        self.widget_Background = MatplotlibWidget(self.verticalLayoutWidget_5)
        self.widget_Background.setObjectName("widget_Background")
        self.verticalLayout_10.addWidget(self.widget_Background)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.spinBox_NOderMomentBackground = QtWidgets.QSpinBox(self.verticalLayoutWidget_5)
        self.spinBox_NOderMomentBackground.setMinimum(1)
        self.spinBox_NOderMomentBackground.setMaximum(20)
        self.spinBox_NOderMomentBackground.setObjectName("spinBox_NOderMomentBackground")
        self.horizontalLayout_3.addWidget(self.spinBox_NOderMomentBackground)
        self.label_13 = QtWidgets.QLabel(self.verticalLayoutWidget_5)
        self.label_13.setObjectName("label_13")
        self.horizontalLayout_3.addWidget(self.label_13)
        self.pushButton_NOderMomentBackground = QtWidgets.QPushButton(self.verticalLayoutWidget_5)
        self.pushButton_NOderMomentBackground.setObjectName("pushButton_NOderMomentBackground")
        self.horizontalLayout_3.addWidget(self.pushButton_NOderMomentBackground)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem3)
        self.verticalLayout_10.addLayout(self.horizontalLayout_3)
        self.pushButton_SaveBackground = QtWidgets.QPushButton(self.verticalLayoutWidget_5)
        self.pushButton_SaveBackground.setObjectName("pushButton_SaveBackground")
        self.verticalLayout_10.addWidget(self.pushButton_SaveBackground)
        self.verticalLayoutWidget_6 = QtWidgets.QWidget(self.splitter_3)
        self.verticalLayoutWidget_6.setObjectName("verticalLayoutWidget_6")
        self.verticalLayout_11 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_6)
        self.verticalLayout_11.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_11.setObjectName("verticalLayout_11")
        self.pushButton_OpenSignal = QtWidgets.QPushButton(self.verticalLayoutWidget_6)
        self.pushButton_OpenSignal.setObjectName("pushButton_OpenSignal")
        self.verticalLayout_11.addWidget(self.pushButton_OpenSignal)
        self.widget_Signal = MatplotlibWidget(self.verticalLayoutWidget_6)
        self.widget_Signal.setObjectName("widget_Signal")
        self.verticalLayout_11.addWidget(self.widget_Signal)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.spinBox_NOderMomentSignal = QtWidgets.QSpinBox(self.verticalLayoutWidget_6)
        self.spinBox_NOderMomentSignal.setEnabled(True)
        self.spinBox_NOderMomentSignal.setMinimum(1)
        self.spinBox_NOderMomentSignal.setMaximum(20)
        self.spinBox_NOderMomentSignal.setProperty("value", 1)
        self.spinBox_NOderMomentSignal.setObjectName("spinBox_NOderMomentSignal")
        self.horizontalLayout_5.addWidget(self.spinBox_NOderMomentSignal)
        self.label_14 = QtWidgets.QLabel(self.verticalLayoutWidget_6)
        self.label_14.setObjectName("label_14")
        self.horizontalLayout_5.addWidget(self.label_14)
        self.pushButton_NOderMomentSignal = QtWidgets.QPushButton(self.verticalLayoutWidget_6)
        self.pushButton_NOderMomentSignal.setObjectName("pushButton_NOderMomentSignal")
        self.horizontalLayout_5.addWidget(self.pushButton_NOderMomentSignal)
        spacerItem4 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_5.addItem(spacerItem4)
        self.verticalLayout_11.addLayout(self.horizontalLayout_5)
        self.pushButton_SaveSignal = QtWidgets.QPushButton(self.verticalLayoutWidget_6)
        self.pushButton_SaveSignal.setObjectName("pushButton_SaveSignal")
        self.verticalLayout_11.addWidget(self.pushButton_SaveSignal)
        self.horizontalLayout_6.addWidget(self.splitter_3)
        self.verticalLayout_9.addWidget(self.widget)
        self.verticalLayout_12.addLayout(self.verticalLayout_9)
        self.tabWidget.addTab(self.tab_3, "")
        self.tab_4 = QtWidgets.QWidget()
        self.tab_4.setObjectName("tab_4")
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout(self.tab_4)
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.verticalLayout_13 = QtWidgets.QVBoxLayout()
        self.verticalLayout_13.setObjectName("verticalLayout_13")
        self.pushButton_OpenBackgroundMC = QtWidgets.QPushButton(self.tab_4)
        self.pushButton_OpenBackgroundMC.setObjectName("pushButton_OpenBackgroundMC")
        self.verticalLayout_13.addWidget(self.pushButton_OpenBackgroundMC)
        self.pushButton_OpenSignalMC = QtWidgets.QPushButton(self.tab_4)
        self.pushButton_OpenSignalMC.setObjectName("pushButton_OpenSignalMC")
        self.verticalLayout_13.addWidget(self.pushButton_OpenSignalMC)
        self.pushButton_ImageProcess = QtWidgets.QPushButton(self.tab_4)
        self.pushButton_ImageProcess.setObjectName("pushButton_ImageProcess")
        self.verticalLayout_13.addWidget(self.pushButton_ImageProcess)
        self.widget_3 = QtWidgets.QWidget(self.tab_4)
        self.widget_3.setObjectName("widget_3")
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout(self.widget_3)
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.splitter_4 = QtWidgets.QSplitter(self.widget_3)
        self.splitter_4.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_4.setObjectName("splitter_4")
        self.verticalLayoutWidget_9 = QtWidgets.QWidget(self.splitter_4)
        self.verticalLayoutWidget_9.setObjectName("verticalLayoutWidget_9")
        self.verticalLayout_14 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_9)
        self.verticalLayout_14.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_14.setObjectName("verticalLayout_14")
        self.widget_NOderImage = MatplotlibWidget(self.verticalLayoutWidget_9)
        self.widget_NOderImage.setObjectName("widget_NOderImage")
        self.verticalLayout_14.addWidget(self.widget_NOderImage)
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        self.spinBox_NOderImage = QtWidgets.QSpinBox(self.verticalLayoutWidget_9)
        self.spinBox_NOderImage.setMinimum(1)
        self.spinBox_NOderImage.setMaximum(20)
        self.spinBox_NOderImage.setObjectName("spinBox_NOderImage")
        self.horizontalLayout_8.addWidget(self.spinBox_NOderImage)
        self.label_15 = QtWidgets.QLabel(self.verticalLayoutWidget_9)
        self.label_15.setObjectName("label_15")
        self.horizontalLayout_8.addWidget(self.label_15)
        spacerItem5 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_8.addItem(spacerItem5)
        self.doubleSpinBox_UplimitImage = QtWidgets.QDoubleSpinBox(self.verticalLayoutWidget_9)
        self.doubleSpinBox_UplimitImage.setProperty("value", 1.0)
        self.doubleSpinBox_UplimitImage.setObjectName("doubleSpinBox_UplimitImage")
        self.horizontalLayout_8.addWidget(self.doubleSpinBox_UplimitImage)
        self.label_12 = QtWidgets.QLabel(self.verticalLayoutWidget_9)
        self.label_12.setObjectName("label_12")
        self.horizontalLayout_8.addWidget(self.label_12)
        spacerItem6 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_8.addItem(spacerItem6)
        self.pushButton_NOderImage = QtWidgets.QPushButton(self.verticalLayoutWidget_9)
        self.pushButton_NOderImage.setObjectName("pushButton_NOderImage")
        self.horizontalLayout_8.addWidget(self.pushButton_NOderImage)
        spacerItem7 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_8.addItem(spacerItem7)
        self.verticalLayout_14.addLayout(self.horizontalLayout_8)
        self.pushButton_SaveCumulant = QtWidgets.QPushButton(self.verticalLayoutWidget_9)
        self.pushButton_SaveCumulant.setObjectName("pushButton_SaveCumulant")
        self.verticalLayout_14.addWidget(self.pushButton_SaveCumulant)
        self.verticalLayoutWidget_10 = QtWidgets.QWidget(self.splitter_4)
        self.verticalLayoutWidget_10.setObjectName("verticalLayoutWidget_10")
        self.verticalLayout_15 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_10)
        self.verticalLayout_15.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_15.setObjectName("verticalLayout_15")
        self.widget_NOderImageNorm = MatplotlibWidget(self.verticalLayoutWidget_10)
        self.widget_NOderImageNorm.setObjectName("widget_NOderImageNorm")
        self.verticalLayout_15.addWidget(self.widget_NOderImageNorm)
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.spinBox_NOderImageNorm = QtWidgets.QSpinBox(self.verticalLayoutWidget_10)
        self.spinBox_NOderImageNorm.setEnabled(True)
        self.spinBox_NOderImageNorm.setMinimum(1)
        self.spinBox_NOderImageNorm.setMaximum(20)
        self.spinBox_NOderImageNorm.setProperty("value", 4)
        self.spinBox_NOderImageNorm.setObjectName("spinBox_NOderImageNorm")
        self.horizontalLayout_9.addWidget(self.spinBox_NOderImageNorm)
        self.label_16 = QtWidgets.QLabel(self.verticalLayoutWidget_10)
        self.label_16.setObjectName("label_16")
        self.horizontalLayout_9.addWidget(self.label_16)
        spacerItem8 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_9.addItem(spacerItem8)
        self.doubleSpinBox_UplimitImageNorm = QtWidgets.QDoubleSpinBox(self.verticalLayoutWidget_10)
        self.doubleSpinBox_UplimitImageNorm.setSingleStep(1e-05)
        self.doubleSpinBox_UplimitImageNorm.setProperty("value", 1.0)
        self.doubleSpinBox_UplimitImageNorm.setObjectName("doubleSpinBox_UplimitImageNorm")
        self.horizontalLayout_9.addWidget(self.doubleSpinBox_UplimitImageNorm)
        self.label_17 = QtWidgets.QLabel(self.verticalLayoutWidget_10)
        self.label_17.setObjectName("label_17")
        self.horizontalLayout_9.addWidget(self.label_17)
        spacerItem9 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_9.addItem(spacerItem9)
        self.pushButton_NOderImageNorm = QtWidgets.QPushButton(self.verticalLayoutWidget_10)
        self.pushButton_NOderImageNorm.setObjectName("pushButton_NOderImageNorm")
        self.horizontalLayout_9.addWidget(self.pushButton_NOderImageNorm)
        spacerItem10 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_9.addItem(spacerItem10)
        self.verticalLayout_15.addLayout(self.horizontalLayout_9)
        self.pushButton_SaveCumulantNorm = QtWidgets.QPushButton(self.verticalLayoutWidget_10)
        self.pushButton_SaveCumulantNorm.setObjectName("pushButton_SaveCumulantNorm")
        self.verticalLayout_15.addWidget(self.pushButton_SaveCumulantNorm)
        self.horizontalLayout_7.addWidget(self.splitter_4)
        self.verticalLayout_13.addWidget(self.widget_3)
        self.horizontalLayout_10.addLayout(self.verticalLayout_13)
        self.tabWidget.addTab(self.tab_4, "")
        self.verticalLayout_7.addWidget(self.tabWidget)
        self.verticalLayout_5.addLayout(self.verticalLayout_7)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1008, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(3)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.groupBox.setTitle(_translate("MainWindow", "Camera Setting"))
        self.pushButton_ConnectCamera.setText(_translate("MainWindow", "Connect Camera"))
        self.label_10.setText(_translate("MainWindow", "              Exposure Time"))
        self.label_7.setText(_translate("MainWindow", "         X shift"))
        self.label_6.setText(_translate("MainWindow", "              Width"))
        self.label_8.setText(_translate("MainWindow", "          Y shift"))
        self.label_9.setText(_translate("MainWindow", "              Hight"))
        self.pushButton_SetCamera.setText(_translate("MainWindow", "Setting Camera"))
        self.label_3.setText(_translate("MainWindow", "Setting and measurement Inf."))
        self.label_2.setText(_translate("MainWindow", "Segment Number"))
        self.label.setText(_translate("MainWindow", "Frame Number"))
        self.pushButton_RunMultiFrames.setText(_translate("MainWindow", "Run to Get Multi-Frames  Data"))
        self.pushButton_Save2Nc.setText(_translate("MainWindow", "SAVE to .nc"))
        self.pushButton_StopImaging.setText(_translate("MainWindow", "Stop"))
        self.pushButton_StartImaging.setText(_translate("MainWindow", "Start"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "Camera Setting & Imaging"))
        self.pushButton_OpenFile.setText(_translate("MainWindow", "OpenFile"))
        self.pushButton_AverageImage.setText(_translate("MainWindow", "Average Image"))
        self.label_4.setText(_translate("MainWindow", "TextLabel"))
        self.pushButton_Distribution.setText(_translate("MainWindow", "Distribution"))
        self.pushButton_CoherentLength.setText(_translate("MainWindow", "Coherent Length"))
        self.label_5.setText(_translate("MainWindow", "X"))
        self.label_11.setText(_translate("MainWindow", "Y"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "Coherent Proper."))
        self.pushButton_OpenBackground.setText(_translate("MainWindow", "Open Background File"))
        self.label_13.setText(_translate("MainWindow", "TextLabel"))
        self.pushButton_NOderMomentBackground.setText(_translate("MainWindow", "PushButton"))
        self.pushButton_SaveBackground.setText(_translate("MainWindow", "Save Background Data"))
        self.pushButton_OpenSignal.setText(_translate("MainWindow", "Open Signal File"))
        self.label_14.setText(_translate("MainWindow", "TextLabel"))
        self.pushButton_NOderMomentSignal.setText(_translate("MainWindow", "PushButton"))
        self.pushButton_SaveSignal.setText(_translate("MainWindow", "Save Signal  Data"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("MainWindow", "Moment"))
        self.pushButton_OpenBackgroundMC.setText(_translate("MainWindow", "Open Background Momtum ( Cumulant )"))
        self.pushButton_OpenSignalMC.setText(_translate("MainWindow", "Open Signal Momtum ( Cumulant )"))
        self.pushButton_ImageProcess.setText(_translate("MainWindow", "Image Processing"))
        self.label_15.setText(_translate("MainWindow", "TextLabel"))
        self.label_12.setText(_translate("MainWindow", "Uplimit"))
        self.pushButton_NOderImage.setText(_translate("MainWindow", "PushButton"))
        self.pushButton_SaveCumulant.setText(_translate("MainWindow", "Save Cumulant"))
        self.label_16.setText(_translate("MainWindow", "TextLabel"))
        self.label_17.setText(_translate("MainWindow", "TextLabel"))
        self.pushButton_NOderImageNorm.setText(_translate("MainWindow", "PushButton"))
        self.pushButton_SaveCumulantNorm.setText(_translate("MainWindow", "Save CumulantNorm"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_4), _translate("MainWindow", "Normalized Image"))

from MatplotlibWidget import MatplotlibWidget
