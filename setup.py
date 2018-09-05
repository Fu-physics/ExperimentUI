from __future__ import unicode_literals
import sys
import os
import random
import time

import numpy as np
import xarray as xr


from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtWidgets import QInputDialog, QPushButton, QMainWindow, QApplication, QSpinBox, QLabel
from PyQt5.QtWidgets import QWidget, QAction, QTabWidget, QVBoxLayout, QHBoxLayout
from PyQt5.QtWidgets import QGroupBox, QDialog, QGridLayout
from PyQt5.QtWidgets import QApplication, QWidget, QLineEdit, QFileDialog
from PyQt5.QtWidgets import QPlainTextEdit
from PyQt5.QtWidgets import QFileDialog


import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
# Make sure that we are using QT5BB
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
#plt.ion()
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib import cm


from scipy.stats import poisson
from statistics import mean
import seaborn as sns

import sys


#home made packages
from ThorlabsCamera import ThorlabsCamer, ReadData

from DataProcess import DataProcess

#from DataProcess import ReadData

#from MatplotlibEmbeddingQt5 import MyMplCanvas, MyStaticMplCanvas, MyDynamicMplCanvas

#Gui package

from test1 import Ui_MainWindow



import pyqtgraph as pg


progname = os.path.basename(sys.argv[0])
progversion = "0.1"


class MyUi(QMainWindow, Ui_MainWindow):
    def __init__(self, parent = None):
        super (MyUi, self).__init__(parent)
        self.setupUi(self)
        self.DataBase()


###    Table 1 ###########################################################################
        self.pushButton_ConnectCamera.clicked.connect(self.ConnectCamera)
        self.pushButton_SetCamera.clicked.connect(self.SetCamera)

        ## Camera imaging plot
        self.plotImaging = plt.figure("Imaging")
        axes_SingImaging = self.plotImaging.add_subplot(111)
        self.add_plotfigure(self.plotImaging, self.verticalLayout_ImageShow)
        ## Start and Stop button
        self.pushButton_StartImaging.clicked.connect(self.StartImaging(axes_SingImaging))
        self.pushButton_StopImaging.clicked.connect(self.StopImaging)

        ## get multi-frames button
        self.pushButton_RunMultiFrames.clicked.connect(self.RunMultiFrames)
        self.pushButton_Save2Nc.clicked.connect(self.Save2Nc)


###    Table 2 ###########################################################################
        
        ## open File to plot Distribution and Coherent Length
        self.pushButton_OpenFile.clicked.connect(self.OpenFile)

        ## plot AverageImage 

        self.figAverageImage = plt.figure("First Order Imaging")
        #axesTopL = self.plotFirstOrd.add_subplot(111)
        self.add_plotfigure(self.figAverageImage, self.verticalLayout_AverageImage)
        self.pushButton_AverageImage.clicked.connect(self.PlotAverageImage(self.figAverageImage))


        self.figDistribution = plt.figure("Distribution Plot")
        #axesBot = self.plotDistri.add_subplot(111)
        self.add_plotfigure(self.figDistribution, self.verticalLayout_Distribution)
         # Just some button connected to `plot` method
        self.pushButton_Distribution.clicked.connect(self.PlotDistribution(self.figDistribution))


        ## plot CoherentLength 
        self.figCoherentLength = plt.figure("CoherentLength Plot")
        #axes_CoherentLength = self.plotCoherentLength.add_subplot(111)
        self.add_plotfigure(self.figCoherentLength, self.verticalLayout_CoherentLength)
        self.pushButton_CoherentLength.clicked.connect(self.CoherentLength(self.figCoherentLength))


###    Table 3 ###########################################################################
    
        ## plot High Order 

        ## Open File
        self.pushButton_OpenBackground.clicked.connect(self.OpenBackground)
        self.pushButton_OpenSignal.clicked.connect(self.OpenSignal)

        #self.figNOderMemontBackground = plt.figure("NOderMemontBackground")
        #axesTopL = self.plotFirstOrd.add_subplot(111)
        #self.add_plotfigure(self.figNOderMemontBackground, self.verticalLayout_NOrdMomentBackground)
        self.pushButton_NOderMomentBackground.clicked.connect(self.NOderMemontBackground)
        self.pushButton_NOderMomentSignal.clicked.connect(self.NOderMomentSignal)

        self.pushButton_SaveBackground.clicked.connect(self.SaveBackground)
        self.pushButton_SaveSignal.clicked.connect(self.SaveSignal)


    def DataBase(self):
        self.MomentBackground = {}
        self.CumulantBackground = {}
        self.MomentSignal = {}
        self.CumulantSignal = {}


    def ConnectCamera(self):
        """
        connect the camera

        :return:
        """
        #try:
        self.Cam = ThorlabsCamer()
        self.Cam.ConnectCamera()
        print(self.Cam.cam.serial)
        self.textBrowser_SetMeasureInf.setTextColor(QtCore.Qt.gray)
        self.textBrowser_SetMeasureInf.append("serial: " + str(self.Cam.cam.serial))
        self.textBrowser_SetMeasureInf.append("model: " + str(self.Cam.cam.model))
        self.textBrowser_SetMeasureInf.setTextColor(QtCore.Qt.green)
        self.textBrowser_SetMeasureInf.append("Camera connected... ")

    def SetCamera(self):

        exposeTime = self.spinBox_ExpTime.value()
        width = self.spinBox_Width.value()
        xshift = self.spinBox_XShift.value()
        hight = self.spinBox_Hight.value()
        yshift = self.spinBox_Yshift.value()

        if (xshift - 16 )%4 == 0 and (yshift - 16 )%4 == 0 and width%4 == 0 and hight%4 == 0:
            self.Cam.SetCamera(yshift = yshift, xshift = xshift, hight = hight, width = width, exposeTime = exposeTime)

            self.textBrowser_SetMeasureInf.setTextColor(QtCore.Qt.gray)
            self.textBrowser_SetMeasureInf.append("xshift, yshift: (" + str(xshift) + "," + str(yshift) +")")
            self.textBrowser_SetMeasureInf.append("width, hight: (" + str(width) + "," + str(hight) +")")
            self.textBrowser_SetMeasureInf.append("exposeTime: " + str(exposeTime) + " us" )

            print("setting the camera ...")
            self.textBrowser_SetMeasureInf.setTextColor(QtCore.Qt.green)
            self.textBrowser_SetMeasureInf.append("camera set down")
        else:

            self.textBrowser_SetMeasureInf.setTextColor(QtCore.Qt.red)
            self.textBrowser_SetMeasureInf.append("Plz check the parameter : width % 4 ==0 && Hight % 4 ==0")
            self.textBrowser_SetMeasureInf.append("(Xshift - 16 ) % 4 ==0 && (Yshift - 16 ) % 4 ==0")

    def add_plotfigure(self, figureName, plot_layout):
        # self.figureName = plt.figure()                      # a figure instance to plot on
        # if put "plt.ion" on the head, which will make two more figures idependently.

        # this is the Canvas Widget that displays the `figure`, it takes the `figure` instance as a parameter to __init__
        canvas_figureName = FigureCanvas(figureName)
        toolbar_figureName = NavigationToolbar(canvas_figureName,
                                               self)  # this is the Navigation widget, it takes the Canvas widget and a parent

        plot_layout.addWidget(toolbar_figureName)  # this also needed to show the Navigation of plot
        plot_layout.addWidget(canvas_figureName)  # add Canvas Widget(plot widget) onto tab_2

    
    def StartImaging(self, ax):

        """
        1) get each frame and plot it
        2) auto-update the plot

        :param ax: plt.addsubplot(111)
        :return: auto-update the plot
        """

        def inner_SingleImage_fun():

            def update_data():
                #print("AutoFresh  the figure !")

                #print("begin to get SingleImageData data ... ")
                singledata = self.Cam.SingleImageData(self.textBrowser_SetMeasureInf)
                #print("The image data shape is: ", np.shape(singledata))
                cax.set_data(singledata)

                #We need to draw *and* flush
                self.plotImaging.canvas.draw()
                self.plotImaging.canvas.flush_events()

            self.timer = QtCore.QTimer(self)
            #print("begin to plot !")
            self.textBrowser_SetMeasureInf.setTextColor(QtCore.Qt.green)
            self.textBrowser_SetMeasureInf.append("Start Imaging . . .")

            try:
                self.cbar.remove()
                #print("clear self.cbar !")
            except:
                pass
                #print("fail to clear self.cbar !")

            ax.cla()
            singledata = self.Cam.SingleImageData(self.textBrowser_SetMeasureInf)
            cax = ax.imshow(singledata, interpolation='nearest')
            ax.set_title('CMOS Camera')
            self.cbar = self.plotImaging.colorbar(cax, orientation='vertical')

            #cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar
            self.timer.timeout.connect(update_data)
            self.timer.start(100)

        return inner_SingleImage_fun


    def StopImaging(self):
        """
        Stop the auto-plot
        """
        self.timer.stop()
        self.textBrowser_SetMeasureInf.setTextColor(QtCore.Qt.red)
        self.textBrowser_SetMeasureInf.append("Stop Imaging!")



    def RunMultiFrames(self):
        """
        get the multi-fram image data, which are saved into *.npy files

        """


        frameNumber = self.spinBox_FrameNum.value()

        segmentNumber = self.spinBox_SegmentNum.value()


        self.textBrowser_SetMeasureInf.setTextColor(QtCore.Qt.gray)
        self.textBrowser_SetMeasureInf.append("frameNumber, segmentNumber is: "+ str(frameNumber)+ ", " + str(segmentNumber))
        self.textBrowser_SetMeasureInf.setTextColor(QtCore.Qt.green)
        self.textBrowser_SetMeasureInf.append("Running to get MultiImageData data ... ")
        print("frameNumber, segmentNumber is: ", frameNumber, segmentNumber)
        print("begin to get MultiImageData data ... ")
        self.Cam.MultiImageData(infoObj = self.textBrowser_SetMeasureInf,  frame_number_expected = frameNumber, segment_frame = segmentNumber)


    def Save2Nc(self):

        """
        1, transfer *.npy to *.nc data, which includes the data dimension info and experiment notes
        2, delete *.npy file
        2, save the *.nc file
        """

        frameNumber = self.spinBox_FrameNum.value()

        segmentNumber = self.spinBox_SegmentNum.value()

        exposeTime = self.spinBox_ExpTime.value()
        width = self.spinBox_Width.value()
        xshift = self.spinBox_XShift.value()
        hight = self.spinBox_Hight.value()
        yshift = self.spinBox_Yshift.value()

        print("frameNumber, segmentNumber, width, high is: ", frameNumber, segmentNumber, width, hight)
        app = ReadData(noteObj = self.textBrowser_SetMeasureInf, frameNumber=frameNumber, segmentFrame=segmentNumber, width=width, hight=hight)
        self.multiFrameData = app.ImageData()

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        # it just provides the name of file that you want to write into
        fileName, _= QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;NC Files (*.nc)", options=options)
       
        if fileName:
            print(fileName)

        self.multiFrameData.to_netcdf(fileName + '.nc')
        self.textBrowser_SetMeasureInf.setTextColor(QtCore.Qt.green)
        self.textBrowser_SetMeasureInf.append("the data has saved as .nc file! ")


    def OpenFile(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.nc)", options=options)
        if self.fileName:
            print(self.fileName)
        try:
            self.APP_dataprocess = DataProcess(self.fileName)
        except:
            print("Load file fall, Plz open file again ... ")




    def PlotAverageImage(self, fig):

        def inner_PlotFirstOrdrfun():
            print("inner_PlotFirstOrdrfun")

            # As self.firstOrdImaging will be used later, so we have to save as self.
            self.firstOrdImaging = self.APP_dataprocess.Average_Fluctuation() 

            
            #try:
            #    self.firstOrdImaging = self.APP_dataprocess.Average_Fluctuation()
            #    print(np.shape(self.firstOrdImaging))
            #except Exception:
            #    print("self.firstOrdImaging error!")

            # plot the Averge imaging
            ax = fig.add_subplot(111)

            try:
                self.cbar_FOrder.remove()
                ax.cla()
                #print("clear self.cbar !")
            except:
                pass
                #print("fail to clear self.cbar !")


            im = ax.imshow(self.firstOrdImaging)
            # create an axes on the right side of ax. The width of cax will be 5%
            # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            self.cbar_FOrder =  plt.colorbar(im, cax=cax)
            #plt.colorbar(im, cax=cax, ticks=[0, 5, 10])
            ax.set_title('1th Order')

            plt.savefig('1th Order Imaing.eps', format='eps', dpi=100)
            plt.close()
        
        return inner_PlotFirstOrdrfun



    def PlotDistribution(self, fig):

        def BoseEinstein(Nbar, n = 51):
            nList = np.linspace(0, n, n+1, dtype = int)
            result = 1/(1+Nbar)*(Nbar/(1+Nbar))**nList
            return result


        def G2(FirstArray1d, SecondArray1d):

            Frames = len(FirstArray1d)
            #print("Frames is: " , Frames)
            FirstAverage = np.sum(FirstArray1d)/Frames
            SecondAverage = np.sum(SecondArray1d)/Frames
            
            result = np.sum(FirstArray1d*SecondArray1d)/ Frames /(FirstAverage * SecondAverage)
            return result

        def inner_PlotDistrifun():

            """
            Plot the photon number distribution ...
            """
                    
            font = {'family': 'serif',
                    'color':  'darkred',
                    'weight': 'normal',
                    'size': 16}

            Nmax = 100
            bins = np.linspace(0, Nmax, Nmax+1)
            nList = np.linspace(0, Nmax, Nmax+1, dtype = int)

            y_location = self.spinBox_PixelY.value()
            x_location = self.spinBox_PixelX.value()

            # get pixel intensity data
            Array1 = self.APP_dataprocess.PixelData(y_location, x_location)
            Array2 = Array1
            g2 = G2(Array1, Array2)
            print("g2 is:", g2)

            arr = []
            rv = poisson(self.firstOrdImaging[y_location, x_location])
            for num in range(0,40):
                arr.append(rv.pmf(num))

            ax = fig.add_subplot(111)

            try:
                ax.cla()
                #print("clear self.cbar !")
            except:
                pass
                #print("fail to clear self.cbar !")
            
            ax.hist(Array1 , bins, normed=True, label = "Data distribution") 
            ax.plot(nList, BoseEinstein(self.firstOrdImaging[y_location, x_location], Nmax), label ="BoseEinstein distribution")
            ax.plot(arr, linewidth=2.0, label ="Possion distribution")
            ax.set_title("Pixel Position({},{}); <$I$>:{}".format(x_location , y_location, self.firstOrdImaging[y_location, x_location]), fontdict=font)
            
            ax.text(22, .08, r"g2:{}".format(g2), fontdict=font)
            ax.legend() 
            
            fig.savefig('PixelPosition({},{})PhotDist.eps'.format(x_location , y_location), format='eps', dpi=300)
            plt.close()

        return inner_PlotDistrifun


    def CoherentLength(self, fig):

        def inner_CoherentLength():

            """
            Plot the photon number distribution ...
            """

            ax = fig.add_subplot(111)

            try:
                ax.cla()
                #print("clear self.cbar !")
            except:
                pass
                #print("fail to clear self.cbar !")
            xCorr, yCorr = self.APP_dataprocess.SpatialCorrelation([self.spinBox_PixelX.value(), self.spinBox_PixelY.value()])
            ax.plot(xCorr)
            ax.set_title("G2 @({}, {})".format(self.spinBox_PixelX.value(), self.spinBox_PixelY.value()))
            fig.savefig("G2 @({}, {}).png".format(self.spinBox_PixelX.value(), self.spinBox_PixelY.value()), format="png", dpi = 100)
            plt.close()

        return inner_CoherentLength
    
    def OpenBackground(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.nc)", options=options)
        if fileName:
            print(fileName)
        try:
            self.Background_dataprocess = DataProcess(fileName)
        except:
            print("Load file fall, Plz open file again ... ")

    def OpenSignal(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.nc)", options=options)
        if fileName:
            print(fileName)
        try:
            self.Signal_dataprocess = DataProcess(fileName)
        except:
            print("Load file fall, Plz open file again ... ")


    def NOderMemontBackground(self):

        #def inner_NOderMemontBackground():
        NorderValue = self.spinBox_NOderMomentBackground.value()
        print(self.Background_dataprocess)

        try:
            self.corbar.remove()   # here "self." must be added !
        except:
            pass

        try:
               ## befor calculate the High Order moment, we need to get the Aaverage and Fluctuation
            self.Background_dataprocess.Average_Fluctuation()
            NorderOrdImaging = self.Background_dataprocess.NOrder(NorderValue)

            self.MomentBackground["{}".format(NorderValue)] = NorderOrdImaging
            #imageData = np.array([[1, 2, 3], [4, 5, 6]], np.int32)  # used for testing


            myplotSt = self.widget_Background.mpl
            cax = myplotSt.axes.imshow(NorderOrdImaging) #, clim=(0.0, 10))
            myplotSt.axes.set_title('Title')

            # Add colorbar, make sure to specify tick locations to match desired ticklabels
            self.corbar = myplotSt.fig.colorbar(cax)  #, ticks=[0, 1])

            #plt.close()
            #return inner_NOderMemontBackground
        except:
            print("Background file not opened, Plz open Background file . . . ")


    def NOderMomentSignal(self):
        pass

            
    def SaveBackground(self):
        pass
        #np.save("MomentBackground", self.MomentBackground)
        #np.save("CumulantBackground", self.CumulantBackground)


    def SaveSignal(self):
        pass
        #np.save("MomentSignal", self.MomentBackground)
        #np.save("CumulantSignal", self.CumulantBackground)








if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    ui = MyUi()
    ui.show()
    sys.exit(app.exec_())

