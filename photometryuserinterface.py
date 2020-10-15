import photometry
from PyQt5 import QtWidgets, uic, QtCore
from PyQt5.QtCore import QThread, pyqtSignal
from PyQt5.QtWidgets import QFileDialog, QDialog, QMessageBox, QApplication, QPushButton, QLabel, QCheckBox
from PyQt5.QtGui import QPixmap, QImage, QPainter, QPen, QColor, QFont
from PySide2 import QtUiTools
from astropy.io import fits
from astropy.visualization import ZScaleInterval, MinMaxInterval, ManualInterval
import sys
import numpy as np
from PIL import Image, ImageEnhance
import glob
import math
import os
from django.template.defaultfilters import slugify
import time

defaultSettings = photometry.Settings()
displayFile = ""
"""
class Thread(QThread):
    signal = QtCore.pyqtSignal(str, list)

    def __init__(self, m, da, bi, fl, i):
        QThread.__init__(self)
        self.mainArray = m
        self.darkArray = da
        self.biasArray = bi
        self.flatArray = fl
        self.i = i

    # run method gets called when we start the thread
    def run(self):
        print(self.i)
        time.sleep(self.i*5)
        # Display
        stars = 0
        displayAns = []
        # run photometry
        c  # self.walkthroughReferenceStars(mainArray[i].fileName, stars=stars)
        # Row 6 start
        #self.processinglbl.setText("Reference Star Selection")

        # for i in range(len(stars)):
        self.stars = stars
        self.file = self.mainArray[self.i].fileName
        #self.signal.emit(self.mainArray[self.i].fileName, stars)
        """

class HomePage(QtWidgets.QMainWindow):
    # Return Signals
    go_newProject = QtCore.pyqtSignal()

    def __init__(self):
        # Initialization
        super(HomePage, self).__init__()
        uic.loadUi('homepage/hompage.ui', self)

        # Buttons
        self.newprojectbtn.clicked.connect(lambda: self.goToNewProject())

    # Functions
    def goToNewProject(self):
        self.close()
        self.go_newProject.emit()

class NewProject(QtWidgets.QMainWindow):
    # Return Signals
    go_home = QtCore.pyqtSignal()
    go_processing = QtCore.pyqtSignal()

    def __init__(self):
        # Initialization
        super(NewProject, self).__init__() # Call the inherited classes __init__ method
        uic.loadUi('homepage/newproject.ui', self) # Load the .ui file

        # Message Boxes
        self.fileFailMessage = QMessageBox()
        self.fileFailMessage.setIcon(QMessageBox.Warning)
        self.fileFailMessage.setText("Invalid File/Folder, must be or include .fits files.")
        self.fileFailMessage.setWindowTitle("File Warning")
        self.fileFailMessage.setStandardButtons(QMessageBox.Ok)

        self.fileUploadMessage = QMessageBox()
        self.fileUploadMessage.setIcon(QMessageBox.Question)
        self.fileUploadMessage.setText("Will you be uploading more than one file?")
        self.fileUploadMessage.setWindowTitle("File Question")
        self.fileUploadMessage.setStandardButtons(QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)

        self.submitMessage = QMessageBox()
        self.submitMessage.setIcon(QMessageBox.Information)
        self.submitMessage.setText("These settings cannot be changed later. Submit anyways?")
        self.submitMessage.setWindowTitle("Are you sure?")
        self.submitMessage.setStandardButtons(QMessageBox.Yes | QMessageBox.Cancel)

        # Buttons
        self.mainFile = ""
        self.darkFile = ""
        self.biasFile = ""
        self.flatFile = ""
        self.referenceFile = ""
        self.homebtn.clicked.connect(lambda: self.goToHome())
        self.folderuploadbtn.clicked.connect(lambda: self.getMain("FOLDER"))
        self.fileuploadbtn.clicked.connect(lambda: self.getMain("FILE"))
        self.darkbtn.clicked.connect(self.getDark)
        self.biasbtn.clicked.connect(self.getBias)
        self.flatbtn.clicked.connect(self.getFlat)
        self.readinrefbtn.clicked.connect(self.getRef)
        self.submitbtn.clicked.connect(self.submit)

        # Info

        # Entry Initialization
        # astrometryDotNetFlag
        if not(defaultSettings.astrometryDotNetFlag == 0):
            self.apikeyentry.text = defaultSettings.astrometryDotNetFlag()
        # calibrationFlag
        if defaultSettings.calibrationFlag == 0:
            self.calibratecombo.setCurrentText("No")
        else:
            self.calibratecombo.setCurrentText("Yes")
        # subtractBiasFromDarkFlag
        if defaultSettings.subtractBiasFromDarkFlag == 0:
            self.subtractbiascombo.setCurrentText("No")
        else:
            self.subtractbiascombo.setCurrentText("Yes")
        # blankPerStarFlag
        if defaultSettings.blankPerStarFlag == 0:
            self.backgroundcombo.setCurrentText("From Entire Image")
        else:
            self.backgroundcombo.setCurrentText("From Individual Stars")
        # catalogChoice
        if defaultSettings.catalogChoice == 0:
            self.catalogentry.setText("SIMBAD")
        else:
            self.catalogentry.setText(defaultSettings.catalogChoice)
        # filterChoice
        if defaultSettings.filterChoice == "V":
            self.filtercombo.setCurrentText("Johnson V")
            self.filterentry.setDisabled(True)
        elif defaultSettings.filterChoice == "B":
            self.filtercombo.setCurrentText("Johnson B")
            self.filterentry.setDisabled(True)
        elif defaultSettings.filterChoice == "g":
            self.filtercombo.setCurrentText("g'")
            self.filterentry.setDisabled(True)
        elif defaultSettings.filterChoice == "i":
            self.filtercombo.setCurrentText("i'")
            self.filterentry.setDisabled(True)
        elif defaultSettings.filterChoice == "r":
            self.filtercombo.setCurrentText("r'")
            self.filterentry.setDisabled(True)
        else:
            self.filtercombo.setCurrentText("Other")
            self.filterentry.setText(str(defaultSettings.filterChoice))
            self.filterentry.setDisabled(False)
        self.filtercombo.currentTextChanged.connect(self.filterEntryDisable)
        # lightCurveLineFlag
        if defaultSettings.lightCurveLineFlag == 1:
            self.lightcurvelinecombo.setCurrentText("Yes")
        else:
            self.lightcurvelinecombo.setCurrentText("No")
        # errorChoice
        if defaultSettings.errorChoice == "STD":
            self.errorcombo.setCurrentText("Standard Deviation")
        elif defaultSettings.errorChoice == "JKF":
            self.errorcombo.setCurrentText("Jack Knife")
        elif defaultSettings.errorChoice == "WMG":
            self.errorcombo.setCurrentText("Weighted Magnitude")
        else:
            self.errorcombo.setCurrentText("No error calculation")
        # fwhmFlag
        if defaultSettings.fwhmFlag == 0:
            self.targetradiuscombo.setCurrentText("Find radius manually")
            self.pixradentry.setDisabled(True)
        elif defaultSettings.fwhmFlag == 1:
            self.targetradiuscombo.setCurrentText("Use full-width half-maximum")
            self.pixradentry.setDisabled(True)
        else:
            self.targetradiuscombo.setCurrentText("Other")
            self.pixradentry.setText(str(defaultSettings.fwhmFlag))
            self.pixradentry.setDisabled(False)
        self.targetradiuscombo.currentTextChanged.connect(self.fwhmEntryDisable)
        # readInRadiusFlag
        if defaultSettings.readInRadiusFlag == 1:
            self.refradiuscombo.setCurrentText("Use radius from file")
        else:
            self.refradiuscombo.setCurrentText("Same as target star radius")
        # removeReferenceOutliersFlag
        self.removeoutlierentry.setText(str(defaultSettings.removeReferenceOutliersFlag))
        # walkthroughMode
        if defaultSettings.walkthroughMode == 1:
            self.walkthroughcombo.setCurrentText("Yes")
        else:
            self.walkthroughcombo.setCurrentText("No")


    # Functions
    def goToHome(self):
        self.close()
        self.go_home.emit()

    def goToProcessing(self):
        self.close()
        self.go_processing.emit()

    def submit(self):
        global defaultSettings
        self.fileuploadbtn.setStyleSheet("")
        self.folderuploadbtn.setStyleSheet("")
        self.projectnameentry.setStyleSheet("")
        self.raentry.setStyleSheet("")
        self.decentry.setStyleSheet("")
        self.darkbtn.setStyleSheet("")
        self.biasbtn.setStyleSheet("")
        self.flatbtn.setStyleSheet("")
        self.apikeyentry.setStyleSheet("")
        self.catalogentry.setStyleSheet("")
        self.filterentry.setStyleSheet("")
        self.pixradentry.setStyleSheet("")
        self.removeoutlierentry.setStyleSheet("")
        errorMessage = self.check()
        if not(errorMessage == ""):
            checkMessage = QMessageBox()
            checkMessage.setIcon(QMessageBox.Warning)
            checkMessage.setText(errorMessage)
            checkMessage.setWindowTitle("Form Incomplete")
            checkMessage.setStandardButtons(QMessageBox.Ok)
            checkMessage.exec()
            return 0
        ans = self.submitMessage.exec()
        if ans == QMessageBox.Cancel:
            return 0
        elif ans == QMessageBox.Yes:
            self.setSettings()
            self.goToProcessing()

    def setSettings(self):
        # subtractBiasFromDarkFlag
        if self.subtractbiascombo.currentText() == "No":
            sbfdf = 0
        else:
            sbfdf = 1
        # calibrationFlag
        if self.calibratecombo.currentText() == "No":
            cal = 0
        else:
            cal = 1
        # blankPerStarFlag
        if self.backgroundcombo.currentText() == "From Entire Image":
            bpsf = 0
        else:
            bpsf = 1
        # catalogChoiceFlag
        if self.catalogentry.text() == "" or 0:
            cc = 0
        else:
            cc = str(self.catalogentry.text())
        # filterChoice
        if self.filtercombo.currentText() == "Johnson V":
            fc = "V"
        elif self.filtercombo.currentText() == "Johnson B":
            fc = "B"
        elif self.filtercombo.currentText() == "g'":
            fc = "g"
        elif self.filtercombo.currentText() == "r'":
            fc = "r"
        elif self.filtercombo.currentText() == "i'":
            fc = "i"
        elif self.filtercombo.currentText() == "Other":
            fc = str(self.filterentry.text())
        # lightCurveLineFlag
        if self.lightcurvelinecombo.currentText() == "No":
            lclf = 0
        else:
            lclf = 1
        # errorChoice
        if self.errorcombo.currentText() == "Standard Deviation":
            ec = "STD"
        elif self.errorcombo.currentText() == "Weighted Magnitude":
            ec = "WMG"
        elif self.errorcombo.currentText() == "Jack Knife":
            ec = "JKF"
        else:
            ec = 0
        # readInReferenceFlag
        if self.referenceFile == "":
            rirf = 0
        else:
            rirf = self.referenceFile
        # fwhmFlag
        if self.targetradiuscombo.currentText() == "Find radius manually":
            fwhm = 0
        elif self.targetradiuscombo.currentText() == "Use full-width half-maximum":
            fwhm = 1
        else:
            fwhm = int(self.pixradentry.text())
            if fwhm <= 1:
                fwhm = 2
        # astrometryDotNetFlag
        if self.apikeyentry.text() == "":
            adnf = 0
        else:
            adnf = self.apikeyentry.text()
        # readInRadiusFlag
        if not(rirf == 0) and self.refradiuscombo.currentText() == "Use radius from file":
            riradf = 1
        else:
            riradf = 0
        # useDarkFlag
        if self.darkFile == "" or cal == 0:
            udf = 0
        else:
            udf = 1
        # coordinateChoiceFlag
        if self.degreecombo.currentText() == "Decimal Degrees":
            ccf = "DEC"
        # useBiasFlag
        if self.biasFile == "" or cal == 0:
            ubf = 0
        else:
            ubf = 1
        # walkthroughMode
        if self.walkthroughcombo.currentText() == "Yes":
            wm = 1
        else:
            wm = 0
        photometry.changeSettings(subtractBiasFromDarkFlag=sbfdf, calibrationOutputFlag=1,
                                  calibrationFlag=cal, blankPerStarFlag=bpsf,
                                  catalogChoice=cc, filterChoice=fc,
                                  lightCurveLineFlag=lclf, showLightCurveFlag=-0,
                                  errorChoice=ec, consolePrintFlag=0,
                                  readInReferenceFlag=rirf, fwhmFlag=fwhm,
                                  printReferenceStarsFlag=1, astrometryDotNetFlag=adnf,
                                  removeReferenceOutliersFlag=int(self.removeoutlierentry.text()), readInRadiusFlag=-riradf,
                                  printLightCurveFlag=0, useDarkFlag=udf,
                                  astrometryTimeOutFlag=0, showCMDFlag=0,
                                  printCMDFlag=0, coordinateChoiceFlag=ccf,
                                  universalBlank=-1, useBiasFlag=ubf,
                                  projectName=self.projectnameentry.text(), rightAscension=float(self.raentry.text()),
                                  declination=float(self.decentry.text()), main=self.mainFile,
                                  dark=self.darkFile, bias=self.biasFile, flat=self.flatFile, interface=1,
                                  walkthroughMode=wm)

    def check(self):
        errorMessage = ""
        errorCount = 1
        # Check Main files
        if self.mainFile == "":
            errorMessage = errorMessage + "\n" + str(errorCount) + ") Main file required for photometry"
            errorCount = errorCount + 1
            self.fileuploadbtn.setStyleSheet("border: 1px solid red;")
            self.folderuploadbtn.setStyleSheet("border: 1px solid red;")
        # Check project name input
        name = self.projectnameentry.text()
        if name == "" or len(name) < 1:
            errorMessage = errorMessage + "\n" + str(errorCount) + ") Project Name cannot be blank"
            errorCount = errorCount + 1
            self.projectnameentry.setStyleSheet("border: 1px solid red;")
        else:
            path = os.path.isdir("PhotPSaveData/" + slugify(name))  # If path is true, need a different file name
            if path == True:
                errorMessage = errorMessage + "\n" + str(errorCount) + ") Project Name already exists, choose a unique project name"
                errorCount = errorCount + 1
        # Check coordinate input
        if self.degreecombo.currentText() == "Decimal Degrees":
            if self.raentry.text() == "" or len(self.raentry.text()) < 1:
                errorMessage = errorMessage + "\n" + str(errorCount) + ") Right Ascension cannot be blank"
                errorCount = errorCount + 1
                self.raentry.setStyleSheet("border: 1px solid red;")
            try:
                radius = float(self.raentry.text())
            except ValueError:
                errorMessage = errorMessage + "\n" + str(errorCount) + ") Right ascension must be a number."
                errorCount = errorCount + 1
                self.raentry.setStyleSheet("border: 1px solid red;")
            if self.decentry.text() == "" or len(self.decentry.text()) < 1:
                errorMessage = errorMessage + "\n" + str(errorCount) + ") Declination cannot be blank"
                errorCount = errorCount + 1
                self.decentry.setStyleSheet("border: 1px solid red;")
            try:
                radius = float(self.decentry.text())
            except ValueError:
                errorMessage = errorMessage + "\n" + str(errorCount) + ") Declination must be a number."
                errorCount = errorCount + 1
                self.decentry.setStyleSheet("border: 1px solid red;")

        """
        elif self.coordinateChoice == "Degrees":
            if self.degreeRA.get() == "" or len(self.degreeRA.get()) < 1:
                errorMessage = errorMessage + "\nRight Ascension degrees place cannot be blank"
            if self.minuteRA.get() == "" or len(self.minuteRA.get()) < 1:
                errorMessage = errorMessage + "\nRight Ascension minutes place cannot be blank"
            if self.secondRA.get() == "" or len(self.secondRA.get()) < 1:
                errorMessage = errorMessage + "\nRight Ascension seconds place cannot be blank"
            if self.degreeDec.get() == "" or len(self.degreeDec.get()) < 1:
                errorMessage = errorMessage + "\nDeclination degrees place cannot be blank"
            if self.minuteDec.get() == "" or len(self.minuteDec.get()) < 1:
                errorMessage = errorMessage + "\nDeclination minutes place cannot be blank"
            if self.secondDec.get() == "" or len(self.secondDec.get()) < 1:
                errorMessage = errorMessage + "\nDeclination seconds place cannot be blank"
        """
        # Check Calibration files
        if self.calibratecombo.currentText() == "Yes":
            if self.flatFile == "":
                errorMessage = errorMessage + "\n" + str(errorCount) + ") Flat field(s) required for calibration"
                errorCount = errorCount + 1
                self.flatbtn.setStyleSheet("border: 1px solid red;")
        if self.subtractbiascombo.currentText() == "Yes" and self.calibratecombo.currentText() == "Yes":
            if self.darkFile == "":
                errorMessage = errorMessage + "\n" + str(errorCount) + ") Select dark frame(s) OR set \"Subtract Bias from Dark?\" Setting to \"No\"."
                errorCount = errorCount + 1
                self.darkbtn.setStyleSheet("border: 1px solid red;")
            if self.biasFile == "":
                errorMessage = errorMessage + "\n" + str(errorCount) + ") Select bias frame(s) OR set \"Subtract Bias from Dark?\" Setting to \"No\"."
                errorCount = errorCount + 1
                self.biasbtn.setStyleSheet("border: 1px solid red;")
        # Check filter
        if self.filtercombo.currentText() == "Other":
            if self.filterentry.text() == "":
                errorMessage = errorMessage + "\n" + str(errorCount) + ") Select a filter from the list of filters or input custom filter."
                errorCount = errorCount + 1
                self.filterentry.setStyleSheet("border: 1px solid red;")
        # Check radius
        if self.targetradiuscombo.currentText() == "Other":
            if self.pixradentry.text() == "":
                errorMessage = errorMessage + "\n" + str(errorCount) + ") Select a target radius choice from the list or input a target radius."
                errorCount = errorCount + 1
                self.pixradentry.setStyleSheet("border: 1px solid red;")
            else:
                try:
                    radius = int(self.pixradentry.text())
                    if radius <=0:
                        errorMessage = errorMessage + "\n" + str(
                            errorCount) + ") Chosen target radius must be a positive integer."
                        errorCount = errorCount + 1
                        self.pixradentry.setStyleSheet("border: 1px solid red;")
                except ValueError:
                    errorMessage = errorMessage + "\n" + str(errorCount) + ") Chosen target radius must be an integer."
                    errorCount = errorCount + 1
                    self.pixradentry.setStyleSheet("border: 1px solid red;")
        # Check z-score
        try:
            if not(self.removeoutlierentry.text() == ""):
                radius = int(self.removeoutlierentry.text())
                if radius <= 0:
                    errorMessage = errorMessage + "\n" + str(
                        errorCount) + ") Z-score for outlier removal must be a positive integer."
                    errorCount = errorCount + 1
                    self.removeoutlierentry.setStyleSheet("border: 1px solid red;")
        except ValueError:
            errorMessage = errorMessage + "\n" + str(errorCount) + ") Z-score for outlier removal must be an integer."
            errorCount = errorCount + 1
            self.removeoutlierentry.setStyleSheet("border: 1px solid red;")
        return errorMessage

    def getMain(self, button):
        global displayFile
        if button == "FILE":
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                      "All Files (*);;FITS Files (*.fits);;FITS Files (*.fts)", options=options)
        elif button == "FOLDER":
            fileName = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        try:
            fileName[len(fileName) - 1]
        except IndexError:
            return 0
        if (fileName[len(fileName)-1] == 's' and fileName[len(fileName)-2] == 't'
                and fileName[len(fileName)-3] == 'i' and fileName[len(fileName)-4] == 'f'
                and fileName[len(fileName)-5] == '.') or (fileName[len(fileName)-1] == 's'
                and fileName[len(fileName)-2] == 't' and fileName[len(fileName)-3] == 'f'
                                                          and fileName[len(fileName)-4] == '.'):
            try:
                hdul = fits.open(fileName)  # hdul is the computer version of
                # the file data
            except OSError:
                self.fileFailMessage.exec()
                return 0
            #columns = shape[0]
            data = hdul[0].data
            image = Image.fromarray(data.astype('uint8'), 'L')
            enhancer = ImageEnhance.Contrast(image)
            img = enhancer.enhance(2)
            dat = img.tobytes("raw", "L")
            qImg = QImage(dat, data.shape[1], data.shape[0], data.shape[0], QImage.Format_Grayscale8)
            self.mainPixmap = QPixmap(qImg)
            self.imagelbl.setPixmap(self.mainPixmap)
            self.imagelbl.setScaledContents(True)
            self.filenamelbl.setText(fileName)
            self.mainFile = fileName
            displayFile = fileName
        else:
            displayFile = ""
            try:
                for filename in glob.glob(os.path.join(fileName, '*.fits')):
                    with open(os.path.join(os.getcwd(), filename), 'r') as f:
                        displayFile = filename
                for filename in glob.glob(os.path.join(fileName, '*.fts')):
                    with open(os.path.join(os.getcwd(), filename), 'r') as f:
                        displayFile = filename
            except:
                self.fileFailMessage.exec()
                return 0
            if displayFile == "":
                self.fileFailMessage.exec()
                return 0
            else:
                try:
                    hdul = fits.open(displayFile)  # hdul is the computer version of
                    # the file data
                except OSError:
                    self.fileFailMessage.exec()
                    return 0
                # columns = shape[0]
                data = hdul[0].data
                image = Image.fromarray(data.astype('uint8'), 'L')
                enhancer = ImageEnhance.Contrast(image)
                img = enhancer.enhance(2)
                dat = img.tobytes("raw", "L")
                qImg = QImage(dat, data.shape[1], data.shape[0], data.shape[0], QImage.Format_Grayscale8)
                self.mainPixmap = QPixmap(qImg)
                self.imagelbl.setPixmap(self.mainPixmap)
                self.imagelbl.setScaledContents(True)
                self.filenamelbl.setText(displayFile)
                self.mainFile = fileName

    def getDark(self):
        ans = self.fileUploadMessage.exec()
        if ans == QMessageBox.Yes:
            fileName = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        elif ans == QMessageBox.No:
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                      "All Files (*);;FITS Files (*.fits);;FITS Files (*.fts)",
                                                      options=options)
        elif ans == QMessageBox.Cancel:
            self.darkFile = ""
            self.darkbtn.setText("Dark Frame or Dark Frame Folder")
            return 0
        try:
            fileName[len(fileName) - 1]
        except IndexError:
            return 0
        if (fileName[len(fileName)-1] == 's' and fileName[len(fileName)-2] == 't'
                and fileName[len(fileName)-3] == 'i' and fileName[len(fileName)-4] == 'f'
                and fileName[len(fileName)-5] == '.') or (fileName[len(fileName)-1] == 's'
                and fileName[len(fileName)-2] == 't' and fileName[len(fileName)-3] == 'f'
                                                          and fileName[len(fileName)-4] == '.'):
            try:
                hdul = fits.open(fileName)  # hdul is the computer version of
                # the file data
            except OSError:
                self.fileFailMessage.exec()
                return 0
            self.darkbtn.setText(fileName)
            self.darkFile = fileName
        else:
            displayFile = ""
            try:
                for filename in glob.glob(os.path.join(fileName, '*.fits')):
                    with open(os.path.join(os.getcwd(), filename), 'r') as f:
                        displayFile = filename
                for filename in glob.glob(os.path.join(fileName, '*.fts')):
                    with open(os.path.join(os.getcwd(), filename), 'r') as f:
                        displayFile = filename
            except:
                self.fileFailMessage.exec()
                return 0
            if displayFile == "":
                self.fileFailMessage.exec()
                return 0
            else:
                self.darkbtn.setText(fileName)
                self.darkFile = fileName

    def getBias(self):
        ans = self.fileUploadMessage.exec()
        if ans == QMessageBox.Yes:
            fileName = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        elif ans == QMessageBox.No:
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                      "All Files (*);;FITS Files (*.fits);;FITS Files (*.fts)",
                                                      options=options)
        elif ans == QMessageBox.Cancel:
            self.darkFile = ""
            self.darkbtn.setText("Bias Frame or Bias Frame Folder")
            return 0
        try:
            fileName[len(fileName) - 1]
        except IndexError:
            return 0
        if (fileName[len(fileName)-1] == 's' and fileName[len(fileName)-2] == 't'
                and fileName[len(fileName)-3] == 'i' and fileName[len(fileName)-4] == 'f'
                and fileName[len(fileName)-5] == '.') or (fileName[len(fileName)-1] == 's'
                and fileName[len(fileName)-2] == 't' and fileName[len(fileName)-3] == 'f'
                                                          and fileName[len(fileName)-4] == '.'):
            try:
                hdul = fits.open(fileName)  # hdul is the computer version of
                # the file data
            except OSError:
                self.fileFailMessage.exec()
                return 0
            self.biasbtn.setText(fileName)
            self.biasFile = fileName
        else:
            displayFile = ""
            try:
                for filename in glob.glob(os.path.join(fileName, '*.fits')):
                    with open(os.path.join(os.getcwd(), filename), 'r') as f:
                        displayFile = filename
                for filename in glob.glob(os.path.join(fileName, '*.fts')):
                    with open(os.path.join(os.getcwd(), filename), 'r') as f:
                        displayFile = filename
            except:
                self.fileFailMessage.exec()
                return 0
            if displayFile == "":
                self.fileFailMessage.exec()
                return 0
            else:
                self.biasbtn.setText(fileName)
                self.biasFile = fileName

    def getFlat(self):
        ans = self.fileUploadMessage.exec()
        if ans == QMessageBox.Yes:
            fileName = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        elif ans == QMessageBox.No:
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                      "All Files (*);;FITS Files (*.fits);;FITS Files (*.fts)",
                                                      options=options)
        elif ans == QMessageBox.Cancel:
            self.flatFile = ""
            self.flatbtn.setText("Flat Field or Flat Field Folder")
            return 0
        try:
            fileName[len(fileName) - 1]
        except IndexError:
            return 0
        if (fileName[len(fileName)-1] == 's' and fileName[len(fileName)-2] == 't'
                and fileName[len(fileName)-3] == 'i' and fileName[len(fileName)-4] == 'f'
                and fileName[len(fileName)-5] == '.') or (fileName[len(fileName)-1] == 's'
                and fileName[len(fileName)-2] == 't' and fileName[len(fileName)-3] == 'f'
                                                          and fileName[len(fileName)-4] == '.'):
            try:
                hdul = fits.open(fileName)  # hdul is the computer version of
                # the file data
            except OSError:
                self.fileFailMessage.exec()
                return 0
            self.flatbtn.setText(fileName)
            self.flatFile = fileName
        else:
            displayFile = ""
            try:
                for filename in glob.glob(os.path.join(fileName, '*.fits')):
                    with open(os.path.join(os.getcwd(), filename), 'r') as f:
                        displayFile = filename
                for filename in glob.glob(os.path.join(fileName, '*.fts')):
                    with open(os.path.join(os.getcwd(), filename), 'r') as f:
                        displayFile = filename
            except:
                self.fileFailMessage.exec()
                return 0
            if displayFile == "":
                self.fileFailMessage.exec()
                return 0
            else:
                self.flatbtn.setText(fileName)
                self.flatFile = fileName

    def getRef(self):
        message = QMessageBox()
        message.setIcon(QMessageBox.Question)
        message.setText("Change reference read in file?")
        message.setWindowTitle("File Question")
        message.setStandardButtons(QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
        if self.referenceFile == "":
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            fileName, _ = QFileDialog.getOpenFileName(self, "Select .csv file", "",
                                                      "All Files (*);;CSV Files (*.csv)",
                                                      options=options)
        else:
            ans = message.exec()
            if ans == QMessageBox.Yes:
                options = QFileDialog.Options()
                options |= QFileDialog.DontUseNativeDialog
                fileName, _ = QFileDialog.getOpenFileName(self, "Select .csv file", "",
                                                          "All Files (*);;CSV Files (*.csv)",
                                                          options=options)
            elif ans == QMessageBox.No:
                return 0
            elif ans == QMessageBox.Cancel:
                self.referenceFile = ""
                self.readinrefbtn.setText("Reference Star Read In File (.csv)")
                return 0
        try:
            fileName[len(fileName) - 1]
        except IndexError:
            return 0
        if fileName[len(fileName) - 1] == 'v' and fileName[len(fileName) - 2] == 's' \
                and fileName[len(fileName) - 3] == 'c' and fileName[len(fileName) - 4] == '.':
            try:
                file = open(fileName, "r")
            except FileNotFoundError:
                self.fileFailMessage.exec()
                return 0
            self.referenceFile = fileName
            self.readinrefbtn.setText(fileName)

    def filterEntryDisable(self):
        if self.filtercombo.currentText() == "Other":
            self.filterentry.setDisabled(False)
        else:
            self.filterentry.setDisabled(True)

    def fwhmEntryDisable(self):
        if self.targetradiuscombo.currentText() == "Other":
            self.pixradentry.setDisabled(False)
        else:
            self.pixradentry.setDisabled(True)

class Processing(QtWidgets.QMainWindow):
    # Return Signals
    go_home = QtCore.pyqtSignal()

    def __init__(self):
        # Initialization
        super(Processing, self).__init__() # Call the inherited classes __init__ method
        uic.loadUi('homepage/processing.ui', self) # Load the .ui file
        self.show()
        self.projectnameentry.setText(photometry.getprojectName())
        self.raentry.setText(str(photometry.getrightAscension()))
        self.decentry.setText(str(photometry.getdeclination()))
        self.answers = []

        # Page setup
        self.homebtn.clicked.connect(lambda: self.goToHome())

        # Dialog box
        self.proceedMessage = QMessageBox()
        self.proceedMessage.setIcon(QMessageBox.Information)
        self.proceedMessage.setText("Proceed with these reference stars?")
        self.proceedMessage.setWindowTitle("Submit Reference Stars")
        self.proceedMessage.setStandardButtons(QMessageBox.Yes | QMessageBox.Cancel)

        self.exitWalkthroughMessage = QMessageBox()
        self.exitWalkthroughMessage.setIcon(QMessageBox.Information)
        self.exitWalkthroughMessage.setText("Would you like to exit walkthrough mode and use these reference stars as the default?\n")
        self.exitWalkthroughMessage.setWindowTitle("Exit Walkthrough")
        self.exitWalkthroughMessage.setStandardButtons(QMessageBox.Yes | QMessageBox.No)

        # Display
        hdul = fits.open(displayFile)
        data = hdul[0].data
        image = Image.fromarray(data.astype('uint8'), 'L')
        enhancer = ImageEnhance.Contrast(image)
        img = enhancer.enhance(2)
        dat = img.tobytes("raw", "L")
        qImg = QImage(dat, data.shape[1], data.shape[0], data.shape[0], QImage.Format_Grayscale8)
        self.mainPixmap = QPixmap(qImg)
        self.imagelbl.setPixmap(self.mainPixmap)
        self.imagelbl.setScaledContents(True)
        self.filenamelbl.setText(displayFile)
        #self.run.emit()
        QApplication.processEvents()
        self.runPhotometry()
        QApplication.processEvents()
        while self.answers == []:
            time.sleep(1)
        photometry.plotResults(self.answers, chartTitle=photometry.getprojectName())
        photometry.printResultsToFile(self.answers)

    # Functions
    def goToHome(self):
        self.close()
        self.go_home.emit()

    def runPhotometry(self):
        # Check if main is a single file or a directory
        if displayFile == photometry.getmain():
            mainArray = []
            temp = photometry.Photometry(photometry.getmain(), 0, 0, 0, 0)
            mainArray.append(temp)
        else:
            mainArray = photometry.getFiles(photometry.getmain())
        # Check if calibration files are single or in a directory
        # dark
        if photometry.getuseDarkFlag() == 0 or photometry.getcalibrationFlag() == 0:
            darkArray = []
        elif (photometry.getdark()[len(photometry.getdark()) - 1] == 's' and photometry.getdark()[
            len(photometry.getdark()) - 2] == 't'
              and photometry.getdark()[len(photometry.getdark()) - 3] == 'i' and photometry.getdark()[
                  len(photometry.getdark()) - 4] == 'f'
              and photometry.getdark()[len(photometry.getdark()) - 5] == '.') or (
                photometry.getdark()[len(photometry.getdark()) - 1] == 's'
                and photometry.getdark()[len(photometry.getdark()) - 2] == 't' and photometry.getdark()[
                    len(photometry.getdark()) - 3] == 'f'
                and photometry.getdark()[len(photometry.getdark()) - 4] == '.'):
            darkArray = []
            temp = photometry.Photometry(photometry.getdark(), 0, 0, 0, 0)
            darkArray.append(temp)
        else:
            darkArray = photometry.getFiles(photometry.getdark())
        # bias
        if photometry.getuseBiasFlag() == 0 or photometry.getcalibrationFlag() == 0:
            biasArray = []
        elif (photometry.getbias()[len(photometry.getbias()) - 1] == 's' and photometry.getbias()[
            len(photometry.getbias()) - 2] == 't'
              and photometry.getbias()[len(photometry.getbias()) - 3] == 'i' and photometry.getbias()[
                  len(photometry.getbias()) - 4] == 'f'
              and photometry.getbias()[len(photometry.getbias()) - 5] == '.') or (
                photometry.getbias()[len(photometry.getbias()) - 1] == 's'
                and photometry.getbias()[len(photometry.getbias()) - 2] == 't' and photometry.getbias()[
                    len(photometry.getbias()) - 3] == 'f'
                and photometry.getbias()[len(photometry.getbias()) - 4] == '.'):
            biasArray = []
            temp = photometry.Photometry(photometry.getbias(), 0, 0, 0, 0)
            biasArray.append(temp)
        else:
            biasArray = photometry.getFiles(photometry.getbias())
        # flat
        if photometry.getcalibrationFlag() == 0:
            flatArray = []
        elif (photometry.getflat()[len(photometry.getflat()) - 1] == 's' and photometry.getflat()[
            len(photometry.getflat()) - 2] == 't'
              and photometry.getflat()[len(photometry.getflat()) - 3] == 'i' and photometry.getflat()[
                  len(photometry.getflat()) - 4] == 'f'
              and photometry.getflat()[len(photometry.getflat()) - 5] == '.') or (
                photometry.getflat()[len(photometry.getflat()) - 1] == 's'
                and photometry.getflat()[len(photometry.getflat()) - 2] == 't' and photometry.getflat()[
                    len(photometry.getflat()) - 3] == 'f'
                and photometry.getflat()[len(photometry.getflat()) - 4] == '.'):
            flatArray = []
            temp = photometry.Photometry(photometry.getflat(), 0, 0, 0, 0)
            flatArray.append(temp)
        else:
            flatArray = photometry.getFiles(photometry.getflat())
        ans = 0
        for i in range(len(mainArray)):
            ans = self.runLetsGo(mainArray, darkArray, biasArray, flatArray, i)
            if not(ans == 0):
                self.answers.append(ans)
                # Check reference stars
                self.drawCircles(ans.fileName, ans.referenceStars)
                QApplication.processEvents()
                if photometry.getwalkthroughMode() == 1:
                    self.processinglbl.setText("Reference Star Options")
                    self.refChecks = []
                    self.refLbls = []
                    for j in range(len(ans.referenceStars)):
                        tempCheck = QCheckBox("\n" + str(ans.referenceStars[j].id) + " Magnitude:" + str(photometry.truncate(ans.referenceStars[j].magnitude,decimals=2))
                                              + " \nCalculated Target Magnitude: " + str(photometry.truncate(ans.referenceStars[j].targetMagnitude, decimals=2))
                                              + "\n")
                        tempCheck.setChecked(True)
                        tempCheck.stateChanged.connect(lambda: self.starUnchecked(ans))
                        self.gridLayout.addWidget(tempCheck, 7+j, 1, 2, 2)
                        QApplication.processEvents()
                        self.refChecks.append(tempCheck)
                        QApplication.processEvents()
                    self.submitBtn = QPushButton("Submit Reference Stars")
                    self.gridLayout.addWidget(self.submitBtn, 8+len(ans.referenceStars), 1, 1, 2)
                    QApplication.processEvents()
                    self.loop = QtCore.QEventLoop()
                    self.submitBtn.clicked.connect(lambda: self.submitReferenceStars(ans, mainArray, darkArray, biasArray, flatArray, i))
                    self.loop.exec_()
                    #time.sleep(60)

        return ans

    def starUnchecked(self, ans):
        tempStars = []
        for i in range(len(self.refChecks)):
            if self.refChecks[i].isChecked() == True:
                tempStars.append(ans.referenceStars[i])
        self.drawCircles(ans.fileName, tempStars)

    def submitReferenceStars(self, ans, mainArray, darkArray, biasArray, flatArray, i):
        response = self.proceedMessage.exec()
        if response == QMessageBox.Cancel:
            return 0
        elif response == QMessageBox.Yes:
            x = len(self.refChecks) - 1
            while x >= 0:
                self.refChecks[x].hide()
                self.submitBtn.hide()
                QApplication.processEvents()
                #del self.refChecks[x]
                x = x - 1
                QApplication.processEvents()
            self.processinglbl.setText("Processing \n(Please do not\nclose this window)")
            QApplication.processEvents()
            tempStars = []
            for a in range(len(self.refChecks)):
                if self.refChecks[a].isChecked() == True:
                    tempStars.append(ans.referenceStars[a])
            refSetting = photometry.getreadInReferenceFlag()
            photometry.changeSettings(readInReferenceFlag=tempStars)
            ans = self.runLetsGo(mainArray, darkArray, biasArray, flatArray, i)
            if not(ans == 0):
                for b in range(len(self.answers)):
                    if self.answers[b].fileName == ans.fileName:
                        self.answers[b] = ans
            response2 = self.exitWalkthroughMessage.exec()
            if response2 == QMessageBox.No:
                photometry.changeSettings(readInReferenceFlag=refSetting)
            elif response2 == QMessageBox.Yes:
                photometry.changeSettings(walkthroughMode=0)
        self.loop.quit()

    def drawCircles(self, file, stars):
        hdul = fits.open(file)  # hdul is the computer version of
        data = hdul[0].data
        image = Image.fromarray(data.astype('uint8'), 'L')
        enhancer = ImageEnhance.Contrast(image)
        img = enhancer.enhance(2)
        dat = img.tobytes("raw", "L")
        qImg = QImage(dat, data.shape[1], data.shape[0], data.shape[0], QImage.Format_Grayscale8)
        self.mainPixmap = QPixmap(qImg)
        w = photometry.worldCoordinateSystem(hdul)
        tarX, tarY = w.all_world2pix(photometry.getrightAscension(), photometry.getdeclination(), 0)

        painter = QPainter()
        painter.begin(self.mainPixmap)
        pen = QPen()
        pen.setColor(QColor(200, 0, 0))
        pen.setWidth(10)
        painter.setPen(pen)

        for i in range(len(stars)):
            X, Y = w.all_world2pix(stars[i].ra, stars[i].dec, 0)
            if not((X < tarX+10) and (X > tarX-10)) and not((Y < tarY+10) and (Y > tarY-10)):
                painter.drawEllipse(X-50, Y-50, 100, 100)
                painter.setFont(QFont('Helvetica', 64))
                painter.drawText(X-100, Y-70, str(stars[i].id))

        pen.setColor(QColor(0, 200, 0))
        painter.setPen(pen)
        painter.drawEllipse(tarX-100, tarY-100, 200, 200)
        painter.drawText(tarX-164, tarY - 120, "Target Star")
        painter.end()
        self.imagelbl.setPixmap(self.mainPixmap)
        self.imagelbl.setScaledContents(True)
        self.filenamelbl.setText(file)
        hdul.close()
        QApplication.processEvents()

    def runLetsGo(self, mainArray, darkArray, biasArray, flatArray, i):
        if photometry.getcalibrationFlag() == 0:
            ans = photometry.letsGo(photometry.getrightAscension(), photometry.getdeclination(),
                                    mainArray[i].fileName)
        elif photometry.getcalibrationFlag() == 1 and photometry.getuseDarkFlag() == 1 and photometry.getuseBiasFlag() == 1:
            d = photometry.matchCal(mainArray[i].fileName, darkArray)
            b = photometry.matchCal(mainArray[i].fileName, biasArray)
            f = photometry.matchCal(mainArray[i].fileName, flatArray)
            ans = photometry.letsGo(photometry.getrightAscension(), photometry.getdeclination(),
                                    mainArray[i].fileName, d, b, f)
        elif photometry.getcalibrationFlag() == 1 and photometry.getuseDarkFlag() == 0 and photometry.getuseBiasFlag() == 1:
            b = photometry.matchCal(mainArray[i].fileName, biasArray)
            f = photometry.matchCal(mainArray[i].fileName, flatArray)
            ans = photometry.letsGo(photometry.getrightAscension(), photometry.getdeclination(),
                                    mainArray[i].fileName, biasFrame=b, flatField=f)
        elif photometry.getcalibrationFlag() == 1 and photometry.getuseDarkFlag() == 1 and photometry.getuseBiasFlag() == 0:
            d = photometry.matchCal(mainArray[i].fileName, darkArray)
            f = photometry.matchCal(mainArray[i].fileName, flatArray)
            ans = photometry.letsGo(photometry.getrightAscension(), photometry.getdeclination(),
                                    mainArray[i].fileName, darkFrame=d, flatField=f)
        else:
            f = photometry.matchCal(mainArray[i].fileName, flatArray)
            ans = photometry.letsGo(photometry.getrightAscension(), photometry.getdeclination(),
                                    mainArray[i].fileName, flatField=f)
        return ans

class Controller:

    def __init__(self):
        pass

    def show_home(self):
        self.homePage = HomePage()
        self.homePage.go_newProject.connect(self.show_new_project)
        self.homePage.show()

    def show_new_project(self):
        self.newProject = NewProject()
        self.newProject.go_home.connect(self.show_home)
        self.newProject.go_processing.connect(self.show_processing)
        self.newProject.show()

    def show_processing(self):
        self.processing = Processing()
        self.processing.go_home.connect(self.show_home)
        #self.processing.runPhotometry()
