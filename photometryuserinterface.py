"""
Program written by Alexis Tudor at the University of Nevada, Reno
Email at alexisrenee1@gmail.com
Copyright and Licensing: GNU @ Alexis Tudor
"""
import photometry
from PyQt5 import QtWidgets, uic, QtCore
from PyQt5.QtCore import QUrl
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QApplication, QPushButton, QLabel, QCheckBox, QFrame
from PyQt5.QtGui import QPixmap, QImage, QPainter, QPen, QColor, QFont, QDesktopServices
from astropy.io import fits
from PIL import Image, ImageEnhance
import glob
import math
import os
from django.template.defaultfilters import slugify
import time
import shutil
from shutil import copyfile
import warnings

defaultSettings = photometry.Settings()
displayFile = ""
currProject = ""

class HomePage(QtWidgets.QMainWindow):
    # Return Signals
    go_home = QtCore.pyqtSignal()
    go_my_projects = QtCore.pyqtSignal()
    go_settings = QtCore.pyqtSignal()
    go_about = QtCore.pyqtSignal()
    go_processing = QtCore.pyqtSignal()
    go_help = QtCore.pyqtSignal()
    go_project_page = QtCore.pyqtSignal()
    go_newProject = QtCore.pyqtSignal()

    def __init__(self):
        # Initialization
        super(HomePage, self).__init__()
        uic.loadUi('homepage/hompage.ui', self)

        # Buttons
        self.newprojectbtn.clicked.connect(lambda: self.goToNewProject())
        self.myprojectsbtn.clicked.connect(lambda: self.goToMyProjects())
        self.aboutbtn.clicked.connect(lambda: self.goToAbout())
        self.settingsbtn.clicked.connect(lambda: self.goToSettings())
        self.helpbtn.clicked.connect(lambda: self.goToHelp())

    # Functions
    # Functions to go to different screens
    def goToHome(self):
        self.close()
        self.go_home.emit()

    def goToNewProject(self):
        self.close()
        self.go_newProject.emit()

    def goToProcessing(self):
        self.close()
        self.go_processing.emit()

    def goToSettings(self):
        self.close()
        self.go_settings.emit()

    def goToAbout(self):
        self.close()
        self.go_about.emit()

    def goToHelp(self):
        self.close()
        self.go_help.emit()

    def goToMyProjects(self):
        self.close()
        self.go_my_projects.emit()

    def goToProjectPage(self):
        self.close()
        self.go_project_page.emit()

class NewProject(QtWidgets.QMainWindow):
    # Return Signals
    go_home = QtCore.pyqtSignal()
    go_my_projects = QtCore.pyqtSignal()
    go_settings = QtCore.pyqtSignal()
    go_about = QtCore.pyqtSignal()
    go_processing = QtCore.pyqtSignal()
    go_help = QtCore.pyqtSignal()
    go_project_page = QtCore.pyqtSignal()

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
        self.fileUploadMessage.setText("Will you be uploading more than one file? Select \"Yes\" to upload a folder, select \"No\" to upload a single file.")
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
        self.errorFile = ""
        self.referenceFile = ""
        self.calibrationFileDirectory = ""
        self.homebtn.clicked.connect(lambda: self.goToHome())
        self.myprojectsbtn.clicked.connect(lambda: self.goToMyProjects())
        self.aboutbtn.clicked.connect(lambda: self.goToAbout())
        self.settingsbtn.clicked.connect(lambda: self.goToSettings())
        self.helpbtn.clicked.connect(lambda: self.goToHelp())
        self.folderuploadbtn.clicked.connect(lambda: self.getMain("FOLDER"))
        self.fileuploadbtn.clicked.connect(lambda: self.getMain("FILE"))
        self.darkbtn.clicked.connect(self.getDark)
        self.biasbtn.clicked.connect(self.getBias)
        self.flatbtn.clicked.connect(self.getFlat)
        self.readinrefbtn.clicked.connect(self.getRef)
        self.errorstarbtn.clicked.connect(self.getError)
        self.submitbtn.clicked.connect(self.submit)
        self.calFileButton.clicked.connect(self.getCalibrationDirectory)

        # Info buttons setup
        self.fileuploadbtn.setToolTip("Flexible Image Transport System (FITS) files contain\n"+
                                      "both image and header data that can be used not only to\n"+
                                      "see images of the sky but also to keep track of what part\n"+
                                      "of the sky they represent. Use this button to upload a\n"+
                                      "single FITS file for calculation.")

        self.folderuploadbtn.setToolTip("Flexible Image Transport System (FITS) files contain\n" +
                                      "both image and header data that can be used not only to\n" +
                                      "see images of the sky but also to keep track of what part\n" +
                                      "of the sky they represent. Use this button to upload a\n" +
                                      "folder containing FITS files for calculation. The folder\n"+
                                        "does not need to contain only FITS files, but only .fits\n"+
                                        "or .fts files will be considered.")

        self.projectnameinfo.setToolTip("Project names cannot match \nany previous project names.")

        self.projectnameinfo.mousePressEvent = self.projectnameMethod

        self.rainfo.setToolTip("Right Ascension is the angular distance\n"
                               + "of a point east along the celestial equator.\n"+
                               "It is similar to longitude on Earth. It is usually\n"
                               "the stellar coordinate listed first.\n\n"+
                               "Decimal Degree Example: 175.91034694\n" +
                               "Sexagesimal Example: 11 43 38.4924789768")

        self.rainfo.mousePressEvent = self.raMethod

        self.decinfo.setToolTip("Declination is the angular distance\n"+
                                "of a point north or south of the celestial equator.\n"+
                                "It is similar to latitude on Earth. It is usually\n"
                                +"the stellar coordinate listed second.\n\n"+
                               "Decimal Degree Example: 71.68911273\n" +
                               "Sexagesimal Example: +71 41 20.559348288")

        self.decinfo.mousePressEvent = self.decMethod

        self.darkbtn.setToolTip("Dark frames are subtracted from a telescope\n"+
                                "image to remove thermal noise. If multiple dark\n"+
                                "frames are uploaded, for each image a frame will be\n"+
                                "selected based on date and exposure time.")
        self.biasbtn.setToolTip("Bias frames are subtracted from a telescope\n" +
                                "image to remove noise that comes from the amplification\n"+
                                "process the images undergo in processing. If multiple bias\n" +
                                "frames are uploaded, for each image a frame will be\n" +
                                "selected based on date and exposure time.")
        self.flatbtn.setToolTip("Flat fields are divided from a telescope image\n" +
                                "to remove noise that comes from sensitivity and gain.\n" +
                                "If multiple flat fields are uploaded, for each image\n" +
                                "a frame will be selected based on date and exposure time.")

        self.astrometryinfo.setToolTip("Astrometry.net is a service to add world coordinate\n"+
                                       "system (wcs) information to an image of the sky. The\n"+
                                       "wcs data allows the program to know where in the sky\n"+
                                       "an image was taken, and locate stars by stellar coordinates.")

        self.astrometryinfo.mousePressEvent = self.astrometryMethod

        self.calibrateinfo.setToolTip("Calibrating telescope images removes noise introduced by\n"+
                                      "heat, amplification, and gain. Only flat fields are required\n"+
                                      "for calibration.")

        self.calibrateinfo.mousePressEvent = self.calibrateMethod

        self.subtractbiasinfo.setToolTip("Some telescopes include the bias frame information in\n"+
                                         "the dark frame as well. In this case, the bias must be\n"+
                                         "subtracted from the dark frame before calibration to avoid\n"+
                                         "subtracting the bias out twice. (For the Great Basin Observatory\n"+
                                         "please set this to yes if using both dark and bias frames)")

        self.subtractbiasinfo.mousePressEvent = self.subtractbiasMethod

        self.backgroundinfo.setToolTip("In telescope images there is always some noise from background\n"+
                                       "radiation. This background radiation must be subtracted out.\n"+
                                       "Select \"From Entire Image\" to take an count of the\n"+
                                       "background radiation in the entirety of the image, or \n"+
                                       "\"From Individual Stars\" to calculate a local background\n"+
                                       "count for every individual star (recommended when there may\n"+
                                       "be a gradient in the image).")

        self.backgroundinfo.mousePressEvent = self.backgroundMethod

        self.cataloginfo.setToolTip("Photometry+ uses Vizier catalogs to find magnitudes for reference\n"+
                                    "stars. Enter the catalog you'd like to use to search for reference star\n"+
                                    "magnitudes. To use SIMBAD enter \"SIMBAD\". If using uploaded\n"+
                                    "reference stars, feel free to leave this entry blank.")

        self.cataloginfo.mousePressEvent = self.catalogMethod

        self.filterinfo.setToolTip("The filter choice is used to find correct magnitudes for\n"+
                                   "reference stars. Please select the filter the telescope image\n"
                                   "was taken in.")

        self.filterinfo.mousePressEvent = self.filterMethod

        self.otherfilterinfo.setToolTip("If your filter is not supported by Photometry+ right now,\n"
                                        +"please enter the filter you would like to use.\n"+
                                        "WARNING: Filter formatting must be correct for catalog\n"+
                                        "being used. Please double check entered filter formatting.")

        self.otherfilterinfo.mousePressEvent = self.otherfilterMethod

        self.lightcurvelineinfo.setToolTip("Select \"Yes\" to add a light curve connecting magnitudes on\n"+
                                           "final graph. Traditional light curves typically do not have lines.")

        self.lightcurvelineinfo.mousePressEvent = self.lclineMethod

        self.magnitudeinfo.mousePressEvent = self.magnitudeMethod

        self.errorinfo.setToolTip("Select how errors for calculations are calculated.\n"+
                                  "\"Standard Deviation\" takes the standard deviation of\n"+
                                  "all reference star magnitude calculations before averaging\n"+
                                  "them together and uses that value as error.\n"+
                                  "\"Weighted Magnitude\" calculates the magnitude based on\n"+
                                  "the photons counted in the telescope image.\n"+
                                  "\"Jack Knife\" calculates error based on the jack knife method\n"+
                                  "by Anderson et al. in \"A Simple and Direct Measure of Photometric Uncertainties\"\n"+
                                  "Or you can opt to not calculate error (not recommended).")

        self.errorinfo.mousePressEvent = self.errorMethod

        self.readinrefbtn.setToolTip("To calculate the brightness of a star, it must be\n"+
                                     "compared to several other stars of known magnitude.\n"+
                                     "Photometry+ will automatically detect such reference stars,\n"+
                                     "or you can upload your reference stars in the form of a .csv file.\n"+
                                     "See Photometry+ documentation for the correct format for the .csv upload.")


        self.targetradiusinfo.setToolTip("Select how the radius of the target object is determined.\n"+
                                         "Use \"Find radius manually\" to have the radius of the star\n"+
                                         "determined by looking at the edge of the star.\n"+
                                         "Use \"Use full-width half-maximum\" to use three times the \n"+
                                         "full-width half-maximum (FWHM) as the radius. The  FWHM represents\n"+
                                         "the width of that star area where the light coming from the star is\n"+
                                         "at half the maximum light that star puts out.\n"+
                                         "You may also put in the radius to use manually.")

        self.targetradiusinfo.mousePressEvent = self.radiusMethod

        self.pixradiusinfo.setToolTip("Enter the radius, in pixels, to use as the target object radius.")

        self.pixradiusinfo.mousePressEvent = self.pixradiusMethod

        self.refradiusinfo.setToolTip("Select how the radius of reference stars are determined.\n"+
                                      "Set\"Same as target object radius\" to use the target object radius\n"+
                                      "for all stars. This will reduce error introduced by differences\n"+
                                      "in the magnitude of reference stars."+
                                      "Set \"Use radius from file\" to use the radius measurements\n"+
                                      "found in uploaded reference star file.")

        self.refradiusinfo.mousePressEvent = self.refradiusMethod

        self.removerefinfo.setToolTip("After reference stars are used to calculate the magnitude\n"+
                                      "of the target object, those magnitudes are assigned a z-score.\n"+
                                      "The z-score of a number is the measurement of the number of\n"+
                                      "standard deviations a number is away from the mean. Removing\n"+
                                      "reference stars more than 3 standard deviations away from\n"+
                                      "the mean can reduce noise caused by flawed reference stars.")

        self.removerefinfo.mousePressEvent = self.outlierMethod

        self.walkthroughinfo.setToolTip("Use walkthrough mode to step through the program one step\n"+
                                        "at a time for more control over the program and a look into\n"+
                                        "what is actually happening inside of it.")

        self.walkthroughinfo.mousePressEvent = self.walkthroughMethod

        # Entry Initialization
        # astrometryDotNetFlag
        if not(defaultSettings.astrometryDotNetFlag == 0 or defaultSettings.astrometryDotNetFlag == "0"):
            self.apikeyentry.setText(defaultSettings.astrometryDotNetFlag)
        # calibrationFlag
        if defaultSettings.calibrationFlag == 0 or defaultSettings.calibrationFlag == "0":
            self.calibratecombo.setCurrentText("No")
        else:
            self.calibratecombo.setCurrentText("Yes")
        # subtractBiasFromDarkFlag
        if defaultSettings.subtractBiasFromDarkFlag == 0 or defaultSettings.subtractBiasFromDarkFlag == "0":
            self.subtractbiascombo.setCurrentText("No")
        else:
            self.subtractbiascombo.setCurrentText("Yes")
        # blankPerStarFlag
        if defaultSettings.blankPerStarFlag == 0 or defaultSettings.blankPerStarFlag == "0":
            self.backgroundcombo.setCurrentText("From Entire Image")
        else:
            self.backgroundcombo.setCurrentText("From Individual Stars")
        # catalogChoice
        if defaultSettings.catalogChoice == 0 or defaultSettings.catalogChoice == "0":
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
        # oneMagnitudeCalculation
        if defaultSettings.oneMagnitudeCalculation == 0:
            self.magnitudecombo.setCurrentText("Multiple Calculations")
        else:
            self.magnitudecombo.setCurrentText("One Calculation")
        # errorChoice
        if defaultSettings.errorChoice == "STD":
            self.errorcombo.setCurrentText("Standard Deviation")
            self.errorstarbtn.setDisabled(True)
        elif defaultSettings.errorChoice == "JKF":
            self.errorcombo.setCurrentText("Jack Knife")
            self.errorstarbtn.setDisabled(True)
        elif defaultSettings.errorChoice == "WMG":
            self.errorcombo.setCurrentText("Weighted Magnitude")
            self.errorstarbtn.setDisabled(True)
        else:
            self.errorcombo.setCurrentText("No error calculation")
            self.errorstarbtn.setDisabled(True)
        self.errorcombo.currentTextChanged.connect(self.errorDisable)
        # fwhmFlag
        if defaultSettings.fwhmFlag == 0 or defaultSettings.fwhmFlag == "0":
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
    # Return functions
    def goToHome(self):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Warning)
        self.message.setText("Warning, leaving this page will delete all unsaved information for this new project.")
        self.message.setWindowTitle("Warning")
        self.message.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        ans = self.message.exec()
        if ans == QMessageBox.No:
            return 0
        else:
            self.close()
            self.go_home.emit()

    def goToProcessing(self):
        self.close()
        self.go_processing.emit()

    def goToSettings(self):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Warning)
        self.message.setText("Warning, leaving this page will delete all unsaved information for this new project.")
        self.message.setWindowTitle("Warning")
        self.message.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        ans = self.message.exec()
        if ans == QMessageBox.No:
            return 0
        else:
            self.close()
            self.go_settings.emit()

    def goToAbout(self):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Warning)
        self.message.setText("Warning, leaving this page will delete all unsaved information for this new project.")
        self.message.setWindowTitle("Warning")
        self.message.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        ans = self.message.exec()
        if ans == QMessageBox.No:
            return 0
        else:
            self.close()
            self.go_about.emit()

    def goToHelp(self):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Warning)
        self.message.setText("Warning, leaving this page will delete all unsaved information for this new project.")
        self.message.setWindowTitle("Warning")
        self.message.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        ans = self.message.exec()
        if ans == QMessageBox.No:
            return 0
        else:
            self.close()
            self.go_help.emit()

    def goToMyProjects(self):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Warning)
        self.message.setText("Warning, leaving this page will delete all unsaved information for this new project.")
        self.message.setWindowTitle("Warning")
        self.message.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        ans = self.message.exec()
        if ans == QMessageBox.No:
            return 0
        else:
            self.close()
            self.go_my_projects.emit()

    def goToProjectPage(self):
        self.close()
        self.go_project_page.emit()

    # Info button click information
    def projectnameMethod(self, event):
        self.projectnameMessage = QMessageBox()
        self.projectnameMessage.setIcon(QMessageBox.Information)
        self.projectnameMessage.setText("Project names cannot match any previous project names.")
        self.projectnameMessage.setWindowTitle("Project Name Information")
        self.projectnameMessage.setStandardButtons(QMessageBox.Ok)

        self.projectnameMessage.exec()

    def raMethod(self, event):
        self.raMessage = QMessageBox()
        self.raMessage.setIcon(QMessageBox.Information)
        self.raMessage.setText("Right Ascension is the angular distance of a point east along the celestial equator. "
                               "It is similar to longitude on Earth. It is usually the stellar coordinate listed first.")
        self.raMessage.setWindowTitle("Right Ascension Information")
        self.raMessage.setStandardButtons(QMessageBox.Ok)

        self.raMessage.exec()

    def decMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Declination is the angular distance of a point north or south of the celestial equator. "
                             "It is similar to latitude on Earth. It is usually the stellar coordinate listed second.")
        self.message.setWindowTitle("Declination Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def astrometryMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Astrometry.net is a service to add world coordinate system (wcs) information to an image of the sky. "
                             "The wcs data allows the program to know where in the sky an image was taken, and locate stars by stellar coordinates.")
        self.message.setWindowTitle("Astrometry.net Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def calibrateMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Calibrating telescope images removes noise introduced by "
                             "heat, amplification, and gain. Only flat fields are required for calibration.")
        self.message.setWindowTitle("Calibration Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def subtractbiasMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Some telescopes include the bias frame information in "
                             "the dark frame as well. In this case, the bias must be "
                             "subtracted from the dark frame before calibration to avoid "
                             "subtracting the bias out twice. (For the Great Basin Observatory "
                             "please set this to yes if using both dark and bias frames)")
        self.message.setWindowTitle("Bias Subtraction Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def backgroundMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("In telescope images there is always some noise from background "
                             "radiation. This background radiation must be subtracted out. "
                             "Select \"From Entire Image\" to take an count of the "
                             "background radiation in the entirety of the image, or "
                             "\"From Individual Stars\" to calculate a local background "
                             "count for every individual star (recommended when there may be a gradient in the image).")
        self.message.setWindowTitle("Background Radiation Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def catalogMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Photometry+ uses Vizier catalogs to find magnitudes for reference "
                             "stars. Enter the catalog you'd like to use to search for reference star "
                             "magnitudes. To use SIMBAD enter \"SIMBAD\". If using uploaded "
                             "reference stars, feel free to leave this entry blank.")
        self.message.setWindowTitle("Catalog Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def filterMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("The filter choice is used to find correct magnitudes for "
                             "reference stars. Please select the filter the telescope image "
                             "was taken in.")
        self.message.setWindowTitle("Filter Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def otherfilterMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("If your filter is not supported by Photometry+ right now, "
                             "please enter the filter you would like to use. "
                             "WARNING: Filter formatting must be correct for catalog "
                             "being used. Please double check entered filter formatting.")
        self.message.setWindowTitle("Filter Input Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def lclineMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Select \"Yes\" to add a light curve connecting magnitudes on "
                             "final graph. Traditional light curves typically do not have lines.")
        self.message.setWindowTitle("Light Curve Line Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def magnitudeMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Select \"Multiple Calculations\" to calculate a target object magnitude for every reference star, and then average those magnitudes together for the final magnitude. Select \"One Calculation\" to average the reference star counts and magnitudes together, and then calculate one final target object magnitude. ")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def errorMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Select how errors for calculations are calculated. "
                             "\"Standard Deviation\" takes the standard deviation of "
                             "all reference star magnitude calculations before averaging "
                             "them together and uses that value as error. "
                             "\"Weighted Magnitude\" calculates the magnitude based on "
                             "the photons counted in the telescope image. "
                             "\"Jack Knife\" calculates error based on the jack knife method "
                             "by Anderson et al. in \"A Simple and Direct Measure of Photometric Uncertainties\" "
                             "Or you can opt to not calculate error (not recommended).")
        self.message.setWindowTitle("Error Calculation Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def radiusMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Select how the radius of the target object is determined. "
                             "Use \"Find radius manually\" to have the radius of the star "
                             "determined by looking at the edge of the star. "
                             "Use \"Use full-width half-maximum\" to use three times the "
                             "full-width half-maximum (FWHM) as the radius. The  FWHM represents "
                             "the width of that star area where the light coming from the star is "
                             "at half the maximum light that star puts out. "
                             "You may also put in the radius to use manually.")
        self.message.setWindowTitle("Target Radius Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def pixradiusMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Enter the radius, in pixels, to use as the target object radius.")
        self.message.setWindowTitle("Pixel Radius Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def refradiusMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Select how the radius of reference stars are determined. "
                             "Set\"Same as target object radius\" to use the target object radius "
                             "for all stars. This will reduce error introduced by differences "
                             "in the magnitude of reference stars. Set \"Use radius from file\" to use the radius measurements "
                             "found in uploaded reference star file.")
        self.message.setWindowTitle("Reference Star Radius Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def outlierMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("After reference stars are used to calculate the magnitude "
                             "of the target star, those magnitudes are assigned a z-score. "
                             "The z-score of a number is the measurement of the number of "
                             "standard deviations a number is away from the mean. Removing "
                             "reference stars more than 3 standard deviations away from "
                             "the mean can reduce noise caused by flawed reference stars.")
        self.message.setWindowTitle("Remove Reference Star Outlier Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def walkthroughMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Use walkthrough mode to step through the program one step "
                             "at a time for more control over the program and a look into "
                             "what is actually happening inside of it.")
        self.message.setWindowTitle("Walkthrough Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    # Method to submit the project
    def submit(self):
        global defaultSettings
        # Reset all style sheets
        self.fileuploadbtn.setStyleSheet("text-align:center; background-color: rgb(84, 240, 255);")
        self.folderuploadbtn.setStyleSheet("text-align:center; background-color: rgb(84, 240, 255);")
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
        self.errorstarbtn.setStyleSheet("")
        # Check if anything is not sufficient
        errorMessage = self.check()
        # If there's a problem, don't submit
        if not(errorMessage == ""):
            checkMessage = QMessageBox()
            checkMessage.setIcon(QMessageBox.Warning)
            checkMessage.setText(errorMessage)
            checkMessage.setWindowTitle("Form Incomplete")
            checkMessage.setStandardButtons(QMessageBox.Ok)
            checkMessage.exec()
            return 0
        # Check if the user truly wants to submit before submitting
        ans = self.submitMessage.exec()
        if ans == QMessageBox.Cancel:
            return 0
        elif ans == QMessageBox.Yes:
            self.setSettings()
            self.goToProcessing()

    # Set the settings based on the current settings on the page
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
        if self.catalogentry.text() == "" or self.catalogentry.text() == 0 or self.catalogentry.text() == "SIMBAD" or self.catalogentry.text() == "simbad" or self.catalogentry.text() == "Simbad":
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
        # oneMagnitudeCalculation
        if self.magnitudecombo.currentText() == "One Calculation":
            omc = 1
        else:
            omc = 0
        # errorChoice
        if self.errorcombo.currentText() == "Standard Deviation":
            ec = "STD"
        elif self.errorcombo.currentText() == "Weighted Magnitude":
            ec = "WMG"
        elif self.errorcombo.currentText() == "Jack Knife":
            ec = "JKF"
        elif self.errorcombo.currentText() == "Upload Error Stars" and not(self.errorFile == ""):
            ec = self.errorFile
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
        else:
            ccf = "DEG"
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
        # RA and Dec
        if ccf == "DEC":
            r = float(self.raentry.text())
            d = float(self.decentry.text())
        else:
            tempRa = (str(self.raentry.text())).split()
            tempDec = (str(self.decentry.text())).split()
            tempRa2 = []
            tempDec2 = []
            for j in range(len(tempRa)):
                tempRa2.append(float(tempRa[j]))
                tempDec2.append(float(tempDec[j]))
            # Conversion into decimal degrees
            r = (tempRa2[0] * 15) + ((tempRa2[1] / 60) * 15) + ((tempRa2[2] / 3600) * 15)
            d = tempDec2[0] + (tempDec2[1] / 60) + (tempDec2[2] / 3600)
        # Calibration Output
        if self.calibrationFileDirectory == "":
            cof = 0
        else:
            cof = self.calibrationFileDirectory
        photometry.changeSettings(subtractBiasFromDarkFlag=sbfdf, calibrationOutputFlag=cof,
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
                                  projectName=self.projectnameentry.text(), rightAscension=r,
                                  declination=d, main=self.mainFile,
                                  dark=self.darkFile, bias=self.biasFile, flat=self.flatFile, interface=1,
                                  walkthroughMode=wm, oneMagnitudeCalculation=omc)

    # Check if everything is valid for submission
    def check(self):
        errorMessage = ""
        errorCount = 1
        # Check Main files
        if self.mainFile == "":
            errorMessage = errorMessage + "\n" + str(errorCount) + ") Main file required for photometry. Upload the main file containing the target object on the left-hand side of the screen."
            errorCount = errorCount + 1
            self.fileuploadbtn.setStyleSheet("color: red; text-align:center; background-color: rgb(84, 240, 255);")
            self.folderuploadbtn.setStyleSheet("color: red; text-align:center; background-color: rgb(84, 240, 255);")
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
            else:
                try:
                    radius = float(self.raentry.text())
                except ValueError:
                    errorMessage = errorMessage + "\n" + str(errorCount) + ") Right ascension must be a valid number."
                    errorCount = errorCount + 1
                    self.raentry.setStyleSheet("border: 1px solid red;")
            if self.decentry.text() == "" or len(self.decentry.text()) < 1:
                errorMessage = errorMessage + "\n" + str(errorCount) + ") Declination cannot be blank"
                errorCount = errorCount + 1
                self.decentry.setStyleSheet("border: 1px solid red;")
            else:
                try:
                    radius = float(self.decentry.text())
                except ValueError:
                    errorMessage = errorMessage + "\n" + str(errorCount) + ") Declination must be a valid number."
                    errorCount = errorCount + 1
                    self.decentry.setStyleSheet("border: 1px solid red;")
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
                zscore = int(self.removeoutlierentry.text())
                if zscore < 0:
                    errorMessage = errorMessage + "\n" + str(
                        errorCount) + ") Z-score for outlier removal must be a positive integer."
                    errorCount = errorCount + 1
                    self.removeoutlierentry.setStyleSheet("border: 1px solid red;")
        except ValueError:
            errorMessage = errorMessage + "\n" + str(errorCount) + ") Z-score for outlier removal must be an integer."
            errorCount = errorCount + 1
            self.removeoutlierentry.setStyleSheet("border: 1px solid red;")
        # Check error stars
        if self.errorcombo.currentText() == "Upload Error Stars" and self.errorFile == "":
            errorMessage = errorMessage + "\n" + str(errorCount) + ") Error Stars error calculation method required an uploaded file containing error stars."
            errorCount = errorCount + 1
            self.errorstarbtn.setStyleSheet("border: 1px solid red;")
        return errorMessage

    # Get main FITS files
    def getMain(self, button):
        global displayFile
        # Check whether user is uploading one file or a folder
        if button == "FILE":
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                      "All Files (*);;FITS Files (*.fits);;FITS Files (*.fts)", options=options)
        elif button == "FOLDER":
            fileName = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        # Ensure it is a valid filename
        try:
            fileName[len(fileName) - 1]
        except IndexError:
            return 0
        # Ensure it is a FITS file
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
            # Set image to uploaded file
            data = hdul[0].data
            image = Image.fromarray(data.astype('uint8'), 'L')
            enhancer = ImageEnhance.Contrast(image)
            img = enhancer.enhance(2)
            dat = img.tobytes("raw", "L")
            qImg = QImage(dat, data.shape[1], data.shape[0], data.shape[0], QImage.Format_Grayscale8)
            self.mainPixmap = QPixmap(qImg)
            self.imagelbl.setPixmap(self.mainPixmap)
            self.imagelbl.setScaledContents(True)
            n = fileName.split("/")
            na = n[len(n) - 1]
            self.filenamelbl.setText(na)
            self.mainFile = fileName
            displayFile = fileName
        # If it's not a FITS file, check if it is a folder with FITS files
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
                # Display one of the images from the folder
                data = hdul[0].data
                image = Image.fromarray(data.astype('uint8'), 'L')
                enhancer = ImageEnhance.Contrast(image)
                img = enhancer.enhance(2)
                dat = img.tobytes("raw", "L")
                qImg = QImage(dat, data.shape[1], data.shape[0], data.shape[0], QImage.Format_Grayscale8)
                self.mainPixmap = QPixmap(qImg)
                self.imagelbl.setPixmap(self.mainPixmap)
                self.imagelbl.setScaledContents(True)
                n = displayFile.split("/")
                na = n[len(n) - 1]
                self.filenamelbl.setText(na)
                self.mainFile = fileName
        # IF there are not FITS, don't do anything

    # Get dark frame(s)
    def getDark(self):
        ans = self.fileUploadMessage.exec()
        # Check if one file or folder
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
        # Try the filename
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
            n = fileName.split("/")
            na = n[len(n) - 1]
            self.darkbtn.setText("Dark: " + na)
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
                n = fileName.split("/")
                na = n[len(n) - 1]
                self.darkbtn.setText("Dark: " + na)
                self.darkFile = fileName

    # Get bias frame(s)
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
            self.biasFile = ""
            self.biasbtn.setText("Bias Frame or Bias Frame Folder")
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
            n = fileName.split("/")
            na = n[len(n) - 1]
            self.biasbtn.setText("Bias: " + na)
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
                n = fileName.split("/")
                na = n[len(n) - 1]
                self.biasbtn.setText("Bias: " + na)
                self.biasFile = fileName

    # Get flat field(s)
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
            n = fileName.split("/")
            na = n[len(n) - 1]
            self.flatbtn.setText("Flat: " + na)
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
                n = fileName.split("/")
                na = n[len(n) - 1]
                self.flatbtn.setText("Flat: " + na)
                self.flatFile = fileName

    # Get Error Stars file
    def getError(self):
        message = QMessageBox()
        message.setIcon(QMessageBox.Question)
        message.setText("Change error star read in file?")
        message.setWindowTitle("File Question")
        message.setStandardButtons(QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
        if self.errorFile == "":
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
                self.errorFile = ""
                self.errorstarbtn.setText("Error Star Read In File (.csv)")
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
            self.errorFile = fileName
            self.errorstarbtn.setText(fileName)

    # Get Reference Star file
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

    # Get directory to output calibrated files to
    def getCalibrationDirectory(self):
        fileName = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        try:
            fileName[len(fileName) - 1]
        except IndexError:
            self.calibrationFileDirectory = ""
            self.calFileButton.setText("Calibrated File Output Location")
            return 0
        n = fileName.split("/")
        na = n[len(n) - 1]
        self.calFileButton.setText(na)
        self.calibrationFileDirectory = fileName

    # Methods to display entries when not needed
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

    def errorDisable(self):
        if self.errorcombo.currentText() == "Upload Error Stars":
            self.errorstarbtn.setDisabled(False)
        else:
            self.errorstarbtn.setDisabled(True)

class Processing(QtWidgets.QMainWindow):
    # Return Signals
    go_home = QtCore.pyqtSignal()
    go_my_projects = QtCore.pyqtSignal()
    go_settings = QtCore.pyqtSignal()
    go_about = QtCore.pyqtSignal()
    go_processing = QtCore.pyqtSignal()
    go_help = QtCore.pyqtSignal()
    go_project_page = QtCore.pyqtSignal()

    def __init__(self):
        # Initialization
        super(Processing, self).__init__() # Call the inherited classes __init__ method
        global currProject
        uic.loadUi('homepage/processing.ui', self) # Load the .ui file
        self.show()
        self.projectnameentry.setText(photometry.getprojectName())
        self.raentry.setText(str(photometry.getrightAscension()))
        self.decentry.setText(str(photometry.getdeclination()))
        self.answers = []
        self.invalidFiles = []
        currProject = photometry.getprojectName()

        # Page setup
        self.homebtn.clicked.connect(lambda: self.goToHome())
        self.myprojectsbtn.clicked.connect(lambda: self.goToMyProjects())
        self.aboutbtn.clicked.connect(lambda: self.goToAbout())
        self.settingsbtn.clicked.connect(lambda: self.goToSettings())
        self.helpbtn.clicked.connect(lambda: self.goToHelp())

        # Dialog box
        self.proceedMessage = QMessageBox()
        self.proceedMessage.setIcon(QMessageBox.Information)
        self.proceedMessage.setText("Proceed with these reference stars?")
        self.proceedMessage.setWindowTitle("Submit Reference Stars")
        self.proceedMessage.setStandardButtons(QMessageBox.Yes | QMessageBox.Cancel)

        self.exitWalkthroughMessage = QMessageBox()
        self.exitWalkthroughMessage.setIcon(QMessageBox.Information)
        self.exitWalkthroughMessage.setText("Would you like to exit walkthrough mode and use these reference stars for all other files?\n")
        self.exitWalkthroughMessage.setWindowTitle("Exit Walkthrough")
        self.exitWalkthroughMessage.setStandardButtons(QMessageBox.Yes | QMessageBox.No)

        self.finishedMessage = QMessageBox()
        self.finishedMessage.setIcon(QMessageBox.Information)
        self.finishedMessage.setText("Computation Finished. Press \"Ok\" to continue.\n")
        self.finishedMessage.setWindowTitle("Computation Finished")
        self.finishedMessage.setStandardButtons(QMessageBox.Ok)

        # Display Setup
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
        n = displayFile.split("/")
        na = n[len(n) - 1]
        self.filenamelbl.setText(na)
        #self.run.emit()
        QApplication.processEvents()
        self.runPhotometry()
        self.createResultsInfo()
        photometry.plotResults(self.answers, chartTitle=photometry.getprojectName())
        photometry.printResultsToFile(self.answers)
        self.processinglbl.setText("Computation Finished")
        self.resultsBtn = QPushButton("Go To Results Page")
        self.gridLayout.addWidget(self.resultsBtn, 8, 0, 1, 2)
        self.resultsBtn.clicked.connect(self.goToProjectPage)

    # Functions
    def goToHome(self):
        self.close()
        self.go_home.emit()

    def goToProcessing(self):
        self.close()
        self.go_processing.emit()

    def goToSettings(self):
        self.close()
        self.go_settings.emit()

    def goToAbout(self):
        self.close()
        self.go_about.emit()

    def goToHelp(self):
        self.close()
        self.go_help.emit()

    def goToMyProjects(self):
        self.close()
        self.go_my_projects.emit()

    def goToProjectPage(self):
        self.close()
        self.go_project_page.emit()

    def runPhotometry(self):
        # Check if main is a single file or a directory
        if displayFile == photometry.getmain():
            mainArray = []
            temp = photometry.Photometry(photometry.getmain(), 0, 0, 0, 0, 0, 0, 0)
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
            temp = photometry.Photometry(photometry.getdark(), 0, 0, 0, 0, 0, 0, 0)
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
            temp = photometry.Photometry(photometry.getbias(), 0, 0, 0, 0, 0, 0, 0)
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
            temp = photometry.Photometry(photometry.getflat(), 0, 0, 0, 0, 0, 0, 0)
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
                        tempCheck = QCheckBox("\n" + str(ans.referenceStars[j].id) + " Magnitude: " + str(photometry.truncate(ans.referenceStars[j].magnitude,decimals=2))
                                              + " \nCalculated Target Magnitude: " + str(photometry.truncate(ans.referenceStars[j].targetMagnitude, decimals=2))
                                              + "\n")
                        tempCheck.setChecked(True)
                        tempCheck.stateChanged.connect(lambda: self.starUnchecked(ans))
                        self.gridLayout.addWidget(tempCheck, 8+j, 0, 2, 2)
                        QApplication.processEvents()
                        self.refChecks.append(tempCheck)
                        QApplication.processEvents()
                    self.submitBtn = QPushButton("Submit Reference Stars")
                    self.gridLayout.addWidget(self.submitBtn, 9+len(ans.referenceStars), 0, 1, 2)
                    QApplication.processEvents()
                    self.loop = QtCore.QEventLoop()
                    self.submitBtn.clicked.connect(lambda: self.submitReferenceStars(ans, mainArray, darkArray, biasArray, flatArray, i))
                    self.loop.exec_()
            else:
                n = mainArray[i].fileName.split("/")
                na = n[len(n) - 1]
                self.invalidFiles.append(na)

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
            self.projectscrll.ensureWidgetVisible(self.projectnameentry)
            QApplication.processEvents()
            x = len(self.refChecks) - 1
            while x >= 0:
                self.refChecks[x].hide()
                QApplication.processEvents()
                self.submitBtn.hide()
                QApplication.processEvents()
                #del self.refChecks[x]
                x = x - 1
                QApplication.processEvents()
            self.processinglbl.setText("Processing \n(Please do not\nclose this window)")
            QApplication.processEvents()
            tempStars = []
            QApplication.processEvents()
            for a in range(len(self.refChecks)):
                if self.refChecks[a].isChecked() == True:
                    tempStars.append(ans.referenceStars[a])
            QApplication.processEvents()
            refSetting = photometry.getreadInReferenceFlag()
            response2 = self.exitWalkthroughMessage.exec()
            photometry.changeSettings(readInReferenceFlag=tempStars)
            ans = self.runLetsGo(mainArray, darkArray, biasArray, flatArray, i)
            QApplication.processEvents()
            thisFile = -1
            if not(ans == 0):
                for b in range(len(self.answers)):
                    if self.answers[b].fileName == ans.fileName:
                        self.answers[b] = ans
                        thisFile = b
            QApplication.processEvents()
            if response2 == QMessageBox.No:
                photometry.changeSettings(readInReferenceFlag=refSetting)
            elif response2 == QMessageBox.Yes:
                photometry.changeSettings(walkthroughMode=0)
                if not(thisFile == -1):
                    self.answers[thisFile].referenceStars = tempStars
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
        painter.drawText(tarX-164, tarY - 120, "Target Object")
        painter.end()
        self.imagelbl.setPixmap(self.mainPixmap)
        self.imagelbl.setScaledContents(True)
        n = file.split("/")
        na = n[len(n) - 1]
        self.filenamelbl.setText(na)
        hdul.close()
        QApplication.processEvents()

    def runLetsGo(self, mainArray, darkArray, biasArray, flatArray, i):
        if photometry.getcalibrationFlag() == 0:
            ans = self.runPhotometryCopy(photometry.getrightAscension(), photometry.getdeclination(),
                                    mainArray[i].fileName)
        elif photometry.getcalibrationFlag() == 1 and photometry.getuseDarkFlag() == 1 and photometry.getuseBiasFlag() == 1:
            d = photometry.matchCal(mainArray[i].fileName, darkArray)
            b = photometry.matchCal(mainArray[i].fileName, biasArray)
            f = photometry.matchCal(mainArray[i].fileName, flatArray)
            ans = self.runPhotometryCopy(photometry.getrightAscension(), photometry.getdeclination(),
                                    mainArray[i].fileName, d, b, f)
        elif photometry.getcalibrationFlag() == 1 and photometry.getuseDarkFlag() == 0 and photometry.getuseBiasFlag() == 1:
            b = photometry.matchCal(mainArray[i].fileName, biasArray)
            f = photometry.matchCal(mainArray[i].fileName, flatArray)
            ans = self.runPhotometryCopy(photometry.getrightAscension(), photometry.getdeclination(),
                                    mainArray[i].fileName, biasFrame=b, flatField=f)
        elif photometry.getcalibrationFlag() == 1 and photometry.getuseDarkFlag() == 1 and photometry.getuseBiasFlag() == 0:
            d = photometry.matchCal(mainArray[i].fileName, darkArray)
            f = photometry.matchCal(mainArray[i].fileName, flatArray)
            ans = self.runPhotometryCopy(photometry.getrightAscension(), photometry.getdeclination(),
                                    mainArray[i].fileName, darkFrame=d, flatField=f)
        else:
            f = photometry.matchCal(mainArray[i].fileName, flatArray)
            ans = self.runPhotometryCopy(photometry.getrightAscension(), photometry.getdeclination(),
                                    mainArray[i].fileName, flatField=f)
        return ans

    def createResultsInfo(self):
        if not os.path.exists("PhotPSaveData"):
            os.mkdir("PhotPSaveData")
        if not os.path.exists("PhotPSaveData/" + slugify(photometry.getprojectName())):
            os.mkdir("PhotPSaveData/" + slugify(photometry.getprojectName()))
        filename = "PhotPSaveData/" + slugify(photometry.getprojectName()) + "/resultssummary.csv"
        file = open(filename, "w")
        aveMag = 0
        aveErr = 0
        for i in range(len(self.answers)):
            aveMag = aveMag + self.answers[i].magnitude
            aveErr = aveErr + self.answers[i].error
        aveMag = aveMag / len(self.answers)
        aveErr = aveErr / len(self.answers)
        aveMag = photometry.truncate(aveMag, decimals=2)
        aveErr = photometry.truncate(aveErr, decimals=6)
        rChi = photometry.reducedChiSquared(self.answers)
        file.write("Info," + str(photometry.getprojectName()).replace(',', '') + "," + str(aveMag) + "," + str(aveErr) + "," + str(rChi) + ","
                   + str(photometry.getrightAscension()) + "," + str(photometry.getdeclination()) + ", \n")
        for i in range(len(self.invalidFiles)):
            file.write(str(self.invalidFiles[i]) + ", \n")
        file.close()

    def runPhotometryCopy(self, targetStarRA, targetStarDec, mainFile, darkFrame="", biasFrame="", flatField=""):
        self.processingBar.setValue(0)
        self.processinglbl.setText("Calibrating files\n(Please do not close this window)")
        QApplication.processEvents()
        if photometry.getwalkthroughMode() == 1:
            time.sleep(2)
        if photometry.getcalibrationFlag() == 1:
            # Calibrate the image
            hdul = photometry.calibrate(mainFile, darkFrame, biasFrame, flatField)
        else:
            # Use the raw image
            hdul = fits.open(mainFile)
        try:
            if hdul == 0:
                return 0
        except ValueError:
            hdul = hdul

        self.processingBar.setValue(1)
        self.processinglbl.setText("Translating from pixels to stellar coordinates\n(Please do not close this window)")
        QApplication.processEvents()
        if photometry.getwalkthroughMode() == 1:
            time.sleep(2)
        # Calculate magnitude
        # w is the reference of world coordinates for this image
        w = photometry.worldCoordinateSystem(hdul)

        # Convert COORDINATE system data into pixel locations for the image
        X, Y = w.all_world2pix(targetStarRA, targetStarDec, 0)

        # Converts pixel values to integers
        try:
            Y = int(Y)  # was pixRA
            X = int(X)  # was pixDec
        except ValueError:
            return 0

        # Check if the file has WCS data
        try:
            hdul[0].header['CRVAL1']
        except KeyError:
            if photometry.getastrometryDotNetFlag() == 0:
                self.processingBar.setValue(2)
                self.processinglbl.setText(
                    "Photometry on" + str(mainFile) + "aborted due to localization error.")
                QApplication.processEvents()
                if photometry.getwalkthroughMode() == 1:
                    time.sleep(2)
                return 0
            else:
                self.processingBar.setValue(2)
                self.processinglbl.setText(
                    "File has no WCS data. Running Astrometry.net\n(Please do not close this window)")
                QApplication.processEvents()
                if photometry.getwalkthroughMode() == 1:
                    time.sleep(2)
                # If it doesn't have WCS, add it
                mainFile = photometry.getWCS(mainFile, ra=targetStarRA, dec=targetStarDec)
                if mainFile == 0:
                    return 0
                if photometry.getcalibrationFlag() == 1:
                    # Calibrate the image
                    hdul = photometry.calibrate(mainFile, darkFrame, biasFrame, flatField)
                else:
                    # Use the raw image
                    hdul = fits.open(mainFile)

        # Calculate magnitude
        # w is the reference of world coordinates for this image
        w = photometry.worldCoordinateSystem(hdul)

        # Convert COORDINATE system data into pixel locations for the image
        X, Y = w.all_world2pix(targetStarRA, targetStarDec, 0)

        # Converts pixel values to integers
        Y = int(Y)  # was pixRA
        X = int(X)  # was pixDec

        self.processingBar.setValue(3)
        self.processinglbl.setText("Finding radius of target object\n(Please do not close this window)")
        QApplication.processEvents()
        if photometry.getwalkthroughMode() == 1:
            time.sleep(2)
        # Center the star (IN PROGRESS)
        # pixRA, pixDec = findCenter(int(Y), int(X), hdul[0].data)

        if photometry.getfwhmFlag() == 1:
            try:
                # Check if the file has FWHM
                radius = int(3 * hdul[0].header['FWHM'])
            except KeyError:
                # If not, find radius manually
                radius = photometry.findRadius(Y, X, hdul[0].data)
        elif photometry.getfwhmFlag() == 0:
            # Set the radius to the distance from the center
            # of the star to the farthest edge of the star
            radius = photometry.findRadius(Y, X, hdul[0].data)
        else:
            radius = photometry.getfwhmFlag()
        # Find the photon counts per pixel of blank sky
        photometry.changeSettings(universalBlank=photometry.findBlank(hdul[0].data, radius))

        self.processingBar.setValue(4)
        self.processinglbl.setText("Counting photons inside target object aperture\n(Please do not close this window)")
        QApplication.processEvents()
        if photometry.getwalkthroughMode() == 1:
            time.sleep(2)
        # Find the photon counts in the target star
        targetStarPhotons, targetStarError, extraStar = photometry.starCount(Y, X, hdul[0].data, radius)

        self.processingBar.setValue(5)
        self.processinglbl.setText("Finding reference stars\n(Please do not close this window)")
        QApplication.processEvents()
        if photometry.getwalkthroughMode() == 1:
            time.sleep(2)
        # Find reference stars
        readInReferenceFilename = "0"
        if not (photometry.getreadInReferenceFlag() == 0):
            if not (isinstance(photometry.getreadInReferenceFlag(), str)):
                stars = photometry.getreadInReferenceFlag()
            elif photometry.getreadInReferenceFlag()[len(photometry.getreadInReferenceFlag()) - 1] == 'v' and \
                    photometry.getreadInReferenceFlag()[len(photometry.getreadInReferenceFlag()) - 2] == 's' and \
                    photometry.getreadInReferenceFlag()[len(photometry.getreadInReferenceFlag()) - 3] == 'c':
                readInReferenceFilename = photometry.getreadInReferenceFlag()
            else:
                for filename in glob.glob(os.path.join(photometry.getreadInReferenceFlag(), '*.csv')):
                    with open(os.path.join(os.getcwd(), filename), 'r') as f:
                        readInReferenceFilename = filename
            # Read in reference stars from file
            if readInReferenceFilename == "0" and isinstance(photometry.getreadInReferenceFlag(), str):
                stars = photometry.findOtherStars(Y, X, hdul[0].data, radius, w)
                if stars == 0:
                    return 0
            elif not (readInReferenceFilename == "0") and isinstance(photometry.getreadInReferenceFlag(), str):
                stars = photometry.readFromFile(readInReferenceFilename, radius, hdul[0].data, w)
                if stars == 0:
                    return 0
                    stars = findOtherStars(Y, X, hdul[0].data, radius, w)
                    if stars == 0:
                        return 0
        else:
            # Finding new stars automatically
            stars = photometry.findOtherStars(Y, X, hdul[0].data, radius, w)
            if stars == 0:
                return 0

        # Calculate magnitudes, average, and error
        self.processingBar.setValue(6)
        self.processinglbl.setText("Calculating magnitude and error of target object\n(Please do not close this window)")
        QApplication.processEvents()
        if photometry.getwalkthroughMode() == 1:
            time.sleep(2)
        a = radius * radius * math.pi
        sigmaBkg = math.sqrt(abs(a * photometry.getuniversalBlank()))
        ave, error, stars = photometry.calculateMagnitudeAndError(targetStarPhotons, stars, targetStarError, sigmaBkg, file=mainFile, dark=darkFrame, bias=biasFrame, flat=flatField)
        if ave == 0 and error == 0:
            return 0

        if not (photometry.getremoveReferenceOutliersFlag() == 0):
            self.processingBar.setValue(7)
            self.processinglbl.setText("Removing outlier reference stars\n(Please do not close this window)")
            QApplication.processEvents()
            if photometry.getwalkthroughMode() == 1:
                time.sleep(2)
            firstAve = ave
            firstError = error
            # Remove outliers
            stars = photometry.removeReferenceOutliers(stars)
            self.processingBar.setValue(8)
            self.processinglbl.setText("Recalculating magnitude and error without outliers\n(Please do not close this window)")
            QApplication.processEvents()
            if photometry.getwalkthroughMode() == 1:
                time.sleep(2)
            # Recalculate average magnitude and standard deviation without outliers
            ave, error, stars = photometry.calculateMagnitudeAndError(targetStarPhotons, stars, targetStarError, sigmaBkg, file=mainFile, dark=darkFrame, bias=biasFrame, flat=flatField)
            if ave == 0 and error == 0:
                ave = firstAve
                error = firstError

        self.processingBar.setValue(9)
        self.processinglbl.setText("Exporting results\n(Please do not close this window)")
        QApplication.processEvents()
        if photometry.getwalkthroughMode() == 1:
            time.sleep(2)

        # Printing reference stars to files
        if photometry.getprintReferenceStarsFlag() == 1:
            n = mainFile.split("/")
            na = n[len(n) - 1]
            photometry.printReferenceToFile(stars, filename=na[0:10] + "_referencestars.csv")

        # Create and return the results of the photometry
        ans = photometry.Photometry(mainFile, hdul[0].header['JD'], ave, error, stars, targetStarPhotons, radius, photometry.getuniversalBlank())
        self.processingBar.setValue(10)
        QApplication.processEvents()
        return ans

class ProjectPage(QtWidgets.QMainWindow):
    # Return Signals
    go_home = QtCore.pyqtSignal()
    go_my_projects = QtCore.pyqtSignal()
    go_settings = QtCore.pyqtSignal()
    go_about = QtCore.pyqtSignal()
    go_processing = QtCore.pyqtSignal()
    go_help = QtCore.pyqtSignal()
    go_project_page = QtCore.pyqtSignal()

    def __init__(self):
        # Initialization
        super(ProjectPage, self).__init__()  # Call the inherited classes __init__ method
        global currProject
        uic.loadUi('homepage/projectpage.ui', self)  # Load the .ui file
        self.show()
        QApplication.processEvents()

        # Buttons
        self.homebtn.clicked.connect(lambda: self.goToHome())
        self.myprojectsbtn.clicked.connect(lambda: self.goToMyProjects())
        self.aboutbtn.clicked.connect(lambda: self.goToAbout())
        self.settingsbtn.clicked.connect(lambda: self.goToSettings())
        self.helpbtn.clicked.connect(lambda: self.goToHelp())
        self.downloadbtn.clicked.connect(self.download)
        self.datadownloadbtn.clicked.connect(self.downloadData)

        # Info
        self.magnitudeinfo.mousePressEvent = self.magnitudeMethod
        self.errorinfo.mousePressEvent = self.errorMethod
        self.chiinfo.mousePressEvent = self.chiMethod
        self.filesinfo.mousePressEvent = self.fileMethod

        # Project to show
        self.projectName = slugify(currProject)
        self.stars = self.getStarFromFile()
        self.answers = self.getPhotometryFromFile()
        self.aveMag, self.aveErr, self.rChi, self.ra, self.dec, self.invalidFiles = self.getSummaryFromFile()

        # Show info
        self.lightCurvePixmap = QPixmap("PhotPSaveData/" + slugify(self.projectName) + "/lightcurve.png")
        self.imagelbl.setPixmap(self.lightCurvePixmap)
        self.projectnameentry.setText(currProject)
        self.raentry.setText(str(self.ra))
        self.decentry.setText(str(self.dec))
        if self.aveMag == 0:
            self.avemaglbl.setText("Average Target Object Magnitude cannot be found")
        else:
            self.avemaglbl.setText("Average Target Object Magnitude: " + str(self.aveMag))
        if self.aveErr == 0:
            self.aveerrlbl.setText("Average Measurement Error cannot be found")
        else:
            self.aveerrlbl.setText("Average Measurement Error: " + str(self.aveErr))
        if self.rChi == 0 and self.aveMag == 0:
            self.rchilbl.setText("Reduced χ² cannot be found")
        else:
            self.rchilbl.setText("Reduced χ²: " + str(photometry.truncate(self.rChi, 2)))
        if self.invalidFiles == 0:
            self.invalidfileslbl.setText("Files Discarded cannot be found")
        elif len(self.invalidFiles) <= 0:
            self.invalidfileslbl.setText("Files Discarded: N/A")
        else:
            inFiles = ""
            for i in range(len(self.invalidFiles)):
                inFiles = inFiles + str(self.invalidFiles[i])
                if not(i == len(self.invalidFiles)-1):
                    inFiles = inFiles + ", \n"
            self.invalidfileslbl.setText("Files Discarded: \n" + inFiles)
            #for i in range(len(self.answers)):


    # Functions
    def goToHome(self):
        self.close()
        self.go_home.emit()

    def goToProcessing(self):
        self.close()
        self.go_processing.emit()

    def goToSettings(self):
        self.close()
        self.go_settings.emit()

    def goToAbout(self):
        self.close()
        self.go_about.emit()

    def goToHelp(self):
        self.close()
        self.go_help.emit()

    def goToMyProjects(self):
        self.close()
        self.go_my_projects.emit()

    def goToProjectPage(self):
        self.close()
        self.go_project_page.emit()

    def magnitudeMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("The magnitude of a star is a measure of its brightness, with larger numbers indicating stars that are less bright as seen from earth.")
        self.message.setWindowTitle("Magnitude Information")
        self.message.setStandardButtons(QMessageBox.Ok)
        self.message.exec()

    def errorMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("The error represents problems with the measured magnitude. A larger error can indicate problems with the telescope, problems with atmosphere in the image, or crowded fields.")
        self.message.setWindowTitle("Error Information")
        self.message.setStandardButtons(QMessageBox.Ok)
        self.message.exec()

    def chiMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("The reduced chi squared measurement represents the variability of an object.")
        self.message.setWindowTitle("Reduced Chi Squared Information")
        self.message.setStandardButtons(QMessageBox.Ok)
        self.message.exec()

    def fileMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Files can be discarded for a variety of reasons. This can include a lack of world coordinate system data for the image, interference in the telescope or atmosphere, or if there are not valid reference stars in the field. We recommend looking at these files to ascertain the problem.")
        self.message.setWindowTitle("Magnitude Information")
        self.message.setStandardButtons(QMessageBox.Ok)
        self.message.exec()

    def getStarFromFile(self):
        fileName = "PhotPSaveData/" + slugify(self.projectName) + "/referencestars.csv"
        try:
            file = open(fileName, "r")
        except FileNotFoundError:
            return 0
        stars = []
        for line in file:
            array = line.split(",")
            if not (array[0] == 'Name') and not (array[0] == '\n') and not (array[0] == '') and not (array[0] == 'ID'):
                stars.append(photometry.Star(array[0], float(array[1]), float(array[2]), float(array[3]), 0, float(array[5]),
                             float(array[6]), 0))
        file.close()
        return stars

    def getPhotometryFromFile(self):
        filename = "PhotPSaveData/" + slugify(self.projectName) + "/results.csv"
        try:
            file = open(filename, "r")
        except FileNotFoundError:
            return 0
        # Reads in data
        results = []
        for line in file:
            array = line.split(",")
            if not(array[0] == "File Name") and not (array[0] == '\n') and not (array[0] == '') and not(array[0] == 'X'):
                results.append(photometry.Photometry(array[0], array[1], array[2], array[3], 0, 0, 0, 0))
        file.close()
        return results

    def getSummaryFromFile(self):
        filename = "PhotPSaveData/" + slugify(self.projectName) + "/resultssummary.csv"
        try:
            file = open(filename, "r")
        except FileNotFoundError:
            return 0,0,0,0,0,0
        # Reads in data
        aveMag = 0
        aveErr = 0
        rChi = 0
        ra = 0
        dec = 0
        invalidFiles = []
        for line in file:
            array = line.split(",")
            if array[0] == "Info" or array[0] == '\n' or array[0] == '':
                # project name = array[1]
                aveMag = float(array[2])
                aveErr = float(array[3])
                rChi = float(array[4])
                ra = float(array[5])
                dec = float(array[6])
            else:
                invalidFiles.append(array[0])
        return aveMag, aveErr, rChi, ra, dec, invalidFiles

    def download(self):
        fileName = str(QFileDialog.getExistingDirectory(self, "Select Download Directory"))
        try:
            fileName[len(fileName) - 1]
        except IndexError:
            return 0
        #file = fileName + "/" + slugify(self.projectName) + "-LightCurve.png"
        #self.lightCurvePixmap.save(file)
        copyfile("PhotPSaveData/" + slugify(self.projectName) + "/lightcurve.png",
                 fileName + "/" + slugify(self.projectName) + "-LightCurve.png")

    def downloadData(self):
        fileName = str(QFileDialog.getExistingDirectory(self, "Select Download Directory"))
        try:
            fileName[len(fileName) - 1]
        except IndexError:
            return 0
        copyfile("PhotPSaveData/" + slugify(self.projectName) + "/results.csv", fileName + "/" + slugify(self.projectName) + "-resultsfile.csv")

class MyProjects(QtWidgets.QMainWindow):
    # Return Signals
    go_home = QtCore.pyqtSignal()
    go_my_projects = QtCore.pyqtSignal()
    go_settings = QtCore.pyqtSignal()
    go_about = QtCore.pyqtSignal()
    go_processing = QtCore.pyqtSignal()
    go_help = QtCore.pyqtSignal()
    go_project_page = QtCore.pyqtSignal()
    go_newProject = QtCore.pyqtSignal()

    def __init__(self):
        # Initialization
        super(MyProjects, self).__init__()  # Call the inherited classes __init__ method
        global currProject
        uic.loadUi('homepage/myprojects.ui', self)  # Load the .ui file
        self.show()
        QApplication.processEvents()

        # Get projects
        self.getProjects()

        # Buttons
        self.homebtn.clicked.connect(lambda: self.goToHome())
        self.aboutbtn.clicked.connect(lambda: self.goToAbout())
        self.settingsbtn.clicked.connect(lambda: self.goToSettings())
        self.helpbtn.clicked.connect(lambda: self.goToHelp())
        self.newprojectbtn.clicked.connect(lambda: self.goToNewProject())
        self.displayProjects()

    # Functions
    def goToHome(self):
        self.close()
        self.go_home.emit()

    def goToNewProject(self):
        self.close()
        self.go_newProject.emit()

    def goToProcessing(self):
        self.close()
        self.go_processing.emit()

    def goToSettings(self):
        self.close()
        self.go_settings.emit()

    def goToAbout(self):
        self.close()
        self.go_about.emit()

    def goToHelp(self):
        self.close()
        self.go_help.emit()

    def goToMyProjects(self):
        self.close()
        self.go_my_projects.emit()

    def goToProjectPage(self):
        self.close()
        self.go_project_page.emit()

    def getProjects(self):
        self.projects = []
        for root, dirs, files in os.walk("PhotPSaveData", topdown=False):
            for name in dirs:
                if not(name == "Calibrated"):
                    self.projects.append(os.path.join(root, name))

    def getSummaryFromFile(self, dirName):
        filename = dirName + "/resultssummary.csv"
        try:
            file = open(filename, "r")
        except FileNotFoundError:
            return 0,0,0,0
        # Reads in data
        projectName = ""
        aveMag = 0
        aveErr = 0
        rChi = 0
        for line in file:
            array = line.split(",")
            if array[0] == "Info" or array[0] == '\n' or array[0] == '':
                projectName = array[1]
                aveMag = float(array[2])
                aveErr = float(array[3])
                rChi = float(array[4])
        return projectName, aveMag, aveErr, rChi

    def displayProjects(self):
        # Display projects
        self.images = []
        self.labels = []
        self.viewbtns = []
        self.deletebtns = []
        self.lines = []
        loc = 0
        for i in range(len(self.projects)):
            projectName, aveMag, aveErr, rChi = self.getSummaryFromFile(self.projects[i])
            if not(projectName == 0 and aveMag == 0 and aveErr == 0 and rChi == 0):
                # Image display
                lightCurvePixmap = QPixmap(self.projects[i] + "/lightcurve.png")
                smaller_pixmap = lightCurvePixmap.scaled(100, 100, QtCore.Qt.KeepAspectRatio, QtCore.Qt.FastTransformation)
                imagelbl = QLabel(self)
                imagelbl.setPixmap(smaller_pixmap)
                self.gridLayout.addWidget(imagelbl, 3+loc, 0, 2, 1)
                self.images.append(imagelbl)
                # Project Name Label
                namelbl = QLabel(str(projectName))
                self.gridLayout.addWidget(namelbl, 3+loc, 1, 1, 1)
                self.labels.append(namelbl)
                # View Project Button
                view = QPushButton("View Project")
                view.clicked.connect(lambda: self.viewProject(projectName))
                self.gridLayout.addWidget(view, 3+loc, 2, 1, 1)
                self.viewbtns.append(view)
                loc = loc + 1
                # Information label
                infolbl = QLabel("Average Magnitude: " + str(aveMag) + "     Average Error: " + str(aveErr) + "     Reduced χ²: " + str(photometry.truncate(rChi, 2)))
                self.gridLayout.addWidget(infolbl, 3 + loc, 1, 1, 1)
                self.labels.append(infolbl)
                # Delete Project Button
                delete = QPushButton("Delete Project")
                delete.clicked.connect(lambda: self.deleteProject(projectName))
                self.gridLayout.addWidget(delete, 3 + loc, 2, 1, 1)
                self.deletebtns.append(delete)
                loc = loc + 1
                line = QFrame()
                line.setFrameShape(QFrame.HLine)
                line.setLineWidth(1)
                self.gridLayout.addWidget(line, 3 + loc, 0, 1, 3)
                self.lines.append(line)
                loc = loc + 1

    def viewProject(self, num):
        global currProject
        for i in range(len(self.viewbtns)):
            if self.sender() == self.viewbtns[i]:
                projectName, aveMag, aveErr, rChi = self.getSummaryFromFile(self.projects[i])
                currProject = projectName
        self.goToProjectPage()

    def deleteProject(self, num):
        loc = 0
        for i in range(len(self.deletebtns)):
            if self.sender() == self.deletebtns[i]:
                projectName, aveMag, aveErr, rChi = self.getSummaryFromFile(self.projects[i])
                loc = i
        self.deleteMessage = QMessageBox()
        self.deleteMessage.setIcon(QMessageBox.Warning)
        self.deleteMessage.setText("Are you sure you want to delete project " + str(projectName) +
                                   "? \nThis action is permanent and cannot be undone.")
        self.deleteMessage.setWindowTitle("File Deletion Warning")
        self.deleteMessage.setStandardButtons(QMessageBox.Yes | QMessageBox.Cancel)
        response = self.deleteMessage.exec()
        if response == QMessageBox.Yes:
            shutil.rmtree(self.projects[loc])
            for i in range(len(self.images)):
                self.images[i].hide()
                QApplication.processEvents()
            for i in range(len(self.labels)):
                self.labels[i].hide()
                QApplication.processEvents()
            for i in range(len(self.viewbtns)):
                self.viewbtns[i].hide()
                QApplication.processEvents()
            for i in range(len(self.deletebtns)):
                self.deletebtns[i].hide()
                QApplication.processEvents()
            for i in range(len(self.lines)):
                self.lines[i].hide()
                QApplication.processEvents()
            self.scrollArea.ensureWidgetVisible(self.projectlbl)
            QApplication.processEvents()
            self.getProjects()
            self.displayProjects()
        else:
            return 0

class About(QtWidgets.QMainWindow):
    # Return Signals
    go_home = QtCore.pyqtSignal()
    go_my_projects = QtCore.pyqtSignal()
    go_settings = QtCore.pyqtSignal()
    go_about = QtCore.pyqtSignal()
    go_processing = QtCore.pyqtSignal()
    go_help = QtCore.pyqtSignal()
    go_project_page = QtCore.pyqtSignal()

    def __init__(self):
        # Initialization
        super(About, self).__init__()  # Call the inherited classes __init__ method
        global currProject
        uic.loadUi('homepage/about.ui', self)  # Load the .ui file
        self.show()
        QApplication.processEvents()
        # Buttons
        self.homebtn.clicked.connect(lambda: self.goToHome())
        self.myprojectsbtn.clicked.connect(lambda: self.goToMyProjects())
        self.settingsbtn.clicked.connect(lambda: self.goToSettings())
        self.helpbtn.clicked.connect(lambda: self.goToHelp())
        self.downloadbtn.clicked.connect(self.download)
        # Website Links
        self.plotkinlbl.mousePressEvent = self.plotkinWebsite
        self.shawlbl.mousePressEvent = self.shawWebsite
        self.dascalulbl.mousePressEvent = self.dascaluWebsite
        self.tudorlbl.mousePressEvent = self.tudorWebsite

    # Functions
    def goToHome(self):
        self.close()
        self.go_home.emit()

    def goToProcessing(self):
        self.close()
        self.go_processing.emit()

    def goToSettings(self):
        self.close()
        self.go_settings.emit()

    def goToAbout(self):
        self.close()
        self.go_about.emit()

    def goToHelp(self):
        self.close()
        self.go_help.emit()

    def goToMyProjects(self):
        self.close()
        self.go_my_projects.emit()

    def goToProjectPage(self):
        self.close()
        self.go_project_page.emit()

    def download(self):
        fileName = str(QFileDialog.getExistingDirectory(self, "Select Download Directory"))
        try:
            fileName[len(fileName) - 1]
        except IndexError:
            return 0
        shutil.copy("assets/PhotometryPlusLogo.png", fileName)

    def plotkinWebsite(self, event):
        QDesktopServices.openUrl(QUrl("https://sites.google.com/site/richplotkin/"))

    def shawWebsite(self, event):
        QDesktopServices.openUrl(QUrl("https://aarranshaw.github.io/research.html"))

    def dascaluWebsite(self, event):
        QDesktopServices.openUrl(QUrl("https://www.cse.unr.edu/~dascalus/"))

    def tudorWebsite(self, event):
        QDesktopServices.openUrl(QUrl("https://www.linkedin.com/in/alexis-tudor-does-cse/"))

class Settings(QtWidgets.QMainWindow):
    # Return Signals
    go_home = QtCore.pyqtSignal()
    go_my_projects = QtCore.pyqtSignal()
    go_settings = QtCore.pyqtSignal()
    go_about = QtCore.pyqtSignal()
    go_processing = QtCore.pyqtSignal()
    go_help = QtCore.pyqtSignal()
    go_project_page = QtCore.pyqtSignal()

    def __init__(self):
        # Initialization
        super(Settings, self).__init__()  # Call the inherited classes __init__ method
        global currProject
        uic.loadUi('homepage/settings.ui', self)  # Load the .ui file
        self.show()
        QApplication.processEvents()

        # Buttons
        self.homebtn.clicked.connect(lambda: self.goToHome())
        self.myprojectsbtn.clicked.connect(lambda: self.goToMyProjects())
        self.aboutbtn.clicked.connect(lambda: self.goToAbout())
        self.helpbtn.clicked.connect(lambda: self.goToHelp())
        self.submitbtn.clicked.connect(self.setSettings)

        # Info
        self.astrometryinfo.setToolTip("Astrometry.net is a service to add world coordinate\n" +
                                       "system (wcs) information to an image of the sky. The\n" +
                                       "wcs data allows the program to know where in the sky\n" +
                                       "an image was taken, and locate stars by stellar coordinates.")

        self.astrometryinfo.mousePressEvent = self.astrometryMethod

        self.calibrateinfo.setToolTip("Calibrating telescope images removes noise introduced by\n" +
                                      "heat, amplification, and gain. Only flat fields are required\n" +
                                      "for calibration.")

        self.calibrateinfo.mousePressEvent = self.calibrateMethod

        self.subtractbiasinfo.setToolTip("Some telescopes include the bias frame information in\n" +
                                         "the dark frame as well. In this case, the bias must be\n" +
                                         "subtracted from the dark frame before calibration to avoid\n" +
                                         "subtracting the bias out twice. (For the Great Basin Observatory\n" +
                                         "please set this to yes if using both dark and bias frames)")

        self.subtractbiasinfo.mousePressEvent = self.subtractbiasMethod

        self.backgroundinfo.setToolTip("In telescope images there is always some noise from background\n" +
                                       "radiation. This background radiation must be subtracted out.\n" +
                                       "Select \"From Entire Image\" to take an count of the\n" +
                                       "background radiation in the entirety of the image, or \n" +
                                       "\"From Individual Stars\" to calculate a local background\n" +
                                       "count for every individual star (recommended when there may\n" +
                                       "be a gradient in the image).")

        self.backgroundinfo.mousePressEvent = self.backgroundMethod

        self.cataloginfo.setToolTip("Photometry+ uses Vizier catalogs to find magnitudes for reference\n" +
                                    "stars. Enter the catalog you'd like to use to search for reference star\n" +
                                    "magnitudes. To use SIMBAD enter \"SIMBAD\". If using uploaded\n" +
                                    "reference stars, feel free to leave this entry blank.")

        self.cataloginfo.mousePressEvent = self.catalogMethod

        self.filterinfo.setToolTip("The filter choice is used to find correct magnitudes for\n" +
                                   "reference stars. Please select the filter the telescope image\n"
                                   "was taken in.")

        self.filterinfo.mousePressEvent = self.filterMethod

        self.otherfilterinfo.setToolTip("If your filter is not supported by Photometry+ right now,\n"
                                        + "please enter the filter you would like to use.\n" +
                                        "WARNING: Filter formatting must be correct for catalog\n" +
                                        "being used. Please double check entered filter formatting.")

        self.otherfilterinfo.mousePressEvent = self.otherfilterMethod

        self.lightcurvelineinfo.setToolTip("Select \"Yes\" to add a light curve connecting magnitudes on\n" +
                                           "final graph. Traditional light curves typically do not have lines.")

        self.lightcurvelineinfo.mousePressEvent = self.lclineMethod

        self.magnitudeinfo.mousePressEvent = self.magnitudeMethod

        self.errorinfo.setToolTip("Select how errors for calculations are calculated.\n" +
                                  "\"Standard Deviation\" takes the standard deviation of\n" +
                                  "all reference star magnitude calculations before averaging\n" +
                                  "them together and uses that value as error.\n" +
                                  "\"Weighted Magnitude\" calculates the magnitude based on\n" +
                                  "the photons counted in the telescope image.\n" +
                                  "\"Jack Knife\" calculates error based on the jack knife method\n" +
                                  "by Anderson et al. in \"A Simple and Direct Measure of Photometric Uncertainties\"\n" +
                                  "Or you can opt to not calculate error (not recommended).")

        self.errorinfo.mousePressEvent = self.errorMethod

        self.targetradiusinfo.setToolTip("Select how the radius of the target object is determined.\n" +
                                         "Use \"Find radius manually\" to have the radius of the star\n" +
                                         "determined by looking at the edge of the star.\n" +
                                         "Use \"Use full-width half-maximum\" to use three times the \n" +
                                         "full-width half-maximum (FWHM) as the radius. The  FWHM represents\n" +
                                         "the width of that star area where the light coming from the star is\n" +
                                         "at half the maximum light that star puts out.\n" +
                                         "You may also put in the radius to use manually.")

        self.targetradiusinfo.mousePressEvent = self.radiusMethod

        self.pixradiusinfo.setToolTip("Enter the radius, in pixels, to use as the target object radius.")

        self.pixradiusinfo.mousePressEvent = self.pixradiusMethod

        self.refradiusinfo.setToolTip("Select how the radius of reference stars are determined.\n" +
                                      "Set\"Same as target object radius\" to use the target object radius\n" +
                                      "for all stars. This will reduce error introduced by differences\n" +
                                      "in the magnitude of reference stars." +
                                      "Set \"Use radius from file\" to use the radius measurements\n" +
                                      "found in uploaded reference star file.")

        self.refradiusinfo.mousePressEvent = self.refradiusMethod

        self.removerefinfo.setToolTip("After reference stars are used to calculate the magnitude\n" +
                                      "of the target object, those magnitudes are assigned a z-score.\n" +
                                      "The z-score of a number is the measurement of the number of\n" +
                                      "standard deviations a number is away from the mean. Removing\n" +
                                      "reference stars more than 3 standard deviations away from\n" +
                                      "the mean can reduce noise caused by flawed reference stars.")

        self.removerefinfo.mousePressEvent = self.outlierMethod

        self.walkthroughinfo.setToolTip("Use walkthrough mode to step through the program one step\n" +
                                        "at a time for more control over the program and a look into\n" +
                                        "what is actually happening inside of it.")

        self.walkthroughinfo.mousePressEvent = self.walkthroughMethod

        # Entry Initialization
        # astrometryDotNetFlag
        if not (defaultSettings.astrometryDotNetFlag == "0" or defaultSettings.astrometryDotNetFlag == 0):
            self.apikeyentry.setText(defaultSettings.astrometryDotNetFlag)
        # calibrationFlag
        if defaultSettings.calibrationFlag == "0" or defaultSettings.calibrationFlag == 0:
            self.calibratecombo.setCurrentText("No")
        else:
            self.calibratecombo.setCurrentText("Yes")
        # subtractBiasFromDarkFlag
        if defaultSettings.subtractBiasFromDarkFlag == 0 or defaultSettings.subtractBiasFromDarkFlag == "0":
            self.subtractbiascombo.setCurrentText("No")
        else:
            self.subtractbiascombo.setCurrentText("Yes")
        # blankPerStarFlag
        if defaultSettings.blankPerStarFlag == 0 or defaultSettings.blankPerStarFlag == "0":
            self.backgroundcombo.setCurrentText("From Entire Image")
        else:
            self.backgroundcombo.setCurrentText("From Individual Stars")
        # catalogChoice
        if defaultSettings.catalogChoice == 0 or defaultSettings.catalogChoice == "0":
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
        # oneMagnitudeCalculation
        if defaultSettings.oneMagnitudeCalculation == 0:
            self.magnitudecombo.setCurrentText("Multiple Calculations")
        else:
            self.magnitudecombo.setCurrentText("One Calculation")
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
        if defaultSettings.fwhmFlag == 0 or defaultSettings.fwhmFlag == "0":
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

    def goToSettings(self):
        self.close()
        self.go_settings.emit()

    def goToAbout(self):
        self.close()
        self.go_about.emit()

    def goToHelp(self):
        self.close()
        self.go_help.emit()

    def goToMyProjects(self):
        self.close()
        self.go_my_projects.emit()

    def goToProjectPage(self):
        self.close()
        self.go_project_page.emit()

    def astrometryMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Astrometry.net is a service to add world coordinate system (wcs) information to an image of the sky. "
                             "The wcs data allows the program to know where in the sky an image was taken, and locate stars by stellar coordinates.")
        self.message.setWindowTitle("Astrometry.net Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def calibrateMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Calibrating telescope images removes noise introduced by "
                             "heat, amplification, and gain. Only flat fields are required for calibration.")
        self.message.setWindowTitle("Calibration Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def subtractbiasMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Some telescopes include the bias frame information in "
                             "the dark frame as well. In this case, the bias must be "
                             "subtracted from the dark frame before calibration to avoid "
                             "subtracting the bias out twice. (For the Great Basin Observatory "
                             "please set this to yes if using both dark and bias frames)")
        self.message.setWindowTitle("Bias Subtraction Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def backgroundMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("In telescope images there is always some noise from background "
                             "radiation. This background radiation must be subtracted out. "
                             "Select \"From Entire Image\" to take an count of the "
                             "background radiation in the entirety of the image, or "
                             "\"From Individual Stars\" to calculate a local background "
                             "count for every individual star (recommended when there may be a gradient in the image).")
        self.message.setWindowTitle("Background Radiation Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def catalogMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Photometry+ uses Vizier catalogs to find magnitudes for reference "
                             "stars. Enter the catalog you'd like to use to search for reference star "
                             "magnitudes. To use SIMBAD enter \"SIMBAD\". If using uploaded "
                             "reference stars, feel free to leave this entry blank.")
        self.message.setWindowTitle("Catalog Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def filterMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("The filter choice is used to find correct magnitudes for "
                             "reference stars. Please select the filter the telescope image "
                             "was taken in.")
        self.message.setWindowTitle("Filter Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def otherfilterMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("If your filter is not supported by Photometry+ right now, "
                             "please enter the filter you would like to use. "
                             "WARNING: Filter formatting must be correct for catalog "
                             "being used. Please double check entered filter formatting.")
        self.message.setWindowTitle("Filter Input Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def lclineMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Select \"Yes\" to add a light curve connecting magnitudes on "
                             "final graph. Traditional light curves typically do not have lines.")
        self.message.setWindowTitle("Light Curve Line Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def magnitudeMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Select \"Multiple Calculations\" to calculate a target object magnitude for every reference star, and then average those magnitudes together for the final magnitude. Select \"One Calculation\" to average the reference star counts and magnitudes together, and then calculate one final target object magnitude. ")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def errorMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Select how errors for calculations are calculated. "
                             "\"Standard Deviation\" takes the standard deviation of "
                             "all reference star magnitude calculations before averaging "
                             "them together and uses that value as error. "
                             "\"Weighted Magnitude\" calculates the magnitude based on "
                             "the photons counted in the telescope image. "
                             "\"Jack Knife\" calculates error based on the jack knife method "
                             "by Anderson et al. in \"A Simple and Direct Measure of Photometric Uncertainties\" "
                             "Or you can opt to not calculate error (not recommended).")
        self.message.setWindowTitle("Error Calculation Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def radiusMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Select how the radius of the target object is determined. "
                             "Use \"Find radius manually\" to have the radius of the star "
                             "determined by looking at the edge of the star. "
                             "Use \"Use full-width half-maximum\" to use three times the "
                             "full-width half-maximum (FWHM) as the radius. The  FWHM represents "
                             "the width of that star area where the light coming from the star is "
                             "at half the maximum light that star puts out. "
                             "You may also put in the radius to use manually.")
        self.message.setWindowTitle("Target Radius Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def pixradiusMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Enter the radius, in pixels, to use as the target object radius.")
        self.message.setWindowTitle("Pixel Radius Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def refradiusMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Select how the radius of reference stars are determined. "
                             "Set\"Same as target object radius\" to use the target object radius "
                             "for all stars. This will reduce error introduced by differences "
                             "in the magnitude of reference stars. Set \"Use radius from file\" to use the radius measurements "
                             "found in uploaded reference star file.")
        self.message.setWindowTitle("Reference Star Radius Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def outlierMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("After reference stars are used to calculate the magnitude "
                             "of the target object, those magnitudes are assigned a z-score. "
                             "The z-score of a number is the measurement of the number of "
                             "standard deviations a number is away from the mean. Removing "
                             "reference stars more than 3 standard deviations away from "
                             "the mean can reduce noise caused by flawed reference stars.")
        self.message.setWindowTitle("Remove Reference Star Outlier Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

    def walkthroughMethod(self, event):
        self.message = QMessageBox()
        self.message.setIcon(QMessageBox.Information)
        self.message.setText("Use walkthrough mode to step through the program one step "
                             "at a time for more control over the program and a look into "
                             "what is actually happening inside of it.")
        self.message.setWindowTitle("Walkthrough Information")
        self.message.setStandardButtons(QMessageBox.Ok)

        self.message.exec()

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

    def setSettings(self):
        global defaultSettings
        # subtractBiasFromDarkFlag
        if self.subtractbiascombo.currentText() == "No":
            defaultSettings.subtractBiasFromDarkFlag = 0
        else:
            defaultSettings.subtractBiasFromDarkFlag = 1
        # calibrationFlag
        if self.calibratecombo.currentText() == "No":
            defaultSettings.calibrationFlag = 0
        else:
            defaultSettings.calibrationFlag = 1
        # blankPerStarFlag
        if self.backgroundcombo.currentText() == "From Entire Image":
            defaultSettings.blankPerStarFlag = 0
        else:
            defaultSettings.blankPerStarFlag = 1
        # catalogChoiceFlag
        if self.catalogentry.text() == "" or self.catalogentry.text() == 0 or self.catalogentry.text() == "0" or self.catalogentry.text() == "SIMBAD" or self.catalogentry.text() == "simbad" or self.catalogentry.text() == "Simbad" or self.catalogentry.text().lower() == "simbad":
            defaultSettings.catalogChoice = 0
        else:
            defaultSettings.catalogChoice = str(self.catalogentry.text())
        # filterChoice
        if self.filtercombo.currentText() == "Johnson V":
            defaultSettings.filterChoice = "V"
        elif self.filtercombo.currentText() == "Johnson B":
            defaultSettings.filterChoice = "B"
        elif self.filtercombo.currentText() == "g'":
            defaultSettings.filterChoice = "g"
        elif self.filtercombo.currentText() == "r'":
            defaultSettings.filterChoice = "r"
        elif self.filtercombo.currentText() == "i'":
            defaultSettings.filterChoice = "i"
        elif self.filtercombo.currentText() == "Other":
            defaultSettings.filterChoice = str(self.filterentry.text())
        # lightCurveLineFlag
        if self.lightcurvelinecombo.currentText() == "No":
            defaultSettings.lightCurveLineFlag = 0
        else:
            defaultSettings.lightCurveLineFlag = 1
        # oneMagnitudeCalculation
        if self.magnitudecombo.currentText() == "One Calculation":
            defaultSettings.oneMagnitudeCalculation = 1
        else:
            defaultSettings.oneMagnitudeCalculation = 0
        # errorChoice
        if self.errorcombo.currentText() == "Standard Deviation":
            defaultSettings.errorChoice = "STD"
        elif self.errorcombo.currentText() == "Weighted Magnitude":
            defaultSettings.errorChoice = "WMG"
        elif self.errorcombo.currentText() == "Jack Knife":
            defaultSettings.errorChoice = "JKF"
        else:
            defaultSettings.errorChoice = 0
        # fwhmFlag
        if self.targetradiuscombo.currentText() == "Find radius manually":
            defaultSettings.fwhmFlag = 0
        elif self.targetradiuscombo.currentText() == "Use full-width half-maximum":
            defaultSettings.fwhmFlag = 1
        else:
            defaultSettings.fwhmFlag = int(self.pixradentry.text())
            if defaultSettings.fwhmFlag <= 1:
                defaultSettings.fwhmFlag = 2
        # astrometryDotNetFlag
        if self.apikeyentry.text() == "":
            defaultSettings.astrometryDotNetFlag = 0
        else:
            defaultSettings.astrometryDotNetFlag = str(self.apikeyentry.text())
        # removeReferenceOutliersFlag
        defaultSettings.removeReferenceOutliersFlag = int(self.removeoutlierentry.text())
        # walkthroughMode
        if self.walkthroughcombo.currentText() == "Yes":
            defaultSettings.walkthroughMode = 1
        else:
            defaultSettings.walkthroughMode = 0
        self.saveSettings()
        saved = QMessageBox()
        saved.setIcon(QMessageBox.Information)
        saved.setText("Default settings saved!")
        saved.setWindowTitle("Success")
        saved.setStandardButtons(QMessageBox.Ok)
        saved.exec()

    def saveSettings(self):
        global defaultSettings
        photometry.saveSettings(set=defaultSettings, default=1)

class Help(QtWidgets.QMainWindow):
    # Return Signals
    go_home = QtCore.pyqtSignal()
    go_my_projects = QtCore.pyqtSignal()
    go_settings = QtCore.pyqtSignal()
    go_about = QtCore.pyqtSignal()
    go_processing = QtCore.pyqtSignal()
    go_help = QtCore.pyqtSignal()
    go_project_page = QtCore.pyqtSignal()

    def __init__(self):
        # Initialization
        super(Help, self).__init__()  # Call the inherited classes __init__ method
        global currProject
        uic.loadUi('homepage/help.ui', self)  # Load the .ui file
        self.show()
        QApplication.processEvents()

        # Buttons
        self.homebtn.clicked.connect(lambda: self.goToHome())
        self.myprojectsbtn.clicked.connect(lambda: self.goToMyProjects())
        self.aboutbtn.clicked.connect(lambda: self.goToAbout())
        self.settingsbtn.clicked.connect(lambda: self.goToSettings())
        QApplication.processEvents()

    # Functions
    def goToHome(self):
        self.close()
        self.go_home.emit()

    def goToProcessing(self):
        self.close()
        self.go_processing.emit()

    def goToSettings(self):
        self.close()
        self.go_settings.emit()

    def goToAbout(self):
        self.close()
        self.go_about.emit()

    def goToHelp(self):
        self.close()
        self.go_help.emit()

    def goToMyProjects(self):
        self.close()
        self.go_my_projects.emit()

    def goToProjectPage(self):
        self.close()
        self.go_project_page.emit()

class Controller:

    def __init__(self):
        global defaultSettings
        warnings.filterwarnings("ignore")
        photometry.changeSettings(interface=1)
        defaultSettings = photometry.readSettings("PhotPSaveData/defaultsettings.csv")
        if defaultSettings == 0:
            defaultSettings = photometry.Settings()
        pass

    def show_home(self):
        self.homePage = HomePage()
        self.homePage.go_newProject.connect(self.show_new_project)
        self.homePage.go_my_projects.connect(self.show_my_project)
        self.homePage.go_settings.connect(self.show_settings)
        self.homePage.go_help.connect(self.show_help)
        self.homePage.go_about.connect(self.show_about)
        self.homePage.show()

    def show_new_project(self):
        self.newProject = NewProject()
        self.newProject.go_home.connect(self.show_home)
        self.newProject.go_processing.connect(self.show_processing)
        self.newProject.go_my_projects.connect(self.show_my_project)
        self.newProject.go_settings.connect(self.show_settings)
        self.newProject.go_help.connect(self.show_help)
        self.newProject.go_about.connect(self.show_about)
        self.newProject.show()

    def show_processing(self):
        self.processing = Processing()
        self.processing.go_home.connect(self.show_home)
        self.processing.go_project_page.connect(self.show_project_page)
        self.processing.go_my_projects.connect(self.show_my_project)
        self.processing.go_settings.connect(self.show_settings)
        self.processing.go_help.connect(self.show_help)
        self.processing.go_about.connect(self.show_about)
        self.processing.show()

    def show_project_page(self):
        self.projectPage = ProjectPage()
        self.projectPage.go_home.connect(self.show_home)
        self.projectPage.go_my_projects.connect(self.show_my_project)
        self.projectPage.go_settings.connect(self.show_settings)
        self.projectPage.go_help.connect(self.show_help)
        self.projectPage.go_about.connect(self.show_about)
        self.projectPage.show()

    def show_my_project(self):
        self.myProjects = MyProjects()
        self.myProjects.go_home.connect(self.show_home)
        self.myProjects.go_newProject.connect(self.show_new_project)
        self.myProjects.go_project_page.connect(self.show_project_page)
        self.myProjects.go_settings.connect(self.show_settings)
        self.myProjects.go_help.connect(self.show_help)
        self.myProjects.go_about.connect(self.show_about)
        self.myProjects.show()

    def show_about(self):
        self.about = About()
        self.about.go_home.connect(self.show_home)
        self.about.go_my_projects.connect(self.show_my_project)
        self.about.go_settings.connect(self.show_settings)
        self.about.go_help.connect(self.show_help)
        self.about.show()

    def show_help(self):
        self.help = Help()
        self.help.go_home.connect(self.show_home)
        self.help.go_my_projects.connect(self.show_my_project)
        self.help.go_about.connect(self.show_about)
        self.help.go_settings.connect(self.show_settings)
        self.help.show()

    def show_settings(self):
        self.settings = Settings()
        self.settings.go_home.connect(self.show_home)
        self.settings.go_my_projects.connect(self.show_my_project)
        self.settings.go_about.connect(self.show_about)
        self.settings.go_help.connect(self.show_help)
        self.settings.show()