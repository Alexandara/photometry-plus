"""
Program written by Alexis Tudor at the University of Nevada, Reno
Email at alexisrenee1@gmail.com
Copyright and Licensing: GNU @ Alexis Tudor
"""
from __future__ import division
from __future__ import print_function

import glob
import math
import os

import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from photutils import Background2D
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import CircularAnnulus
from scipy.interpolate import make_interp_spline, BSpline
from scipy import stats
import statistics
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from django.template.defaultfilters import slugify

from client import Client

"""
Star objects contain everything that makes a star a star, including
the name of the star, its location in decimal degrees, the radius 
being used for the star, the counts calculated over that radius,
the magnitude of the star and the magnitude calculated for the 
target star using this star, and the error of the counts in the star.
"""
class Star:
    def __init__(self, id, ra, dec, radius, counts, magnitude, targetMagnitude, error):
        self.id = id
        self.ra = ra
        self.dec = dec
        self.radius = radius
        self.counts = counts
        self.magnitude = magnitude
        self.targetMagnitude = targetMagnitude
        self.error = error #Error of the star with background counts included

"""
Photometry objects contain the results of differential photometry.
It includes the name of the file being operated on, the julian
date of the observation, the magnitude of the target star, and the
error of the magnitude calculation. 
"""
class Photometry:
    def __init__(self, fileName, JD, magnitude, error, referenceStars, counts, radius, blanks):
        self.fileName = fileName
        self.JD = JD
        self.magnitude = magnitude
        self.error = error
        self.referenceStars = referenceStars
        self.targetStarCounts = counts
        self.targetStarRadius = radius
        self.fileBlanks = blanks

class Settings:
    def __init__(self):
        # projectName
        # Contains the name for the project
        self.projectName = ""

        # rightAscension
        # Contains the right ascension for the project
        self.rightAscension = 0

        # declination
        # Contains the declination for the project
        self.declination = 0

        # main
        # Contains main filepath for the project
        self.main = ""

        # dark
        # Contains dark filepath for the project
        self.dark = ""

        # bias
        # Contains bias filepath for the project
        self.bias = ""

        # flat
        # Contains main filepath for the project
        self.flat = ""

        # subtractBiasFromDarkFlag:
        # 0 = use dark and bias as-is
        # 1 = subtract bias from dark before calibrating
        #   NOTE: Dark files from the Great Basin Observatory have bias included and the bias
        #         will need to be subtracted out
        self.subtractBiasFromDarkFlag = 1 # INTERFACE

        # useDarkFlag:
        # 0 = do not use dark file for calibration
        # 1 = subtract dark file from main file
        self.useDarkFlag = 1 # INTERFACE

        # useBiasFlag:
        # 0 = do not use bias file for calibration
        # 1 = subtract bias file from main file
        self.useBiasFlag = 1  # INTERFACE

        # calibrationOutputFlag:
        # 0 = not outputting calibrated files
        # 1 = output calibrated files
        # directory = where to output calibrated files
        self.calibrationOutputFlag = 0 # REDUNDANT IN INTERFACE

        # calibrationFlag:
        # 0 = not performing calibration in main methods
        # 1 = performing calibration in main methods
        self.calibrationFlag = 1 # INTERFACE

        # blankPerStarFlag:
        # 0 = calculate blank sky counts for the entire image
        # 1 = calculate blank sky counts for each star
        self.blankPerStarFlag = 0 # INTERFACE

        # catalogChoice:
        # II/336/apass9 = use APASS catalog
        # Vizier catalog designation = use user-chosen catalog
        # 0 = use SIMBAD
        #   WARNING: Using SIMBAD may result in errors being introduced into the calculation
        #   as different sources are used for magnitudes in SIMBAD
        self.catalogChoice = "II/336/apass9" # INTERFACE

        # filterChoice:
        # V = Johnson V filter
        # B = Johnson B filter
        # g = g' filter
        # r = r' filter
        # i = i' filter
        # Other filter = user-chosen filter
        #   WARNING: using a non-supported filter requires that the user sets this flag
        #   to the correct formatting for accessing that filter in the catalog of their
        #   choice
        self.filterChoice = "V" # INTERFACE

        # lightCurveLineFlag:
        # 0 = generate light curve plots without lines (traditional)
        # 1 = generate light curve plots with smooth lines
        self.lightCurveLineFlag = 0 # INTERFACE

        # showLightCurveFlag:
        # 0 = do not show light curve plots on the screen
        # 1 = show light curve plots on the screen
        self.showLightCurveFlag = 0 # Not viable for interface

        # printLightCurveFlag:
        # 0 = do not print light curve to file
        # 1 = print light curve to file
        self.printLightCurveFlag = 1 # Not viable for interface

        # showCMDFlag:
        # 0 = do not show CMD plot on the screen
        # 1 = show CMD plot on the screen
        self.showCMDFlag = 0 # Not viable for interface

        # printCMDFlag:
        # 0 = do not print CMD to file
        # 1 = print CMD to file
        self.printCMDFlag = 1 # Not viable for interface

        # errorChoice:
        # STD = use standard deviation method for error management
        # JKF = use jack knife method for error management
        # WMG = use weighted magnitude for error management
        # file name = calculate error based on deviations in calculated magnitude of uploaded error stars
        #             vs. given magnitude
        # 0 = do not calculate error
        self.errorChoice = "STD" # INTERFACE

        # consolePrintFlag:
        # 0 = do not print to console
        # 1 = print updates to console
        #   WARNING: printing to console may slightly slow down the program
        self.consolePrintFlag = 1 # Not viable for interface

        # readInReferenceFlag:
        # 0 = do not read in reference stars from file
        # file name = use this reference file
        # directory name = detect reference file in directory
        self.readInReferenceFlag = 0 # INTERFACE

        # readInRadiusFlag:
        # 0 = use target star radius for all reference stars
        # 1 = use radius in reference file for reference stars
        self.readInRadiusFlag = 0 # INTERFACE

        # fwhmFlag:
        # 0 = detect radius for stars manually
        # 1 = use Full-Width Half-Maximum if available as radius for apertures
        # any other number = use this as the default radius
        self.fwhmFlag = 1 # INTERFACE

        # printReferenceStarsFlag:
        # 0 = do not print reference stars out
        # 1 = print out reference stars used
        self.printReferenceStarsFlag = 0 # Not viable for interface

        # astrometryDotNetFlag:
        # 0 = do not use astrometry.net
        # astrometry.net API key = use astrometry.net with this API key
        self.astrometryDotNetFlag = 0 # INTERFACE

        # astrometryTimeOutFlag:
        # 0 = run astrometry.net until it finds the job ID for the submission
        # any positive number = run astrometry.net until this many iterations
        self.astrometryTimeOutFlag = 0 # INTERFACE

        # universalBlank:
        # 0 = blank not calculated
        # any number = median blank counts
        #   NOTE: This number is for the program to use, it is not
        #         intended for use as a setting
        self.universalBlank = 0 # Not viable for interface

        # removeReferenceOutliersFlag:
        # 0 = do not remove reference star outliers
        # 3 = remove reference star outliers with absolute value z-score higher than 3
        # any other number = remove reference stars with absolute value z-scores above this number
        self.removeReferenceOutliersFlag = 3

        # coordinateChoiceFlag:
        # DEC = decimal degrees
        # DEG = degrees
        self.coordinateChoiceFlag = "DEC" # INTERFACE

        # interface:
        # 0 = use program as code
        # 1 = save program for interface
        self.interface = 0

        # walkthroughMode:
        # 0 = do not walk through
        # 1 = walk through
        self.walkthroughMode = 0

        # outputDetailedInformation
        # 0 = do not output detailed information
        # 1 = output detailed information
        self.outputDetailedInformation = 0

        # oneMagnitudeCalculation
        # 0 = calculate magnitude for each reference stars
        # 1 = calculate magnitude once using average counts and magnitudes for reference stars
        # WARNING: Not compatible with weighted magnitude error setting
        self.oneMagnitudeCalculation = 0

        # currentErrorStarsCalc
        # 0 = error stars based error not being calculated
        # 1 = error stars based error currently being calculated
        self.currentErrorStarsCalc = 0



settings = Settings()

"""
NAME:       changeSettings
RETURNS:    N/A
PARAMETERS: All available settings
PURPOSE:    To change program settings.  
"""
def changeSettings(subtractBiasFromDarkFlag=-1, calibrationOutputFlag=-1,
                   calibrationFlag=-1, blankPerStarFlag=-1,
                   catalogChoice=-1, filterChoice=-1,
                   lightCurveLineFlag=-1, showLightCurveFlag=-1,
                   errorChoice=-1, consolePrintFlag=-1,
                   readInReferenceFlag=-1, fwhmFlag=-1,
                   printReferenceStarsFlag=-1, astrometryDotNetFlag=-1,
                   removeReferenceOutliersFlag=-1, readInRadiusFlag=-1,
                   printLightCurveFlag=-1, useDarkFlag=-1,
                   astrometryTimeOutFlag=-1, showCMDFlag=-1,
                   printCMDFlag=-1, coordinateChoiceFlag=-1,
                   universalBlank=-1, useBiasFlag=-1,
                   projectName=-1, rightAscension=-1,
                   declination=-1, main=-1,
                   dark=-1, bias=-1, flat=-1, interface=-1,
                   walkthroughMode=-1, outputDetailedInformation=-1,
                   oneMagnitudeCalculation=-1):
    global settings
    if not(projectName == -1):
        settings.projectName = projectName
    if not (interface == -1):
        settings.interface = interface
    if not(rightAscension == -1):
        settings.rightAscension = rightAscension
    if not(declination == -1):
        settings.declination = declination
    if not(main == -1):
        settings.main = main
    if not(dark == -1):
        settings.dark = dark
    if not(bias == -1):
        settings.bias = bias
    if not(flat == -1):
        settings.flat = flat
    if not (subtractBiasFromDarkFlag == -1):
        settings.subtractBiasFromDarkFlag = subtractBiasFromDarkFlag
    if not (calibrationOutputFlag == -1):
        settings.calibrationOutputFlag = calibrationOutputFlag
    if not (calibrationFlag == -1):
        settings.calibrationFlag = calibrationFlag
    if not (blankPerStarFlag == -1):
        settings.blankPerStarFlag = blankPerStarFlag
    if not (catalogChoice == -1):
        settings.catalogChoice = catalogChoice
    if not (filterChoice == -1):
        settings.filterChoice = filterChoice
    if not (lightCurveLineFlag == -1):
        settings.lightCurveLineFlag = lightCurveLineFlag
    if not (showLightCurveFlag == -1):
        settings.showLightCurveFlag = showLightCurveFlag
    if not (errorChoice == -1):
        settings.errorChoice = errorChoice
    if not (consolePrintFlag == -1):
        settings.consolePrintFlag = consolePrintFlag
    if not (readInReferenceFlag == -1):
        settings.readInReferenceFlag = readInReferenceFlag
    if not (fwhmFlag == -1):
        settings.fwhmFlag = fwhmFlag
    if not (printReferenceStarsFlag == -1):
        settings.printReferenceStarsFlag = printReferenceStarsFlag
    if not (astrometryDotNetFlag == -1):
        settings.astrometryDotNetFlag = astrometryDotNetFlag
    if not (removeReferenceOutliersFlag == -1):
        settings.removeReferenceOutliersFlag = removeReferenceOutliersFlag
    if not (readInRadiusFlag == -1):
        settings.readInRadiusFlag = readInRadiusFlag
    if not (printLightCurveFlag == -1):
        settings.printLightCurveFlag = printLightCurveFlag
    if not (useDarkFlag == -1):
        settings.useDarkFlag = useDarkFlag
    if not (useBiasFlag == -1):
        settings.useBiasFlag = useBiasFlag
    if not (astrometryTimeOutFlag == -1):
        settings.astrometryTimeOutFlag = astrometryTimeOutFlag
    if not (showCMDFlag == -1):
        settings.showCMDFlag = showCMDFlag
    if not (printCMDFlag == -1):
        settings.printCMDFlag = printCMDFlag
    if not (coordinateChoiceFlag == -1):
        settings.coordinateChoiceFlag = coordinateChoiceFlag
    if not (universalBlank == -1):
        settings.universalBlank = universalBlank
    if not(walkthroughMode == -1):
        settings.walkthroughMode = walkthroughMode
    if not(outputDetailedInformation == -1):
        settings.outputDetailedInformation = outputDetailedInformation
    if not(oneMagnitudeCalculation == -1):
        settings.oneMagnitudeCalculation = oneMagnitudeCalculation

"""
NAME:       get[SETTING]
RETURNS:    Requested setting
PARAMETERS: N/A
PURPOSE:    To access settings from different files  
"""
def getsubtractBiasFromDarkFlag():
    return settings.subtractBiasFromDarkFlag
def getcalibrationOutputFlag():
    return settings.calibrationOutputFlag
def getcalibrationFlag():
    return settings.calibrationFlag
def getblankPerStarFlag():
    return settings.blankPerStarFlag
def getcatalogChoice():
    return settings.catalogChoice
def getfilterChoice():
    return settings.filterChoice
def getlightCurveLineFlag():
    return settings.lightCurveLineFlag
def getshowLightCurveFlag():
    return settings.showLightCurveFlag
def geterrorChoice():
    return settings.errorChoice
def getconsolePrintFlag():
    return settings.consolePrintFlag
def getreadInReferenceFlag():
    return settings.readInReferenceFlag
def getfwhmFlag():
    return settings.fwhmFlag
def getprintReferenceStarsFlag():
    return settings.printReferenceStarsFlag
def getastrometryDotNetFlag():
    return settings.astrometryDotNetFlag
def getremoveReferenceOutliersFlag():
    return settings.removeReferenceOutliersFlag
def getreadInRadiusFlag():
    return settings.readInRadiusFlag
def getprintLightCurveFlag():
    return settings.printLightCurveFlag
def getuseDarkFlag():
    return settings.useDarkFlag
def getuseBiasFlag():
    return settings.useBiasFlag
def getastrometryTimeOutFlag():
    return settings.astrometryTimeOutFlag
def getshowCMDFlag():
    return settings.showCMDFlag
def getprintCMDFlag():
    return settings.printCMDFlag
def getcoordinateChoiceFlag():
    return settings.coordinateChoiceFlag
def getprojectName():
    return settings.projectName
def getrightAscension():
    return settings.rightAscension
def getdeclination():
    return settings.declination
def getmain():
    return settings.main
def getdark():
    return settings.dark
def getbias():
    return settings.bias
def getflat():
    return settings.flat
def getinterface():
    return settings.interface
def getwalkthroughMode():
    return settings.walkthroughMode
def getuniversalBlank():
    return settings.universalBlank
def getoutputDetailedInformation():
    return settings.outputDetailedInformation
def getoneMagnitudeCalculation():
    return settings.oneMagnitudeCalculation

"""
NAME:       saveSettings
RETURNS:    N/A
PARAMETERS: The settings object to save (set)
            Whether this is the UI default settings or not (default)
                NOTE: default = 1 means the default program settings are being saved
PURPOSE:    To save project settings for future use
"""
def saveSettings(set=settings, default=0):
    if set.projectName == "" and default == 0:
        return 0
    # Create new output directory if none exists
    if default == 0:
        # Saving chosen settings for a specific project
        if not os.path.exists("PhotPSaveData"):
            os.mkdir("PhotPSaveData")
        if not os.path.exists("PhotPSaveData/" + slugify(set.projectName)):
            os.mkdir("PhotPSaveData/" + slugify(set.projectName))
        # Create new file
        filename = "PhotPSaveData/" + slugify(set.projectName) + "/settings.csv"
        try:
            file = open(filename, "w")
        except FileNotFoundError:
           return -1
    else:
        # Saving default settings
        if not os.path.exists("PhotPSaveData"):
            os.mkdir("PhotPSaveData")
        # Create new file
        filename = "PhotPSaveData/defaultsettings.csv"
        try:
            file = open(filename, "w")
        except FileNotFoundError:
           return -1
    # Write
    file.write("Project Name," +  str(set.projectName) + ",\n")
    file.write("Right Ascension," +  str(set.rightAscension) + ",\n")
    file.write("Declination," + str(set.declination) + ",\n")
    file.write("Main File Location," + str(set.main) + ",\n")
    file.write("Dark File Location," + str(set.dark) + ",\n")
    file.write("Bias File Location," + str(set.bias) + ",\n")
    file.write("Flat File Location," + str(set.flat) + ",\n")
    file.write("subtractBiasFromDarkFlag," + str(set.subtractBiasFromDarkFlag) + ",\n")
    file.write("calibrationOutputFlag," + str(set.calibrationOutputFlag) + ",\n")
    file.write("calibrationFlag," + str(set.calibrationFlag) + ",\n")
    file.write("blankPerStarFlag," + str(set.blankPerStarFlag) + ",\n")
    file.write("catalogChoice," + str(set.catalogChoice) + ",\n")
    file.write("filterChoice," + str(set.filterChoice) + ",\n")
    file.write("lightCurveLineFlag," + str(set.lightCurveLineFlag) + ",\n")
    file.write("showLightCurveFlag," + str(set.showLightCurveFlag) + ",\n")
    file.write("errorChoice," + str(set.errorChoice) + ",\n")
    file.write("consolePrintFlag," + str(set.consolePrintFlag) + ",\n")
    file.write("readInReferenceFlag," + str(set.readInReferenceFlag) + ",\n")
    file.write("fwhmFlag," + str(set.fwhmFlag) + ",\n")
    file.write("printReferenceStarsFlag," + str(set.printReferenceStarsFlag) + ",\n")
    file.write("astrometryDotNetFlag," + str(set.astrometryDotNetFlag) + ",\n")
    file.write("removeReferenceOutliersFlag," + str(set.removeReferenceOutliersFlag) + ",\n")
    file.write("readInRadiusFlag," + str(set.readInRadiusFlag) + ",\n")
    file.write("printLightCurveFlag," + str(set.printLightCurveFlag) + ",\n")
    file.write("useDarkFlag," + str(set.useDarkFlag) + ",\n")
    file.write("astrometryTimeOutFlag," + str(set.astrometryTimeOutFlag) + ",\n")
    file.write("showCMDFlag," + str(set.showCMDFlag) + ",\n")
    file.write("printCMDFlag," + str(set.printCMDFlag) + ",\n")
    file.write("coordinateChoiceFlag," + str(set.coordinateChoiceFlag) + ",\n")
    file.write("universalBlank," + str(set.universalBlank) + ",\n")
    file.write("useBiasFlag," + str(set.useBiasFlag) + ",\n")
    file.write("interface," + str(set.interface) + ",\n")
    file.write("walkthroughMode," + str(set.walkthroughMode) + ",\n")
    file.write("outputDetailedInformation," + str(set.outputDetailedInformation) + ",\n")
    file.write("oneMagnitudeCalculation," + str(set.oneMagnitudeCalculation) + ",\n")
    file.flush()
    file.close()

"""
NAME:       readSettings
RETURNS:    A settings object
PARAMETERS: Settings file to read in from (fileName)
PURPOSE:    To retrieve project settings from file
"""
def readSettings(fileName):
    # Make sure the file is valid
    try:
        file = open(fileName, "r")
    except FileNotFoundError:
        if settings.consolePrintFlag == 1 and settings.interface == 0:
            print("Error: Invalid file for reading in settings")
        return 0
    newSettings = Settings()
    # On every line check what setting is being read in
    for line in file:
        array = line.split(",")
        if array[0] == "Project Name":
            newSettings.projectName = array[1]
        elif array[0] == "Right Ascension":
            newSettings.rightAscension = float(array[1])
        elif array[0] == "Declination":
            newSettings.declination = float(array[1])
        elif array[0] == "Main File Location":
            newSettings.main = array[1]
        elif array[0] == "Dark File Location":
            newSettings.dark = array[1]
        elif array[0] == "Bias File Location":
            newSettings.bias = array[1]
        elif array[0] == "Flat File Location":
            newSettings.flat = array[1]
        elif array[0] == "subtractBiasFromDarkFlag":
            newSettings.subtractBiasFromDarkFlag = int(array[1])
        elif array[0] == "calibrationOutputFlag":
            newSettings.calibrationOutputFlag = int(array[1])
        elif array[0] == "calibrationFlag":
            newSettings.calibrationFlag = int(array[1])
        elif array[0] == "blankPerStarFlag":
            newSettings.blankPerStarFlag = int(array[1])
        elif array[0] == "catalogChoice":
            if array[1] == "0":
                newSettings.catalogChoice = 0
            else:
                newSettings.catalogChoice = array[1]
        elif array[0] == "filterChoice":
            newSettings.filterChoice = str(array[1])
        elif array[0] == "lightCurveLineFlag":
            newSettings.lightCurveLineFlag = int(array[1])
        elif array[0] == "showLightCurveFlag":
            newSettings.showLightCurveFlag = int(array[1])
        elif array[0] == "errorChoice":
            if array[1] == "0":
                newSettings.errorChoice = 0
            else:
                newSettings.errorChoice = array[1]
        elif array[0] == "consolePrintFlag":
            newSettings.consolePrintFlag = int(array[1])
        elif array[0] == "readInReferenceFlag":
            if array[1] == "0":
                newSettings.readInReferenceFlag = 0
            else:
                newSettings.readInReferenceFlag = array[1]
        elif array[0] == "fwhmFlag":
            newSettings.fwhmFlag = int(array[1])
        elif array[0] == "printReferenceStarsFlag":
            newSettings.printReferenceStarsFlag = int(array[1])
        elif array[0] == "astrometryDotNetFlag":
            if array[1] == "0":
                newSettings.astrometryDotNetFlag = 0
            else:
                newSettings.astrometryDotNetFlag = array[1]
        elif array[0] == "removeReferenceOutliersFlag":
            newSettings.removeReferenceOutliersFlag = int(array[1])
        elif array[0] == "readInRadiusFlag":
            newSettings.readInRadiusFlag = int(array[1])
        elif array[0] == "printLightCurveFlag":
            newSettings.printLightCurveFlag = int(array[1])
        elif array[0] == "useDarkFlag":
            newSettings.useDarkFlag = int(array[1])
        elif array[0] == "astrometryTimeOutFlag":
            newSettings.astrometryTimeOutFlag = int(array[1])
        elif array[0] == "showCMDFlag":
            newSettings.showCMDFlag = int(array[1])
        elif array[0] == "printCMDFlag":
            newSettings.printCMDFlag = int(array[1])
        elif array[0] == "coordinateChoiceFlag":
            newSettings.coordinateChoiceFlag = array[1]
        elif array[0] == "universalBlank":
            newSettings.universalBlank = float(array[1])
        elif array[0] == "useBiasFlag":
            newSettings.useBiasFlag = int(array[1])
        elif array[0] == "interface":
            newSettings.interface = int(array[1])
        elif array[0] == "walkthroughMode":
            newSettings.walkthroughMode = int(array[1])
        elif array[0] == "outputDetailedInformation":
            newSettings.outputDetailedInformation = str(array[1])
            if newSettings.outputDetailedInformation == "0" or newSettings.outputDetailedInformation == "1":
                newSettings.outputDetailedInformation = int(newSettings.outputDetailedInformation)
        elif array[0] == "oneMagnitudeCalculation":
            newSettings.oneMagnitudeCalculation = int(array[1])
    file.close()
    return newSettings

"""
NAME:       calibrate
RETURNS:    calibrated data
PARAMETERS: the file to be calibrated (filename)
            the dark frame (dark)
            the bias frame (bias)
            the flat frame (flat)
PURPOSE:    To calibrate raw .fits files into a form that can 
            be used to calculate accurate magnitude data. 
"""
def calibrate(filename, dark="", bias="", flat=""):
    global settings
    if flat == "":
        if settings.consolePrintFlag == 1:
            print("Error in calibrate: No flat field")
            return 0
    # Open the files
    try:
        hdul = fits.open(filename)  # hdul is the computer version of
                                    # the file data
    except OSError:
        if settings.consolePrintFlag == 1:
            print("Error in calibrate: Invalid filename:", filename)
        return 0
    # The following hdul files are computer readable versions of the dark, and flat .fit files
    # Each file is checked for an OSError when opening
    if settings.useDarkFlag == 1:
        try:
            hdulDark = fits.open(dark)
            dataDark = hdulDark[0].data
        except OSError:
            if settings.consolePrintFlag == 1:
                print("Error in calibrate: Invalid filename:", dark)
            return 0
    try:
        hdulFlat = fits.open(flat)
        dataFlat = hdulFlat[0].data
    except OSError:
        if settings.consolePrintFlag == 1:
            print("Error in calibrate: Invalid filename:", flat)
        return 0
    if settings.useBiasFlag == 1:
        try:
            hdulBias = fits.open(bias)
            dataBias = hdulBias[0].data
        except OSError:
            if settings.consolePrintFlag == 1:
                print("Error in calibrate: Invalid filename:", bias)
            return 0

    # Extract data from the files
    data = hdul[0].data
    # Check if bias needs to be subtracted from dark
    if settings.subtractBiasFromDarkFlag == 1 and settings.useDarkFlag == 1 and settings.useBiasFlag == 1:
        dataDark = dataDark - dataBias

    # Calibrate the files, and then reassign the calibrated data to the
    # original data
    if settings.useDarkFlag == 1:
        expTime = hdul[0].header['EXPTIME']/hdulDark[0].header['EXPTIME']
        dataDark = dataDark * expTime
        data = data - dataDark #Subtract the dark frame from the data
    if settings.useBiasFlag == 1:
        data = data - dataBias #Subtract the bias frame from the data

    dataFlatNorm = dataFlat / np.median(dataFlat) #Normalize the flat frame information
    data = data // dataFlatNorm #Divide out the normalized flat frame data
    hdul[0].data = data #This sets the calibrated data to be a part of the
                        #hdul object

    # If outputting calibrated file, create folder and file
    if settings.calibrationOutputFlag == 1:
        n = filename.split("/")
        na = n[len(n) - 1]
        if settings.interface == 0:
            if not os.path.exists("Output"):
                os.mkdir("Output")
            calibrationFile = "Output/CALIBRATED_" + na
            hdul.writeto(calibrationFile, overwrite=True)
    elif not(settings.calibrationOutputFlag == 0):
        try:
            n = filename.split("/")
            na = n[len(n) - 1]
            calibrationFile = str(settings.calibrationOutputFlag)+"/CALIBRATED_"+na
            print(calibrationFile)
            hdul.writeto(calibrationFile, overwrite=True)
        except OSError:
            if settings.consolePrintFlag == 1:
                print("Writing calibration files failed due to inaccessible directory.")
    return hdul

"""
NAME:       findBlankNoRad
RETURNS:    Average counts in blank sky
PARAMETERS: Array of image data for the sky (data)
PURPOSE:    In the case the FWHM cannot be used as the radius
            of the star, this function can be used to find the 
            counts in a blank portion of the sky with an assumed
            star radius of 25
"""
def findBlankNoRad(data):
    # Using Background2d to get a median background count
    background = Background2D(data, (100, 100))
    return background.background_median

"""
NAME:       findRadius
RETURNS:    Detected max radius
PARAMETERS: Center of star pixel Y location (Y)
            Center of star pixel X location (X)
            Image array (data)
PURPOSE:    This finds the radius of a target star by detecting
            the transition to blank sky. This is only for use in 
            cases where using the FWHM is not possible.
"""
def findRadius(Y, X, data):
    # Find the length of the image in pixels
    x, y = data.shape
    # Start at star center
    currY = Y
    currX = X
    if Y > y or X > x:
        return 20

    #Blank is the number of counts per pixel in empty sky
    blank = findBlankNoRad(data) * 1.5
    #Find the distance it takes each direction to get to blank sky
    try:
        while data[currX][currY] > blank:
            currY = currY + 1
            if currY >= y:
                currY = Y
            break
        r1 = currY - Y
        currY = Y
    except IndexError:
        r1 = 0
    try:
        while data[currX][currY] > blank:
            currY = currY - 1
            if currY <= 1:
                currY = Y
                break
        r2 = -1 * (currY - Y)
        currY = Y
    except IndexError:
        r2 = 0
    try:
        while data[currX][currY] > blank:
            currX = currX - 1
            if currX <= 1:
                currX = X * -1
                break
        r3 = -1 * (currX - X)
        currX = X
    except IndexError:
        r3 = 0
    try:
        while data[currX][currY] > blank:
            currX = currX + 1
            if currX >= x:
                currX = X
                break
        r4 = currX - X
    except IndexError:
        r4 = 0

    # Take the maximum out of the four found radius
    max = 100
    if r1 >= r2 and r1 >= r3 and r1 >= r4 and r1 < 50:
        max = r1
    elif r2 >= r1 and r2 >= r3 and r2 >= r4 and r2 < 50:
        max = r2
    elif r3 >= r1 and r3 >= r2 and r3 >= r4 and r3 < 50:
        max = r3
    elif r4 >= r1 and r4 >= r2 and r4 >= r3 and r4 < 50:
        max = r4
    if max == 100 or max == 0:
        max = 30
    return max

"""
NAME:       findBlank
RETURNS:    Average photon count of blank sky
PARAMETERS: Image array (data)
            Radius of primary target star (r)
PURPOSE:    Using the radius of the primary target star, this 
            method averages out the blank counts over an area 4 times 
            the radius of the primary target star and returns it.
"""
def findBlank(data, r):
    # Uses Background2D to find the median blanks over an area the side
    # of the star times four
    r = int(r) # make sure the radius is an integer
    try:
        background = Background2D(data, ((r * 4), (r * 4)))
    except ValueError:
        background = findBlankNoRad(data)
        return background
    return background.background_median

"""
NAME:       starCount
RETURNS:    The amount of photons in a star's image
PARAMETERS: Center of star pixel Y location (Y)
            Center of star pixel X location (X)
            Image data (data)
            Radius of the star (r)
PURPOSE:    Taking in information about any star, this method counts
            the photons in a circular area around the center of the star.
"""
def starCount(Y, X, data, r):
    global settings
    extraStar = 0
    blankMode = 0
    # Creates an aperture centered around the target star of radius r
    if r == 0: # make sure the radius is not 0
        r = 20
    if r < 0: # make sure the radius is a positive number
        r = r * -1
    targetAperture = CircularAperture((X, Y), r=r)
    targetStarTable = aperture_photometry(data, targetAperture)

    # Counts the sum in that aperture
    targetStarPhotons = targetStarTable['aperture_sum'][0]

    # Calculate Error in this count
    error = math.sqrt(abs(targetStarPhotons))

    # Check if there are stars near this star:
    blankAperture = CircularAnnulus((X, Y), r_in=r, r_out=r * 4)
    try:
        # Get array containing each pixel of blank
        blankMask = blankAperture.to_mask(method='center')
        blankData = blankMask.multiply(data)
        mask = blankMask.data
        blankArray = blankData[mask > 0]
        # Get the median
        blankMedian = np.median(blankArray)
        # Get the mode
        try:
            blankMode = statistics.mode(blankArray)
        except statistics.StatisticsError:
            blankMode = blankMedian
        # Get the mean
        blankTable = aperture_photometry(data, blankAperture)
        blankPhotons = blankTable['aperture_sum'][0]
        annulusArea = (math.pi * (r * 4) * (r * 4)) - (math.pi * r * r)
        blankMean = blankPhotons / annulusArea
        if abs(blankMean-blankMedian) > 50 or abs(blankMean-blankMode) > 50 or abs(blankMode-blankMedian) > 50:
            extraStar = 1
    except:
        blank = settings.universalBlank
    # If set to calculate per star
    if (settings.blankPerStarFlag == 1 or settings.universalBlank == 0) and not(blankMode == 0):
        blank = blankMode
    elif settings.blankPerStarFlag == 1 and blankMode == 0 and settings.universalBlank == 0:
        blank = findBlank(data, r)
    else:
        blank = settings.universalBlank

    # Subtracts blank counts per every pixel in the star
    targetStarPhotons = targetStarPhotons - ((math.pi * r * r) * blank)
    return targetStarPhotons, error, extraStar

"""
NAME:       formatFilter
RETURNS:    The filter format for Vizier, Simbad add Votable, and Simbad
PARAMETERS: None
PURPOSE:    This method locates stars nearby the target stars
            and then calls a helper method to construct a star
            object for the located stars.
"""
def formatFilter():
    global settings
    if settings.filterChoice == "V" or settings.filterChoice == "v":
        return 'Vmag', 'flux(V)', 'FLUX_V'
    elif settings.filterChoice == "B":
        return 'Bmag', 'flux(B)', 'FLUX_B'
    elif settings.filterChoice == "g":
        return 'g_mag', 'flux(g)', 'FLUX_g'
    elif settings.filterChoice == "r":
        return 'r_mag', 'flux(r)', 'FLUX_r'
    elif settings.filterChoice == "i":
        return 'i_mag', 'flux(i)', 'FLUX_i'
    else:
        return settings.filterChoice, "flux(V)", 0

"""
NAME:       findOtherStars
RETURNS:    Array of other star objects
PARAMETERS: Center of target star RA (Y)
            Center of target star DEC (X)
            Image data (data)
            Radius to examine (rad)
            Parameters for conversion from world coordinates to
            pixel location (w)
PURPOSE:    This method locates stars nearby the target stars
            and then calls a helper method to construct a star
            object for the located stars.
"""
def findOtherStars(Y, X, data, rad, w):
    global settings
    fCat, fVot, fSim = formatFilter()
    sbad = 0
    stars = []
    if settings.consolePrintFlag == 1:
        print("Starting findOtherStars")
    # Initialize SIMBAD
    customSimbad = Simbad()
    customSimbad.add_votable_fields(fVot)

    # Get sky coordinate object
    skyC = SkyCoord.from_pixel(X, Y, w)

    if not(settings.catalogChoice == 0):
        if settings.consolePrintFlag == 1:
            print("Querying Vizier Catalog: ", settings.catalogChoice)
        # Query all stars in the image
        try:
            result = Vizier.query_region(skyC, radius='0d10m0s', catalog=settings.catalogChoice)
        except requests.exceptions.ConnectionError:
            if settings.consolePrintFlag == 1:
                print("Connection Error: Check your internet connection before searching for reference stars")
            return 0
        # Check if there are results
        if not (result is None):
            if settings.consolePrintFlag == 1:
                print("Reformatting Vizier Data")
            ra = result[0]['RAJ2000']
            dec = result[0]['DEJ2000']
            mag = result[0][fCat]
            # Create star objects from the found reference stars
            for i in range(1, len(mag)):
                # This line checks for a non-null magnitude value
                if type(mag[i]) is np.float32:
                    starX, starY = w.all_world2pix(ra[i], dec[i], 0)
                    if not(starX < 0 or starX >= len(data) or starY < 0 or starY >= len(data[0])):
                        # Calculate counts in the star
                        c, sigmaSrc, extraStar = starCount(starY, starX, data, rad)
                        # Addition of this star to the array of stars objects
                        if extraStar == 0:
                            stars.append(Star("Ref"+str(i), ra[i], dec[i], rad, c, mag[i], 0, sigmaSrc))
        else:
            if settings.consolePrintFlag == 1:
                print("Catalog has no valid entries. Switching to SIMBAD.")
                print("WARNING: Using SIMBAD may introduce error into the calculations due to discrepancies in where magnitudes were obtained.")
            sbad = 1
    else:
        if settings.consolePrintFlag == 1:
            print("WARNING: Using SIMBAD may introduce error into the calculations due to discrepancies in where magnitudes were obtained.")
        sbad = 1
    if sbad == 1 and not(fSim == 0):
        if settings.consolePrintFlag == 1:
            print("Querying SIMBAD")
        # Query all stars in the image
        # Query in an area about the size of the image from GBO
        try:
            result = customSimbad.query_region(skyC, radius='0d10m0s')  # This line creates warnings
        except requests.exceptions.ConnectionError:
            if settings.consolePrintFlag == 1:
                print("Connection Error: Check your internet connection before searching for reference stars")
            return 0
        # If there are any results:
        if not (result is None):
            ra = result['RA']
            dec = result['DEC']
            mag = result[fSim]
            for i in range(1, len(mag)):
                # This line checks for a non-null magnitude value
                if type(mag[i]) is np.float32:
                    tempRa = (str(ra[i])).split()
                    tempDec = (str(dec[i])).split()
                    tempRa2 = []
                    tempDec2 = []
                    for j in range(len(tempRa)):
                        tempRa2.append(float(tempRa[j]))
                        tempDec2.append(float(tempDec[j]))
                    # Conversion into decimal degrees
                    degreeRa = (tempRa2[0]*15) + ((tempRa2[1]/60)*15) + ((tempRa2[2]/3600)*15)
                    degreeDec = tempDec2[0] + (tempDec2[1]/60) + (tempDec2[2]/3600)
                    starX, starY = w.all_world2pix(degreeRa, degreeDec, 0)
                    if not (starX < 0 or starX >= len(data) or starY < 0 or starY >= len(data[0])):
                        # Calculate counts in the star
                        c, sigmaSrc, extraStar = starCount(starY, starX, data, rad)
                        # Addition of this star to the array of stars objects
                        if extraStar == 0:
                            stars.append(Star("Ref"+str(i), degreeRa, degreeDec, rad, c, mag[i], 0, sigmaSrc))
    elif sbad == 1 and len(stars) == 0 and not(fSim == 0):
        if settings.consolePrintFlag == 1:
            print("No reference stars found.")
        return 0
    return stars

# Truncation function from https://realpython.com/python-rounding/#truncation
# Documentation by Alexis Tudor
"""
NAME:       truncate
RETURNS:    A number truncated to a certain place
PARAMETERS: A number to truncate (n)
            The decimal place to round to (decimals)
                NOTE: Default is to round to integer value
PURPOSE:    Truncates a number to a certain place value.
"""
def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier

"""
NAME:       printReferenceToFile
RETURNS:    Nothing
PARAMETERS: Array of star objects (stars)
            .csv Filename to output (filename)
            Flag for if the method is printing error stars (1) or reference stars (0) (errorStars)
PURPOSE:    This method prints out a file  in a new directory named "Output"
            containing the stars used for calculating magnitude
"""
def printReferenceToFile(stars, filename="stars.csv", errorStars=0):
    global settings
    if settings.interface == 0:
        # Create new output directory if none exists
        if not os.path.exists("Output"):
            os.mkdir("Output")
        # Create new file
        filename = "Output/" + filename
    elif settings.interface == 1:
        if not os.path.exists("PhotPSaveData"):
            os.mkdir("PhotPSaveData")
        if not os.path.exists("PhotPSaveData/" + slugify(settings.projectName)):
            os.mkdir("PhotPSaveData/" + slugify(settings.projectName))
        filename = "PhotPSaveData/" + slugify(settings.projectName) + "/referencestars.csv"
    try:
        file = open(filename, "w")
    except FileNotFoundError:
        n = filename.split("/")
        na = n[len(n) - 1]
        filename = "Output/" + na
        try:
            file = open(filename, "w")
        except FileNotFoundError:
            file = open("Output/INVALIDFILENAMESTARS.csv", "w")
    # Create the heading
    if errorStars == 1:
        # Error Stars
        file.write(
            "ID,Right Ascension (Decimal Degrees),Declination (Decimal Degrees),Radius (pixels),Net Counts,Magnitude,Delta Magnitude,Final Error,\n")
    else:
        # Reference Stars
        file.write("ID,Right Ascension (Decimal Degrees),Declination (Decimal Degrees),Radius (pixels),Net Counts,Magnitude,Calculated Target Magnitude,Error,\n")
    # Make sure stars is an array
    try:
        len(stars)
    except TypeError:
        stars = [stars]
    # For each star, print all categories
    for i in range(len(stars)):
        trRA = truncate(stars[i].ra, 6)
        trDec = truncate(stars[i].dec, 6)
        trC = truncate(stars[i].counts)
        trTM = truncate(stars[i].targetMagnitude, 2)
        trE = truncate(stars[i].error, 6)
        file.write(str(stars[i].id) + "," + str(trRA) + "," + str(trDec)
        + "," + str(stars[i].radius) + "," + str(trC) +
                   "," + str(stars[i].magnitude) + "," + str(trTM) + "," + str(trE) + ", \n")
    file.flush()
    file.close()

"""
NAME:       printDetailedInfo
RETURNS:    Nothing
PARAMETERS: Array of Photometry objects (results)
            Name of the folder to put detailed info in (dirName)
PURPOSE:    This method prints out a file  in a new directory named "Output"
            containing detailed information on the calculations
"""
def printDetailedInfo(results, dirName="DetailedDirectory"):
    if settings.interface == 0 and settings.outputDetailedInformation == 1:
        # Create new output directory if none exists
        if not os.path.exists("Output"):
            os.mkdir("Output")
        # Create new file
        n = dirName.split("/")
        na = n[len(n) - 1]
        if len(na) > 6:
            filename = na[0:6] + "_detailedinformation.csv"
        else:
            filename = na + "_detailedinformation.csv"
        filename = "Output/" + filename
    elif settings.interface == 0 and not(settings.outputDetailedInformation == 1):
        # Create new output directory if none exists
        if not os.path.exists("Output"):
            os.mkdir("Output")
        # Create new file
        settings.outputDetailedInformation = str(settings.outputDetailedInformation)
        if len(settings.outputDetailedInformation) > 4:
            if settings.outputDetailedInformation[len(settings.outputDetailedInformation) - 1] == 'v' and settings.outputDetailedInformation[len(settings.outputDetailedInformation) - 2] == 's' and settings.outputDetailedInformation[len(settings.outputDetailedInformation) - 3] == 'c' and settings.outputDetailedInformation[len(settings.outputDetailedInformation) - 4] == '.':
                filename = "Output/" + str(settings.outputDetailedInformation)
            else:
                filename = "Output/" + str(settings.outputDetailedInformation) + ".csv"
        else:
            "Output/INVALIDFILENAMEDETAILEDRESULTS.csv"
    # In the off chance this is being used by the interface
    elif settings.interface == 1:
        if not os.path.exists("PhotPSaveData"):
            os.mkdir("PhotPSaveData")
        if not os.path.exists("PhotPSaveData/" + slugify(settings.projectName)):
            os.mkdir("PhotPSaveData/" + slugify(settings.projectName))
        filename = "PhotPSaveData/" + slugify(settings.projectName) + "/detailedresults.csv"
    # Check the file's validity
    try:
        file = open(filename, "w")
    except FileNotFoundError:
        n = dirName.split("/")
        na = n[len(n) - 1]
        filename = "Output/" + na
        try:
            file = open(filename, "w")
        except FileNotFoundError:
            file = open("Output/INVALIDFILENAMERESULTS.csv", "w")
    # Output data
    # Make sure the results is an array
    try:
        len(results)
    except TypeError:
        results = [results]
    for i in range(len(results)):
        if not (results[i] == 0):
            file.write(
                "File Name,JD,Magnitude,Measurement Error,Target Net Counts,Target Radius (pixel),File Background Counts Per Pixel \n")
            trJD = truncate(results[i].JD, 6)
            trMag = truncate(results[i].magnitude, 2)
            trErr = truncate(results[i].error, 6)
            trCounts = truncate(results[i].targetStarCounts)
            trBack = truncate(results[i].fileBlanks)
            n = results[i].fileName.split("/")
            na = n[len(n) - 1]
            file.write(na + "," + str(trJD) + "," + str(trMag) + ","
                       + str(trErr) + "," + str(trCounts) + "," + str(results[i].targetStarRadius)
            + "," + str(trBack) + ", \n")
            stars = results[i].referenceStars
            file.write("ID,Right Ascension (Decimal Degrees),Declination (Decimal Degrees),Radius (pixels),Star Net Counts,Magnitude,Calculated Target Magnitude,Error Total/Net Counts,\n")
            try:
                len(stars)
            except TypeError:
                stars = [stars]
            # For each reference star, print all categories
            for j in range(len(stars)):
                trRA = truncate(stars[j].ra, 6)
                trDec = truncate(stars[j].dec, 6)
                trC = truncate(stars[j].counts)
                trTM = truncate(stars[j].targetMagnitude, 2)
                trE = truncate(stars[j].error)
                file.write(str(stars[j].id) + "," + str(trRA) + "," + str(trDec)
                           + "," + str(stars[j].radius) + "," + str(trC) +
                           "," + str(stars[j].magnitude) + "," + str(trTM) + "," + str(trE) + ", \n")
    trChi = truncate(reducedChiSquared(results), 2)
    file.write("X,Reduced Chi Square: " + str(trChi) + ", \n")
    file.close()

"""
NAME:       readFromFile
RETURNS:    An array of star objects
PARAMETERS: The name of the file to read from (fileName)
            The radius to use in calculating star counts (radius)
            Image array of the star (data)
            World coordinate system translation constant (w)
            Error star flag which indicates whether we are reading in an error star record (errorStarFlag)
PURPOSE:    This method reads in data on stars from a
            file to use to calculate the magnitude.
"""
def readFromFile(fileName, radius, data, w, errorStarFlag=0):
    global settings
    try:
        file = open(fileName, "r")
    except FileNotFoundError:
        if settings.consolePrintFlag == 1:
            print("Error: Invalid file for reading in reference stars")
        return 0
    stars = []
    for line in file:
        array = line.split(",")
        if not(array[0] == 'Name') and not(array[0] == '\n') and not(array[0] == '') and not(array[0] == 'ID') and errorStarFlag == 0 and not(array[0] == ' ID'):
            # Pull the star location from file
            X, Y = w.all_world2pix(float(array[1]), float(array[2]), 0)
            if settings.readInRadiusFlag == 0:
                # Calculate counts for this image
                starPhotons, e, extraStar = starCount(Y, X, data, radius)
                # Create a star object
                stars.append(Star(array[0], float(array[1]), float(array[2]), radius, starPhotons, float(array[5]), float(array[6]), e))
            elif settings.readInRadiusFlag == 1:
                # Calculate counts for this image
                starPhotons, e, extraStar = starCount(Y, X, data, float(array[3]))
                # Create a star object
                stars.append(Star(array[0], float(array[1]), float(array[2]), float(array[3]), starPhotons, float(array[5]),
                                  float(array[6]), e))
        elif not(array[0] == 'Name') and not(array[0] == '\n') and not(array[0] == '') and not(array[0] == 'ID') and errorStarFlag == 1:
            stars.append(Star(array[0], float(array[1]), float(array[2]), float(array[3]), 0, float(array[5]),
                              float(array[6]), 0))
    return stars

"""
NAME:       reducedChiSquared
RETURNS:    The reduced chi squared value of the data 
PARAMETERS: An array of Photometry objects (info)
PURPOSE:    Calculate the reduced chi squared value of
            a set of results in order to measure variability
"""
def reducedChiSquared(info):
    global settings
    # Make sure an array is passed in
    try:
        len(info)
    except TypeError:
        info = [info]
    if len(info) <= 1:
        # if invalid data is passed in, inform user and exit
        if settings.consolePrintFlag == 1:
            print("Reduced Chi Squared error: Reduced chi squared cannot be calculated for", len(info), "value(s).")
        return 0
    # Calculate the average magnitude of all items
    mBar = 0
    for i in range(len(info)):
        mBar = mBar + info[i].magnitude
        if info[i].error == 0:
            # if the error values are 0, reduced chi squared cannot be calculated
            if settings.consolePrintFlag == 1:
                print("Reduced Chi Squared error: Reduced chi squared calculation needs valid error values.")
            return 0
    mBar = mBar/len(info)
    # Calculate the reduced chi squared value
    chi = 0.0
    for i in range(len(info)):
        chi = chi + ((info[i].magnitude - mBar) / info[i].error) ** 2
        if info[i].error == 0:
            # if the error values are 0, reduced chi squared cannot be calculated
            if settings.consolePrintFlag == 1:
                print("Reduced Chi Squared error: Reduced chi squared calculation needs valid error values.")
            return 0
    chi = chi / (len(info)-1)
    return chi

"""
NAME:       printResultsToFile
RETURNS:    nothing
PARAMETERS: An array of Photometry objects (info)
            .csv file for output (filename)
PURPOSE:    Output a file in a new Output directory with the file name,
            the JD (Julian Date) of the observation, the magnitude 
            reported, and the error of that reported magnitude
"""
def printResultsToFile(info, filename="output.csv"):
    if settings.interface == 0:
        # Create new output directory if none exists
        if not os.path.exists("Output"):
            os.mkdir("Output")
        # Create new file
        filename = "Output/" + filename
    elif settings.interface == 1:
        if not os.path.exists("PhotPSaveData"):
            os.mkdir("PhotPSaveData")
        if not os.path.exists("PhotPSaveData/" + slugify(settings.projectName)):
            os.mkdir("PhotPSaveData/" + slugify(settings.projectName))
        filename = "PhotPSaveData/" + slugify(settings.projectName) + "/results.csv"
    try:
        file = open(filename, "w")
    except FileNotFoundError:
        n = filename.split("/")
        na = n[len(n) - 1]
        filename = "Output/" + na
        try:
            file = open(filename, "w")
        except FileNotFoundError:
            file = open("Output/INVALIDFILENAMERESULTS.csv", "w")
    file.write("File Name,JD,Magnitude,Error, \n")
    # Output data

    try:
        len(info)
    except TypeError:
        info = [info]
    for i in range(len(info)):
        if not(info[i] == 0):
            trJD = truncate(info[i].JD, 6)
            trMag = truncate(info[i].magnitude, 2)
            trErr = truncate(info[i].error, 6)
            n = info[i].fileName.split("/")
            na = n[len(n) - 1]
            file.write(na + "," + str(trJD) + "," + str(trMag) + ","
                       + str(trErr) + ", \n")
            """\X"""
    trChi = truncate(reducedChiSquared(info), 2)
    file.write("X,Reduced Chi Square: " + str(trChi) + ", \n")
    file.close()

"""
NAME:       plotResultsFile
RETURNS:    nothing
PARAMETERS: An input filename for results file (filename)
            An output chart file name (chartname)
            An output chart title (chartTitle)
PURPOSE:    Creates a light curve with error bars
"""
def plotResultsFile(filename, chartname="chart.png", chartTitle="Light Curve"):
    global settings
    # Create new output directory if none exists
    if not os.path.exists("Output"):
        os.mkdir("Output")
    # Read in output values from file
    try:
        file = open(filename, "r")
    except FileNotFoundError:
        if settings.consolePrintFlag == 1:
            print("Error: Invalid file for plotting results")
        return 0
    jd = []
    mag = []
    err = []
    for line in file:
        array = line.split(",")
        if not(array[0] == 'File Name') and not(array[0] == '\n') and not(array[0] == '') and not(array[0] == "X"):
            jd.append(float(array[1]))
            mag.append(float(array[2]))
            err.append(float(array[3]))
    # Smooth line print
    if settings.lightCurveLineFlag == 1:
        xnew = np.linspace(min(jd), max(jd), 300)
        spl = make_interp_spline(jd, mag, k=3)  # type: BSpline
        power_smooth = spl(xnew)
        plt.plot(xnew, power_smooth)
    # Chart Title
    plt.title(chartTitle)
    # Error bars
    plt.errorbar(jd, mag, err, fmt='ko')
    # X axis label
    plt.xlabel('Julian Day')
    # Y axis label
    plt.ylabel('Magnitude')
    # Inverting the y axis because a smaller magnitude is a brighter object
    plt.gca().invert_yaxis()
    if settings.interface == 0:
        if not os.path.exists("Output"):
            os.mkdir("Output")
        chartname = "Output/" + chartname
    elif settings.interface == 1:
        if not os.path.exists("PhotPSaveData"):
            os.mkdir("PhotPSaveData")
        if not os.path.exists("PhotPSaveData/" + slugify(settings.projectName)):
            os.mkdir("PhotPSaveData/" + slugify(settings.projectName))
        filename = "PhotPSaveData/" + slugify(settings.projectName) + "/lightcurve.png"
        plt.savefig(filename)
    if settings.printLightCurveFlag == 1:
        try:
            plt.savefig(chartname)  # to save to file
        except FileNotFoundError:
            n = chartname.split("/")
            na = n[len(n) - 1]
            chartname = "Output/" + na
            try:
                plt.savefig(chartname)
            except FileNotFoundError:
                plt.savefig("Output/INVALIDFILENAMECHART.png")
    if settings.showLightCurveFlag == 1:
        plt.show()  # to print to screen
    plt.clf()

"""
NAME:       plotResults
RETURNS:    nothing
PARAMETERS: An array of Photometry objects (ans)
            An output chart file name (chartname)
            An output chart title (chartTitle)
PURPOSE:    Creates a light curve with error bars
"""
def plotResults(ans, chartname="chart.png", chartTitle="Light Curve"):
    global settings
    # Create new output directory if none exists
    if not os.path.exists("Output"):
        os.mkdir("Output")
    # Arrange output values in arrays
    jd = []
    mag = []
    err = []
    try:
        len(ans)
    except TypeError:
        ans = [ans]
    for i in range(len(ans)):
        jd.append(ans[i].JD)
        mag.append(ans[i].magnitude)
        err.append(ans[i].error)
    # Smooth line print
    if settings.lightCurveLineFlag == 1 and len(ans) > 1:
        xnew = np.linspace(min(jd), max(jd), 300)
        spl = make_interp_spline(jd, mag, k=3)  # type: BSpline
        power_smooth = spl(xnew)
        plt.plot(xnew, power_smooth)
    # Chart Title
    plt.title(chartTitle)
    # Error bars
    plt.errorbar(jd, mag, err, fmt='ko')
    # X axis label
    plt.xlabel('Julian Day')
    # Y axis label
    plt.ylabel('Magnitude')
    # Inverting the y axis because a smaller magnitude is a brighter object
    plt.gca().invert_yaxis()
    if settings.interface == 0:
        if not os.path.exists("Output"):
            os.mkdir("Output")
        chartname = "Output/" + chartname
    elif settings.interface == 1:
        if not os.path.exists("PhotPSaveData"):
            os.mkdir("PhotPSaveData")
        if not os.path.exists("PhotPSaveData/" + slugify(settings.projectName)):
            os.mkdir("PhotPSaveData/" + slugify(settings.projectName))
        filename = "PhotPSaveData/" + slugify(settings.projectName) + "/lightcurve.png"
        plt.savefig(filename)
    if settings.printLightCurveFlag == 1:
        try:
            plt.savefig(chartname)  # to save to file
        except FileNotFoundError:
            n = chartname.split("/")
            na = n[len(n) - 1]
            chartname = "Output/" + na
            try:
                plt.savefig(chartname)
            except FileNotFoundError:
                plt.savefig("Output/INVALIDFILENAMECHART.png")
    if settings.showLightCurveFlag == 1:
        plt.show()  # to print to screen
    plt.clf()

"""
NAME:       calculateMagnitudeAndError
RETURNS:    Average magnitude, error of the magnitude measurement, and array of stars
            with valid magnitude calculations (i.e. visible stars)
PARAMETERS: A number values representing the counts in the target star (targetStarCounts)
            An array of Star objects representing reference stars (stars)
            A number representing the error in the target star counts (targetStarError)
            The error in the background pixel measurement (sigmaBkg)
                NOTE: The following are only needed for error star calculations
            Main photometry file (file)
            Dark frame (dark)
            Bias frame (bias)
            Flat frame (flat)
PURPOSE:    Calculates the average magnitude of the target star given the 
            reference stars, and the error in that measurement.
"""
def calculateMagnitudeAndError(targetStarPhotons, stars, targetStarError, sigmaBkg, file="", dark="", bias="", flat=""):
    global settings
    # Find the magnitude relative to each other star, and then averages them
    ave = 0
    aveCounts = 0
    aveMag = 0
    error = 0
    # Calculate average magnitude
    toRemove = []
    for i in range(len(stars)):
        if (targetStarPhotons > 0) and (stars[i].counts > 0):
            if settings.currentErrorStarsCalc == 1:
                stars[i].targetMagnitude = ((-2.5) * math.log10(targetStarPhotons / stars[i].counts))
                ave = ave + ((-2.5) * math.log10(targetStarPhotons / stars[i].counts))
            else:
                stars[i].targetMagnitude = ((-2.5) * math.log10(targetStarPhotons / stars[i].counts)
                                            + float(stars[i].magnitude))
                ave = ave + ((-2.5) * math.log10(targetStarPhotons / stars[i].counts)
                             + float(stars[i].magnitude))
            aveCounts = aveCounts + stars[i].counts
            aveMag = aveMag + stars[i].magnitude
        else:
            stars[i].targetMagnitude = 0
            toRemove.append(stars[i])
    for i in range(len(toRemove)):
        stars.remove(toRemove[i])
    if len(stars) == 0:
        # If the magnitude cannot be calculated for any reference star, exit
        if settings.consolePrintFlag == 1:
            print("Calculate magnitude error: Magnitude cannot be calculated for any reference star")
        return 0, 0, 0
    if settings.oneMagnitudeCalculation == 1:
        aveCounts = aveCounts / len(stars)
        aveMag = aveMag / len(stars)
        if settings.currentErrorStarsCalc == 1:
            ave = ((-2.5) * math.log10(targetStarPhotons / aveCounts))
        else:
            ave = ((-2.5) * math.log10(targetStarPhotons / aveCounts)
                                    + float(aveMag))
    else:
        ave = ave / len(stars)


    # Calculate  error
    w = []
    if settings.errorChoice == "STD":
        #Standard deviation as error
        for i in range(len(stars)):
            error = error + ((stars[i].targetMagnitude - ave) * (stars[i].targetMagnitude - ave))
        if len(stars) <= 0 or (len(stars)-1) <= 0:
            # If the magnitude cannot be calculated for any reference star, exit
            if settings.consolePrintFlag == 1:
                print("Calculate magnitude error: Magnitude cannot be calculated for any reference star")
            return 0, 0, 0
        error = error / (len(stars)-1)
        error = math.sqrt(error/len(stars))
    elif settings.errorChoice == "WMG":
        if settings.oneMagnitudeCalculation == 1:
            if settings.consolePrintFlag == 1:
                print("WARNING: oneMagnitudeCalculation setting and weighted magnitude error calculation incompatible. oneMagnitudeCalculation set to 0.")
            settings.oneMagnitudeCalculation = 0
        # Error for weighted magnitudes
        sigmaS = targetStarError
        sigmaS = math.sqrt(sigmaS**2 + sigmaBkg**2)
        for i in range(len(stars)):
            sigmaR = math.sqrt(stars[i].error**2 + sigmaBkg**2)
            sigmaMag = 0
            S = sigmaS**2 * ((2.5/np.log(10))/targetStarPhotons)**2
            R = sigmaR**2 * ((2.5/np.log(10)) / targetStarPhotons) **2
            M = sigmaMag**2
            sigmaI = math.sqrt(S + R + M)
            w.append(1/(sigmaI**2)) # Weights for the magnitudes

        # Calculate actual weighted magnitudes
        weightedAve = 0
        weightSum = 0
        magnitudeSum = 0
        for i in range(len(stars)):
            weightedAve = weightedAve + (stars[i].targetMagnitude * w[i])
            weightSum = weightSum + w[i]
        # Weighted magnitude
        ave = weightedAve/weightSum
        for i in range(len(stars)):
            magnitudeSum = magnitudeSum + ((stars[i].targetMagnitude - ave)**2)
        # Final Error
        if ((len(stars)-1) * weightSum) == 0:
            # If the magnitude cannot be calculated for any reference star, exit
            if settings.consolePrintFlag == 1:
                print("Calculate magnitude error: Magnitude cannot be calculated for any reference star")
            return 0, 0, 0
        error = math.sqrt(magnitudeSum/((len(stars)-1) * weightSum))
        #error = math.sqrt(1/weightSum)
    elif settings.errorChoice == "JKF":
        # Calculate Jack Knife error formula
        mM0 = []
        aveMM0 = 0
        for i in range(len(stars)):
            aveTemp = 0
            for j in range(len(stars)):
                if not(j == i):
                    aveTemp = aveTemp + ((-2.5) * math.log10(targetStarPhotons / stars[j].counts))#+ float(stars[j].magnitude))
            aveMM0 = aveMM0 + ((-2.5) * math.log10(targetStarPhotons / stars[i].counts))# + float(stars[i].magnitude))
            mM0.append(aveTemp/(len(stars)-1))
        aveMM0 = aveMM0/(len(mM0))
        error = 0
        for i in range(len(mM0)):
            error = error + ((mM0[i]-aveMM0)**2)
        #error = error/(len(mM0)-1)
        error = math.sqrt(((len(mM0)-1)/len(mM0))*error)
        #error = math.sqrt(error/len(mM0))
    elif settings.errorChoice == 0 or settings.currentErrorStarsCalc == 1:
        error = 0
    else:
        # Error stars method of error calculation
        # If there's no main file, exit
        if file == "":
            if settings.consolePrintFlag == 1:
                print("Calculate error problem: Error stars method of error requires main file to be passed into calculateMagnitudeAndError. Aborting error calculation.")
            return ave, 0, stars
        # Get error stars from the error choice, which is presumably a file at this point
        errorStars = readFromFile(settings.errorChoice, 0, 0, 0, errorStarFlag=1)
        # If it's not a file, exit
        if errorStars == 0:
            if settings.consolePrintFlag == 1:
                print("Calculate error problem: Setting errorChoice not set to a valid option or error star file. Aborting error calculation.")
            return ave, 0, stars
        count = 0
        error = 0
        # Make sure the errorStars is an array, otherwise return
        try:
            len(errorStars)
        except TypeError:
            if settings.consolePrintFlag == 1:
                print("No valid error stars input, returning error of 0.")
                return ave, 0, stars
        # For every error star, run photometry on it
        for i in range(len(errorStars)):
            settings.currentErrorStarsCalc = 1
            if settings.consolePrintFlag == 1:
                print("Starting error calculations for error star " + str(errorStars[i].id))
                settings.consolePrintFlag = -1
            ans = runPhotometry(errorStars[i].ra, errorStars[i].dec, file, darkFrame=dark, biasFrame=bias, flatField=flat)
            refAve = 0
            # Average magnitude of reference stars used to calculate the error
            for j in range(len(ans.referenceStars)):
                refAve = refAve + ans.referenceStars[j].magnitude
            refAve = refAve/len(ans.referenceStars)
            # If photometry was successful
            if not (ans == 0):
                delta_mag = ans.magnitude - (errorStars[i].magnitude - refAve)
                error = error + delta_mag
                errorStars[i].targetMagnitude = delta_mag
                errorStars[i].counts = ans.targetStarCounts
                count = count + 1
            settings.currentErrorStarsCalc = 0
            if settings.consolePrintFlag == -1:
                settings.consolePrintFlag = 1
        aveError = error/count
        error = 0
        for i in range(len(errorStars)):
            error = error + ((errorStars[i].targetMagnitude - aveError)*(errorStars[i].targetMagnitude - aveError))
        if count <= 0 or (count-1) <= 0:
            # If the magnitude cannot be calculated for any reference star, exit
            if settings.consolePrintFlag == 1:
                print("Calculate magnitude error: Magnitude cannot be calculated for any error stars. Returning error of zero.")
            return ave, 0, stars
        error = error / (count-1)
        error = math.sqrt(error)
        for i in range(len(errorStars)):
            errorStars[i].error = error
        # If user wants detailed info printed
        if settings.outputDetailedInformation == 1:
            if settings.interface == 0 and settings.interface == 0:
                # Create new output directory if none exists
                if not os.path.exists("Output"):
                    os.mkdir("Output")
                # Create new file
                n = file.split("/")
                na = n[len(n) - 1]
                if len(na) > 6:
                    filename = na[0:6] + "_detailederrorstarinformation.csv"
                else:
                    filename = na + "_detailederrorstarinformation.csv"
                filename = "Output/" + filename
                printReferenceToFile(errorStars, filename=filename, errorStars=1)
    return ave, error, stars

"""
NAME:       worldCoordinateSystem
RETURNS:    World coordinate system constant
PARAMETERS: A .fits container (hdul)
PURPOSE:    Uses the .fits data to return the world coordinate system.
"""
def worldCoordinateSystem(hdul):
    try:
        w = wcs.WCS(hdul[0].header)
    except TypeError:
        if settings.consolePrintFlag == 1:
            print("Invalid target for word coordinate system acquisition")
        return 0
    return w

"""
NAME:       getWCS
RETURNS:    New .fits file with WCS included
PARAMETERS: Name of the file to be submitted (file)
            Starting guess RA (ra)
            Starting guess declination (dec)
PURPOSE:    Add WCS header to .fits file without one using astrometry.net
"""
def getWCS(file, ra=0, dec=0):
    global settings
    if settings.astrometryDotNetFlag == 0:
        if settings.consolePrintFlag == 1:
            print("Astrometry.net unavailable: Please set astrometryDotNetFlag setting to a valid Astrometry.net API Key")
        return 0
    # Log in to astrometry.net
    if settings.consolePrintFlag == 1:
        print("Starting getWCS")
    astrometry = Client()
    try:
        astrometry.login(settings.astrometryDotNetFlag)
    except:
        if settings.consolePrintFlag == 1:
            print("Astrometry.net login error: Check if login is accurate and internet connection is secure")
        return 0
    if settings.consolePrintFlag == 1:
        print("Logged in")

    #Get previous jobid
    try:
        oldJobid = astrometry.myjobs()[0]
    except IndexError:
        oldJobid = 0

    # Upload file to astrometry.net
    if settings.consolePrintFlag == 1:
        print("Uploading File")
    if ra == 0 or dec == 0:
        response = astrometry.upload(fn=file)
    else:
        response = astrometry.upload(fn=file, center_ra=ra, center_dec=dec, radius=1)
    if settings.consolePrintFlag == 1:
        print("File Uploaded")

    splitFile = file.split('/')
    if settings.interface == 0:
        # Create new output directory if none exists
        if not os.path.exists("Output"):
            os.mkdir("Output")
        # Create filename to save to
        filename = "Output/CORRECTED_" + splitFile[len(splitFile) - 1]
    elif settings.interface == 1:
        if not os.path.exists("PhotPSaveData"):
            os.mkdir("PhotPSaveData")
        if not os.path.exists("PhotPSaveData/" + slugify(settings.projectName)):
            os.mkdir("PhotPSaveData/" + slugify(settings.projectName))
        if not os.path.exists("PhotPSaveData/" + slugify(settings.projectName) + "/Corrected"):
            os.mkdir("PhotPSaveData/" + slugify(settings.projectName) + "/Corrected")
        filename = "PhotPSaveData/" + slugify(settings.projectName) + "/Corrected" + "/CORRECTED_" + splitFile[len(splitFile) - 1]

    # Prevent program break from timeout
    try:
        # Get the job ID for uploaded file
        # It is compared to the old on to make sure it is the correct job ID
        try:
            jobid = astrometry.myjobs()[0]
        except ConnectionResetError:
            if settings.consolePrintFlag == 1:
                print("Astrometry.net submission job ID cannot be retrieved.")
            return 0
        jCheck = 0
        while oldJobid == jobid:
            jCheck = jCheck + 1
            try:
                jobid = astrometry.myjobs()[0]
            except ConnectionResetError:
                if settings.consolePrintFlag == 1:
                    print("Astrometry.net submission job ID cannot be retrieved.")
                return 0
            if settings.consolePrintFlag == 1:
              if jCheck%50 == 0:
                   print("Preparing job ID for", jCheck, "iterations. Please hold.")
            if jCheck > settings.astrometryTimeOutFlag and not(settings.astrometryTimeOutFlag == 0):
                if settings.consolePrintFlag == 1:
                    print("Astrometry.net submission job ID cannot be retrieved.")
                return 0
        if settings.consolePrintFlag == 1:
            print("Job ID retreived: ", jobid)

        # Request final product from website
        joburl = 'http://nova.astrometry.net/new_fits_file/' + str(jobid)
        if settings.consolePrintFlag == 1:
            print("File requested at ", joburl)

        # Wait for the file to finish downloading before progessing
        problem = 112
        while problem == 112:
            status = astrometry.job_status(jobid)
            solve = 0
            while(status == 'solving'):
                status = astrometry.job_status(jobid)
                if settings.consolePrintFlag == 1:
                    solve = solve + 1
                    if solve % 50 == 0:
                        print("File solving for", solve, "iterations. Please hold.")
            results = requests.get(joburl, allow_redirects=False)
            if astrometry.job_status(jobid) == 'failure':
                if settings.consolePrintFlag == 1:
                    print("Astrometry.net failed to solve submission.")
                return 0
            problem = results.content[0]
    except TimeoutError:
        if settings.consolePrintFlag == 1:
            print("Astrometry.net submission timed out.")
        return 0
    #Writing results to file
    file = open(filename, 'wb')
    file.write(results.content)
    file.close()
    if settings.consolePrintFlag == 1:
        print("File written to ", filename)
    return filename

"""
NAME:       removeReferenceOutliers
RETURNS:    Reference stars without outliers for a target star calculation
PARAMETERS: Array of reference stars (stars)
PURPOSE:    Remove outlier target magnitudes to create better results
"""
def removeReferenceOutliers(stars):
    global settings
    # This is sigma-clipping
    if settings.removeReferenceOutliersFlag == 0:
        return stars
    # Calculate Z-score
    tMags = []
    for i in range(len(stars)):
        tMags.append(stars[i].targetMagnitude)
    z = np.abs(stats.zscore(tMags))
    # Remove the stars that are higher than the specified z-score
    toRemoveStars = []
    for i in range(len(stars)):
        if z[i] >= settings.removeReferenceOutliersFlag:
            toRemoveStars.append(stars[i])
    for i in range(len(toRemoveStars)):
        stars.remove(toRemoveStars[i])
    return stars

"""
NAME:       runPhotometry
RETURNS:    A Photometry object with the data gained during 
            this process OR 0 if the magnitude cannot be
            computed
PARAMETERS: The right ascension of the target star in Decimal Degrees (targetStarRA)
            The declination of the target star in Decimal Degrees (targetStarDec)
            The main .fits file with the image data name (mainFile)
                NOTE: The calibration files are not necessary if calibration is turned off
            The dark frame .fits file name (darkFrame)
            The bias frame .fits file name (biasFrame)
            The flat field .fits file name (flatField)
PURPOSE:    This method combines the methods in this file to perform 
            full differential photometry on one image file. 
"""
def runPhotometry(targetStarRA, targetStarDec, mainFile, darkFrame="", biasFrame="", flatField=""):
    global settings
    if settings.calibrationFlag == 1:
        # Calibrate the image
        hdul = calibrate(mainFile, darkFrame, biasFrame, flatField)
    else:
        try:
            # Use the raw image
            hdul = fits.open(mainFile)
        except OSError:
            if settings.consolePrintFlag == 1:
                print("Error in runPhotometry: Invalid filename:", str(mainFile))
            return 0
    try:
        if hdul == 0:
            if settings.consolePrintFlag == 1:
                print("runPhotometry on", mainFile, "aborted due to calibration error")
            return 0
    except ValueError:
        hdul = hdul

    # Calculate magnitude
    # w is the reference of world coordinates for this image
    w = worldCoordinateSystem(hdul)

    # Convert COORDINATE system data into pixel locations for the image
    X, Y = w.all_world2pix(targetStarRA, targetStarDec, 0)

    # Converts pixel values to integers
    try:
        Y = int(Y)  # was pixRA
        X = int(X)  # was pixDec
    except ValueError:
        if settings.consolePrintFlag == 1:
            print("runPhotometry on", mainFile, "aborted due to invalid target star coordinates.")
        return 0

    # Check if the file has WCS data
    try:
        hdul[0].header['CRVAL1']
    except KeyError:
        if settings.astrometryDotNetFlag == 0:
            if settings.consolePrintFlag == 1:
                print("runPhotometry on", mainFile, "aborted due to localization error.")
                print("Suggestion: Set the astrometryDotNetFlag to a valid Astrometry.net API Key to localize files without world coordinate system data.")
            return 0
        else:
            # If it doesn't have WCS, add it
            mainFile = getWCS(mainFile, ra=targetStarRA, dec=targetStarDec)
            if mainFile == 0:
                if settings.consolePrintFlag == 1:
                    print("runPhotometry on", mainFile, "aborted due to localization error.")
                return 0
            if settings.calibrationFlag == 1:
                # Calibrate the image
                hdul = calibrate(mainFile, darkFrame, biasFrame, flatField)
            else:
                # Use the raw image
                hdul = fits.open(mainFile)

    # Calculate magnitude
    # w is the reference of world coordinates for this image
    w = worldCoordinateSystem(hdul)

    # Convert COORDINATE system data into pixel locations for the image
    X, Y = w.all_world2pix(targetStarRA, targetStarDec, 0)

    # Converts pixel values to integers
    Y = int(Y) #was pixRA
    X = int(X) #was pixDec

    if settings.fwhmFlag == 1:
        try:
            # Check if the file has FWHM
            radius = int(3 * hdul[0].header['FWHM'])
        except KeyError:
            # If not, find radius manually
            radius = findRadius(Y, X, hdul[0].data)
    elif settings.fwhmFlag == 0:
        # Set the radius to the distance from the center
        # of the star to the farthest edge of the star
        radius = findRadius(Y, X, hdul[0].data)
    else:
        radius = settings.fwhmFlag
    # Find the photon counts per pixel of blank sky
    settings.universalBlank = findBlank(hdul[0].data, radius)

    # Find the photon counts in the target star
    targetStarPhotons, targetStarError, extraStar = starCount(Y, X, hdul[0].data, radius)
    if extraStar == 1 and settings.consolePrintFlag == 1:
        print("WARNING: There may be another star near the target star.")

    # Find reference stars
    readInReferenceFilename = "0"
    if not(settings.readInReferenceFlag == 0):
        if not(isinstance(settings.readInReferenceFlag, str)):
            stars = settings.readInReferenceFlag
        elif settings.readInReferenceFlag[len(settings.readInReferenceFlag)-1] == 'v' and settings.readInReferenceFlag[len(settings.readInReferenceFlag)-2] == 's' and settings.readInReferenceFlag[len(settings.readInReferenceFlag)-3] == 'c':
            readInReferenceFilename = settings.readInReferenceFlag
        else:
            for filename in glob.glob(os.path.join(settings.readInReferenceFlag, '*.csv')):
                with open(os.path.join(os.getcwd(), filename), 'r') as f:
                    readInReferenceFilename = filename
        # Read in reference stars from file
        if readInReferenceFilename == "0" and isinstance(settings.readInReferenceFlag, str):
            stars = findOtherStars(Y, X, hdul[0].data, radius, w)
            if stars == 0:
                if settings.consolePrintFlag == 1:
                    print("runPhotometry on", mainFile, "aborted due to problem finding reference stars in field")
                return 0
        elif not(readInReferenceFilename == "0") and isinstance(settings.readInReferenceFlag, str):
            stars = readFromFile(readInReferenceFilename, radius, hdul[0].data, w)
            if stars == 0:
                if settings.consolePrintFlag == 1:
                    print("Reading in reference stars from file failed, finding stars automatically:")
                return 0
                stars = findOtherStars(Y, X, hdul[0].data, radius, w)
                if stars == 0:
                    if settings.consolePrintFlag == 1:
                        print("runPhotometry on", mainFile, "aborted due to problem finding reference stars in field")
                    return 0
    else:
        # Finding new stars automatically
        stars = findOtherStars(Y, X, hdul[0].data, radius, w)
        if stars == 0:
            if settings.consolePrintFlag == 1:
                print("runPhotometry on", mainFile, "aborted due to problem finding reference stars in field")
            return 0

    # Calculate magnitudes, average, and error
    a = radius * radius * math.pi
    sigmaBkg = math.sqrt(abs(a * settings.universalBlank))
    ave, error, stars = calculateMagnitudeAndError(targetStarPhotons, stars, targetStarError, sigmaBkg, file=mainFile, dark=darkFrame, bias=biasFrame, flat=flatField)
    if ave == 0 and error == 0:
        if settings.consolePrintFlag == 1:
            print("runPhotometry on", mainFile, "aborted due to problem calculating magnitude")
        return 0

    if not(settings.removeReferenceOutliersFlag == 0):
        firstAve = ave
        firstError = error
        # Remove outliers
        stars = removeReferenceOutliers(stars)

        # Recalculate average magnitude and standard deviation without outliers
        ave, error, stars = calculateMagnitudeAndError(targetStarPhotons, stars, targetStarError, sigmaBkg, file=mainFile, dark=darkFrame, bias=biasFrame, flat=flatField)
        if ave == 0 and error == 0:
            ave = firstAve
            error = firstError

    # Console output
    if settings.consolePrintFlag == 1:
        print("The magnitude of the star is ", ave )
        print("The error of this calculation is ", error)

    # Printing reference stars to files
    if settings.printReferenceStarsFlag == 1:
        n = mainFile.split("/")
        na = n[len(n)-1]
        printReferenceToFile(stars, filename=na[0:10]+"_referencestars.csv")

    # Create and return the results of the photometry
    ans = Photometry(mainFile, hdul[0].header['JD'], ave, error, stars, targetStarPhotons, radius, settings.universalBlank)
    return ans

"""
NAME:       getFiles
RETURNS:    Array of photometry objects containing names of files and their Julian date
            (other parameters left blank)
PARAMETERS: Directory to search (dirName)
PURPOSE:    Creates an array with the names of files and their Julian date
"""
def getFiles(dirName):
    global settings
    files = []
    # Check all files
    for filename in glob.glob(os.path.join(dirName, '*.fits')):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            # Open the files
            hdul = fits.open(filename)  # hdul is the computer version of
                                        # the file data
            files.append(Photometry(filename, hdul[0].header['JD'], hdul[0].header['EXPTIME'], 0, [], 0, 0, 0))
            hdul.close()
    for filename in glob.glob(os.path.join(dirName, '*.fts')):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            # Open the files
            hdul = fits.open(filename)  # hdul is the computer version of
                                        # the file data
            files.append(Photometry(filename, hdul[0].header['JD'], hdul[0].header['EXPTIME'], 0, [], 0, 0, 0))
            hdul.close()
    if len(files) <= 0:
        if settings.consolePrintFlag == 1:
            print("No files found in directory ", dirName)
        return 0
    return files

"""
NAME:       matchCal
RETURNS:    File name of matched calibration file
PARAMETERS: Name of the file to match a calibration file to (filename)
            Array of photometry objects containing calibration file names and dates (cals)
PURPOSE:    Match a data file to its closest calibration file by date
"""
def matchCal(filename, cals):
    global settings
    calName = "BLANK"
    possCals = []
    minTimeDiff = 10000000000000000000
    if len(cals) <= 0:
        if settings.consolePrintFlag == 1:
            print("No valid calibration files for match to ", str(filename))
        return 0
    try:
        hdul = fits.open(filename)
    except OSError:
        if settings.consolePrintFlag == 1:
            print("Error in matching calibration files: Invalid main .fits file ", str(filename))
        return 0
    currDate = hdul[0].header['JD']
    currExp = hdul[0].header['EXPTIME']
    hdul.close()
    # Find calibration file with smallest time difference from image file
    for i in range(len(cals)):
        if i == 0:
            minTimeDiff = abs(currDate - cals[i].JD)
        if abs(currDate - cals[i].JD) < minTimeDiff:
            minTimeDiff = abs(currDate - cals[i].JD)
            calName = cals[i].fileName
    # Check if multiple files on same day
    for i in range(len(cals)):
        if abs(minTimeDiff - abs(currDate - cals[i].JD)) < 1.5:
            temp = []
            temp.append(cals[i].fileName)
            temp.append(cals[i].magnitude)
            possCals.append(temp)
    # If only one calibration was taken on the same day, just return that one
    if (len(possCals)) <= 1:
        if calName == "BLANK":
            calName = cals[0].fileName
        return calName
    else:
        # Otherwise find the one with the closest exposure time
        minExpDiff = 100000
        calName = ""
        for i in range(len(possCals)):
            if abs(currExp - possCals[i][1]) < minExpDiff:
                minExpDiff = abs(currExp - possCals[i][1])
                calName = possCals[i][0]
        if calName == "BLANK":
            calName = cals[0].fileName
        return calName


"""
NAME:       runFiles
RETURNS:    Array of photometry objects containing results for photometry on every file
            in the directory
PARAMETERS: The right ascension of the target star in Decimal Degrees (targetStarRA)
            The declination of the target star in Decimal Degrees (targetStarDec)
            Directory containing files to process (dirName)
                NOTE: The following directories are only needed if calibration will be done
            Directory containing dark files (darkDirName)
            Directory containing bias files (biasDirName)
            Directory containing flat files (flatDirName)
PURPOSE:    Calculates the average magnitude of the target star given the 
            reference stars, and the error in that measurement by calculating
            standard deviation of the calculated magnitudes.
"""
def runFiles(targetStarRA, targetStarDec,
            dirName, darkDirName="", biasDirName="", flatDirName=""):
    global settings
    path = dirName
    # Make arrays of files and julian dates
    darkArray = getFiles(darkDirName)
    flatArray = getFiles(flatDirName)
    biasArray = getFiles(biasDirName)
    if darkArray == 0 and settings.consolePrintFlag == 1:
        print("No dark frame calibration files found")
    if flatArray == 0 and settings.consolePrintFlag == 1:
        print("No flat field calibration files found")
    if biasArray == 0 and settings.consolePrintFlag == 1:
        print("No bias frame calibration files found")
    if flatArray == 0:
        if settings.consolePrintFlag == 1:
            print("Aborting run files due to error finding calibration files")
        return 0
    if biasArray == 0:
        if settings.consolePrintFlag == 1:
            print("Running files without bias frame calibration.")
        settings.useBiasFlag = 0
        settings.subtractBiasFromDarkFlag = 0
    if darkArray == 0:
        if settings.consolePrintFlag == 1:
            print("Running files without dark frame calibration.")
        settings.useDarkFlag = 0
        settings.subtractBiasFromDarkFlag = 0

    results = []

    # Check and run all .fits files
    for filename in glob.glob(os.path.join(path, '*.fits')):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            # Get Correct calibration files:
            dark = matchCal(filename, darkArray)
            flat = matchCal(filename, flatArray)
            bias = matchCal(filename, biasArray)
            if not(dark == 0 and settings.useDarkFlag == 1) and not(bias == 0 and settings.useBiasFlag) and not(flat == 0):
                if settings.consolePrintFlag == 1:
                    print("Starting Photometry on Filename:", filename)
                    print("Chosen Calibration files:", dark, flat)
                x = runPhotometry(targetStarRA, targetStarDec, filename, dark, bias, flat)
                # If photometry done correctly:
                if x != 0:
                    if x.error < .1:
                        results.append(x)

    #Check and run all .fts files
    for filename in glob.glob(os.path.join(path, '*.fts')):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            # Get calibration files:
            dark = matchCal(filename, darkArray)
            flat = matchCal(filename, flatArray)
            bias = matchCal(filename, biasArray)
            if not(dark == 0 and settings.useDarkFlag == 1) and not(bias == 0 and settings.useBiasFlag) and not(flat == 0):
                if settings.consolePrintFlag == 1:
                    print("Starting Photometry on Filename:", filename)
                    print("Chosen Calibration files:", dark, bias, flat)
                x = runPhotometry(targetStarRA, targetStarDec, filename, dark, bias, flat)
                # If photometry done correctly:
                if x != 0:
                    results.append(x)

    # Check for outputting detailed information
    if not(settings.outputDetailedInformation == 0):
        printDetailedInfo(results, dirName=dirName)

    return results

"""
NAME:       createReferenceData
RETURNS:    Create light curves and Results files for individual reference star
            target magnitude calculations (i.e. one light curve for every reference star)
PARAMETERS: Photometry object containing all results (info)
            Name of the target star for file printing (targetStarName)
PURPOSE:    Allow the user to check each reference star for errors individually
            in a simple and visual way. 
"""
def createReferenceData(info, targetStarName="Star"):
    arr = []
    flag = 0
    arr.append([])
    arr[0].append(info[0].referenceStars[0].id)
    arr[0].append([])
    for i in range(len(info)):
        for j in range(len(info[i].referenceStars)):
            for k in range(len(arr)):
                if info[i].referenceStars[j].id == arr[k][0]:
                    flag = 1
            if flag == 0:
                arr.append([])
                arr[len(arr)-1].append(info[i].referenceStars[j].id)
                arr[len(arr) - 1].append([])
            flag = 0
    for i in range(len(info)):
        for j in range(len(info[i].referenceStars)):
            for k in range(len(arr)):
                if info[i].referenceStars[j].id == arr[k][0]:
                    arr[k][1].append(Photometry(info[i].fileName, info[i].JD, info[i].referenceStars[j].targetMagnitude, 0, [], 0, 0, 0))
    for i in range(len(arr)):
        ans = []
        for j in range(len(arr[i][1])):
            ans.append(arr[i][1][j])
        plotResults(ans, chartname=targetStarName+" Option A with "+str(arr[i][0])+".png", chartTitle=targetStarName+" Option A with "+str(arr[i][0]))
        printResultsToFile(ans, targetStarName+" Option A with "+str(arr[i][0])+".csv")

"""
NAME:       createCMD
RETURNS:    Two star object arrays representing the stars that were used in the plot
PARAMETERS: Array of star objects with V magnitudes (vstars)
            Array of star objects with B magnitudes (bstars)
            Name of the chart file to be printed (chartname)
            Name of the chart to be printed on the chart (charttitle)
PURPOSE:    To create a CMD plot of stars from a field.
"""
def createCMD(vstars, bstars, chartname="chart.png", charttitle="CMD"):
    global settings
    # Create arrays for V on y axis and B-V on x axis
    xAxis = []
    yAxis = []
    usedVStars = []
    usedBStars = []
    try:
        len(vstars)
    except TypeError:
        vstars = [vstars]
    try:
        vstars[0].ra
    except:
        if settings.consolePrintFlag == 1:
            print("CMD aborted for invalid V filter stars.")
        return 0,0
    try:
        len(bstars)
    except TypeError:
        bstars = [bstars]
    try:
        bstars[0].ra
    except:
        if settings.consolePrintFlag == 1:
            print("CMD aborted for invalid B filter stars.")
        return 0,0
    for i in range(len(vstars)):
        for j in range(len(bstars)):
            if abs(vstars[i].ra - bstars[j].ra) < .0001 and abs(vstars[i].dec - bstars[j].dec) < .0001:
                xAxis.append(abs(bstars[j].magnitude - vstars[i].magnitude))
                yAxis.append(vstars[i].magnitude)
                usedVStars.append(vstars[i])
                usedBStars.append(bstars[j])
    plt.plot(xAxis, yAxis, 'o', color='black')
    # Chart Title
    plt.title(charttitle)
    # X axis label
    plt.xlabel('B-V')
    # Y axis label
    plt.ylabel('V')
    # Inverting the y axis because a smaller magnitude is a brighter object
    plt.gca().invert_yaxis()
    if settings.printCMDFlag == 1:
        if settings.interface == 0:
            # Create new output directory if none exists
            if not os.path.exists("Output"):
                os.mkdir("Output")
            chartname = "Output/" + chartname
        elif settings.interface == 1:
            if not os.path.exists("PhotPSaveData"):
                os.mkdir("PhotPSaveData")
            if not os.path.exists("PhotPSaveData/" + slugify(settings.projectName)):
                os.mkdir("PhotPSaveData/" + slugify(settings.projectName))
            chartname = "PhotPSaveData/" + slugify(settings.projectName) + "/cmd.csv"
            plt.savefig(chartname)
        if settings.printLightCurveFlag == 1:
            try:
                plt.savefig(chartname)  # to save to file
            except FileNotFoundError:
                n = chartname.split("/")
                na = n[len(n) - 1]
                chartname = "Output/" + na
                try:
                    plt.savefig(chartname)
                except FileNotFoundError:
                    plt.savefig("Output/INVALIDFILENAMECHART.png")
    if settings.showCMDFlag == 1:
        plt.show()
    return usedVStars, usedBStars

"""
NAME:       findAllStars
RETURNS:    Array of Star objects representing all stars in the image
PARAMETERS: The main .fits file with the image data name (mainFile)
                NOTE: The following files are only required if calibration is to be performed
            The dark frame .fits file name (darkFrame)
            The bias frame .fits file name (biasFrame)
            The flat field .fits file name (flatField)
PURPOSE:    To find and calculate magnitudes of all stars in a field
"""
def findAllStars(mainFile, darkFrame="", biasFrame="", flatField=""):
    global settings
    # Set up the data for search
    hdul = calibrate(mainFile, darkFrame, biasFrame, flatField)
    settings.universalBlank = findBlank(hdul[0].data, int(hdul[0].header['FWHM'])*3)
    w = worldCoordinateSystem(hdul)
    # Find all stars in an image
    mean, median, std = sigma_clipped_stats(hdul[0].data, sigma=3.0)
    daofind = DAOStarFinder(fwhm=3, threshold=5.*std)
    sources = daofind(hdul[0].data)
    print(sources)
    ans = []
    stars = []
    # Get X, Y
    pixX = sources['xcentroid']
    pixY = sources['ycentroid']
    for i in range(len(pixX)):
        world = w.all_pix2world([[pixX[i], pixY[i]]], 0)
        ra = world[0][0]
        dec = world[0][1]
        # Get magnitude for found star
        x = runPhotometry(ra, dec, mainFile, darkFrame, biasFrame, flatField)
        # If magnitude gotten, create a star object
        if not(x == 0):
            ans.append(x)
            stars.append(Star("Ref"+str(i), ra, dec, 0, 0, ans[len(ans)-1].magnitude, 0, ans[len(ans)-1].error))
    return stars