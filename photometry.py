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
import json

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
        self.ra = ra #Pixel location
        self.dec = dec #Pixel location
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
    def __init__(self, fileName, JD, magnitude, error, referenceStars):
        self.fileName = fileName
        self.JD = JD
        self.magnitude = magnitude
        self.error = error
        self.referenceStars = referenceStars

class Settings:
    def __init__(self):
        # subtractBiasFromDarkFlag:
        # 0 = use dark and bias as-is
        # 1 = subtract bias from dark before calibrating
        #   NOTE: Dark files from the Great Basin Observatory have bias included and the bias
        #         will need to be subtracted out
        self.subtractBiasFromDarkFlag = 1

        # calibrationOutputFlag:
        # 0 = not outputting calibrated files
        # 1 = output calibrated files
        self.calibrationOutputFlag = 0

        # calibrationFlag:
        # 0 = not performing calibration in main methods
        # 1 = performing calibration in main methods
        self.calibrationFlag = 1

        # blankPerStarFlag:
        # 0 = calculate blank sky counts for the entire image
        # 1 = calculate blank sky counts for each star
        self.blankPerStarFlag = 0

        # catalogChoice:
        # II/336/apass9 = use APASS catalog
        # Vizier catalog designation = use user-chosen catalog
        # 0 = use SIMBAD
        #   WARNING: Using SIMBAD may result in errors being introduced into the calculation
        #   as different sources are used for magnitudes in SIMBAD
        self.catalogChoice = "II/336/apass9"

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
        self.filterChoice = "V"

        # lightCurveLineFlag:
        # 0 = generate light curve plots without lines (traditional)
        # 1 = generate light curve plots with smooth lines
        self.lightCurveLineFlag = 0

        # showLightCurveFlag:
        # 0 = do not print light curve plots to the screen
        # 1 = print light curve plots to the screen
        self.showLightCurveFlag = 0

        # errorChoice:
        # STD = use standard deviation method for error management
        # JKF = use jack knife method for error management
        # WMG = use weighted magnitude for error management
        # 0 = do not calculate error
        self.errorChoice = "STD"

        # consolePrintFlag:
        # 0 = do not print to console
        # 1 = print updates to console
        #   WARNING: printing to console may slightly slow down the program
        self.consolePrintFlag = 1

        # readInReferenceFlag:
        # 0 = do not read in reference stars from file
        # file name = use this reference file
        # directory name = detect reference file in directory
        self.readInReferenceFlag = 0

        # fwhmFlag:
        # 0 = detect radius for stars manually
        # 1 = use Full-Width Half-Maximum if available as radius for apertures
        self.fwhmFlag = 1

        # printReferenceStarsFlag:
        # 0 = do not print reference stars out
        # 1 = print out reference stars used
        self.printReferenceStarsFlag = 0

        # astrometryDotNetFlag:
        # 0 = do not use astrometry.net
        # astrometry.net API key = use astrometry.net with this API key
        self.astrometryDotNetFlag = 0

        # universalBlank:
        # 0 = blank not calculated
        # any number = median blank counts
        self.universalBlank = 0


settings = Settings()

def changeSettings(subtractBiasFromDarkFlag=-1, calibrationOutputFlag=-1,
                   calibrationFlag=-1, blankPerStarFlag=-1,
                   catalogChoice=-1, filterChoice=-1,
                   lightCurveLineFlag=-1, showLightCurveFlag=-1,
                   errorChoice=-1, consolePrintFlag=-1,
                   readInReferenceFlag=-1, fwhmFlag=-1,
                   printReferenceStarsFlag=-1, astrometryDotNetFlag=-1):
    global settings
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

"""
NAME:       calibrate
RETURNS:    calibrated data
PARAMETERS: the file to be calibrated (filename)
            the dark frame (dark)
            the flat frame (flat)
            the bias frame (bias)
PURPOSE:    To calibrate raw .fits files into a form that can 
            be used to calculate accurate magnitude data. 
"""
def calibrate(filename, dark, bias, flat):
    global settings
    # Open the files
    try:
        hdul = fits.open(filename)  # hdul is the computer version of
                                    # the file data
    except OSError:
        return 0
    #the following hdul files are computer readable versions of
    #the dark, and flat .fit files
    hdulDark = fits.open(dark)
    hdulFlat = fits.open(flat)
    hdulBias = fits.open(bias)

    # Extract data from the files
    data = hdul[0].data
    dataDark = hdulDark[0].data
    dataFlat = hdulFlat[0].data
    dataBias = hdulBias[0].data

    # Check if bias needs to be subtracted from dark
    if settings.subtractBiasFromDarkFlag == 1:
        dataDark = dataDark - dataBias

    # Calibrate the files, and then reassign the calibrated data to the
    # original data
    expTime = hdul[0].header['EXPTIME']/hdulDark[0].header['EXPTIME']
    data = data - (dataDark * expTime) #Subtract the dark frame from the data
    data = data - dataBias #Subtract the bias frame from the data

    dataFlatNorm = dataFlat / np.mean(dataFlat) #Normalize the flat frame information
    data = data // dataFlatNorm #Divide out the normalized flat frame data
    hdul[0].data = data #This sets the calibrated data to be a part of the
                        #hdul object

    # If outputting calibrated file, create folder and file
    if settings.calibrationOutputFlag == 1:
        if not os.path.exists("Output"):
            os.mkdir("Output")
        calibrationFile = "Output/CALIBRATED_" + filename
        hdul.writeto(calibrationFile, overwrite=True)
    return hdul

"""
NAME:       findCenter
RETURNS:    RA and DEC of the magnitude center of the star
PARAMETERS: Center of star pixel Y location (Y)
            Center of star pixel X location (X)
            Sky image data (data)
PURPOSE:    To take the recorded RA and DEC of a star and 
            center it for the image itself at the point of 
            highest photon counts
"""
def findCenter(Y, X, data):
    high = 0
    highY = Y
    highX = X
    #Search an area 100 pixels around the reported center
    for i in range(100):
        for j in range(100):
            #If the pixel being looked at has the highest
            #photon count so far, high should be set to that
            #amount and the RA and DEC should be set to that
            #location
            if data[X + 50 - j][Y + 50 - i] > high:
                high = data[X + 50 - j][Y + 50 - i]
                highY = Y + 50 - i
                highX = X + 50 - j
    #Return the true center of the star
    return highY, highX

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
    # Creates an aperture centered around the target star of radius r
    if r <= 0:
        r = 20
    targetAperture = CircularAperture((X, Y), r=r)
    targetStarTable = aperture_photometry(data, targetAperture)
    # Counts the sum in that aperture
    targetStarPhotons = targetStarTable['aperture_sum'][0]
    # Calculate Error in this count
    error = math.sqrt(targetStarPhotons)
    # If set to calculate per star
    if settings.blankPerStarFlag == 1 or settings.universalBlank == 0:
        blankAperture = CircularAnnulus((X, Y), r_in=r, r_out=r*4)
        blankTable = aperture_photometry(data, blankAperture)
        blankPhotons = blankTable['aperture_sum'][0]
        annulusArea = (math.pi * (r*4) * (r*4)) - (math.pi * r * r)
        blank = blankPhotons/annulusArea
    else:
        blank = settings.universalBlank

    # Subtracts blank counts per every pixel in the star
    targetStarPhotons = targetStarPhotons - ((math.pi * r * r) * blank)
    return targetStarPhotons, error

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
        result = Vizier.query_region(skyC, radius='0d10m0s', catalog=settings.catalogChoice)
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
                    # Calculate counts in the star
                    c, sigmaSrc = starCount(starY, starX, data, rad)
                    # Addition of this star to the array of stars objects
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
        result = customSimbad.query_region(skyC, radius='0d10m0s')  # This line creates warnings
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
                    # Calculate counts in the star
                    c, sigmaSrc = starCount(starY, starX, data, rad)
                    # Addition of this star to the array of stars objects
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
PURPOSE:    This method prints out a file  in a new directory named "Output"
            containing the stars used for calculating magnitude
"""
def printReferenceToFile(stars, filename="stars.csv"):
    # Create new output directory if none exists
    if not os.path.exists("Output"):
        os.mkdir("Output")
    # Create new file
    filename = "Output/" + filename
    file = open(filename, "w")
    # Create the heading
    file.write("ID,Right Ascension (Decimal Degrees),Declination (Decimal Degrees),Radius (pixels),Photons,Magnitude,Calculated Target Magnitude,Error,\n")
    # For each reference star, print all categories
    for i in range(len(stars)):
        trRA = truncate(stars[i].ra, 6)
        trDec = truncate(stars[i].dec, 6)
        trC = truncate(stars[i].counts)
        trTM = truncate(stars[i].targetMagnitude, 2)
        trE = truncate(stars[i].error)
        file.write(str(stars[i].id) + "," + str(trRA) + "," + str(trDec)
        + "," + str(stars[i].radius) + "," + str(trC) +
                   "," + str(stars[i].magnitude) + "," + str(trTM) + "," + str(trE) + ", \n")
    file.flush()
    file.close()

"""
NAME:       readFromFile
RETURNS:    An array of star objects
PARAMETERS: The name of the file to read from (fileName)
            The radius to use in calculating star counts (radius)
            Image array of the star (data)
            World coordinate system translation constant (w)
PURPOSE:    This method reads in data on stars from a
            file to use to calculate the magnitude.
"""
def readFromFile(fileName, radius, data, w):
    file = open(fileName, "r")
    stars = []
    for line in file:
        array = line.split(",")
        if not(array[0] == 'Name') and not(array[0] == '\n') and not(array[0] == '') and not(array[0] == 'ID'):
            # Pull the star location from file
            X, Y = w.all_world2pix(float(array[1]), float(array[2]), 0)
            # Calculate counts for this image
            starPhotons, e = starCount(Y, X, data, radius)
            # Create a star object
            stars.append(Star(array[0], float(array[1]), float(array[2]), radius, starPhotons, float(array[5]), float(array[6]), e))
    return stars

"""
NAME:       reducedChiSquared
RETURNS:    The reduced chi squared value of the data 
PARAMETERS: An array of Photometry objects (info)
PURPOSE:    Calculate the reduced chi squared value of
            a set of results in order to measure variability
"""
def reducedChiSquared(info):
    if len(info) == 1:
        return 0
    mBar = 0
    for i in range(len(info)):
        mBar = mBar + info[i].magnitude
        if info[i].error == 0:
            return 0
    mBar = mBar/len(info)
    chi = 0.0
    for i in range(len(info)):
        chi = chi + ((info[i].magnitude - mBar) / info[i].error) ** 2
        if info[i].error == 0:
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
    # Create new output directory if none exists
    if not os.path.exists("Output"):
        os.mkdir("Output")
    # Create new file
    filename = "Output/" + filename
    file = open(filename, "w")
    file.write("File Name,JD,Magnitude,Error, \n")
    # Output data)
    for i in range(len(info)):
        trJD = truncate(info[i].JD, 6)
        trMag = truncate(info[i].magnitude, 2)
        trErr = truncate(info[i].error, 6)
        file.write(info[i].fileName + "," + str(trJD) + "," + str(trMag) + ","
                   + str(trErr) + ", \n")
    trChi = truncate(reducedChiSquared(info), 2)
    file.write("X,Reduced Chi Square: " + str(trChi) + ", \n")
    file.close()

"""
NAME:       plotResultsFile
RETURNS:    nothing
PARAMETERS: An input filename for results file (filename)
            An output chartname (chartname)
            An output chart title (chartTitle)
PURPOSE:    Creates a light curve with error bars
"""
def plotResultsFile(filename, chartname="chart.pdf", chartTitle="Light Curve"):
    global settings
    # Create new output directory if none exists
    if not os.path.exists("Output"):
        os.mkdir("Output")
    # Read in output values from file
    file = open(filename, "r")
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
    chartname = "Output/" + chartname
    plt.savefig(chartname)  # to save to file
    if settings.showLightCurveFlag == 1:
        plt.show()  # to print to screen
    plt.clf()

"""
NAME:       plotResults
RETURNS:    nothing
PARAMETERS: An array of Photometry objects (ans)
            An output chartname (chartname)
            An output chart title (chartTitle)
PURPOSE:    Creates a light curve with error bars
"""
def plotResults(ans, chartname="chart.pdf", chartTitle="Light Curve"):
    global settings
    # Create new output directory if none exists
    if not os.path.exists("Output"):
        os.mkdir("Output")
    # Arrange output values in arrays
    jd = []
    mag = []
    err = []
    for i in range(len(ans)):
        jd.append(ans[i].JD)
        mag.append(ans[i].magnitude)
        err.append(ans[i].error)
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
    chartname = "Output/" + chartname
    plt.savefig(chartname)  # to save to file
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
PURPOSE:    Calculates the average magnitude of the target star given the 
            reference stars, and the error in that measurement.
"""
def calculateMagnitudeAndError(targetStarPhotons, stars, targetStarError, sigmaBkg):
    global settings
    # Find the magnitude relative to each other star, and then averages them
    ave = 0
    error = 0
    # Calculate average magnitude
    toRemove = []
    for i in range(len(stars)):
        if (targetStarPhotons > 0) and (stars[i].counts > 0):
            stars[i].targetMagnitude = ((-2.5) * math.log10(targetStarPhotons / stars[i].counts)
                              + float(stars[i].magnitude))
            ave = ave + ((-2.5) * math.log10(targetStarPhotons / stars[i].counts)
                         + float(stars[i].magnitude))
        else:
            stars[i].targetMagnitude = 0
            toRemove.append(stars[i])
    for i in range(len(toRemove)):
        stars.remove(toRemove[i])
    if len(stars) == 0:
        # If the magnitude cannot be calculated for any reference star, exit
        return 0, 0, 0
    ave = ave / len(stars)

    # Calculate  error
    w = []
    if settings.errorChoice == "STD":
        #Standard deviation as error
        for i in range(len(stars)):
            error = error + ((stars[i].targetMagnitude - ave) * (stars[i].targetMagnitude - ave))
        if len(stars) <= 0 or (len(stars)-1) <= 0:
            # If the magnitude cannot be calculated for any reference star, exit
            return 0, 0, 0
        error = error / (len(stars)-1)
        error = math.sqrt(error/len(stars))
    elif settings.errorChoice == "WMG":
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
    else:
        error = 0
    return ave, error, stars

"""
NAME:       worldCoordinateSystem
RETURNS:    World coordinate system constant
PARAMETERS: A .fits container (hdul)
PURPOSE:    Uses the .fits data to return the world coordinate system.
"""
def worldCoordinateSystem(hdul):
    w = wcs.WCS(hdul[0].header)
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
    # Log in to astrometry.net
    if settings.consolePrintFlag == 1:
        print("Starting getWCS")
    astrometry = Client()
    astrometry.login(settings.astrometryDotNetFlag)
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

    # Create new output directory if none exists
    if not os.path.exists("Output"):
        os.mkdir("Output")
    # Create filename to save to
    splitFile = file.split('/')
    filename = "Output/CORRECTED_" + splitFile[len(splitFile) - 1]

    # Prevent program break from timeout
    try:
        # Get the job ID for uploaded file
        # It is compared to the old on to make sure it is the correct job ID
        jobid = astrometry.myjobs()[0]
        jCheck = 0
        while oldJobid == jobid:
            jCheck = jCheck + 1
            jobid = astrometry.myjobs()[0]
            if settings.consolePrintFlag == 1:
              if jCheck%50 == 0:
                   print("Preparing job ID for", jCheck, "iterations. Please hold.")
            if jCheck > 1500:
                if settings.consolePrintFlag == 1:
                    print("Job failed.")
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
                    print("Job failed.")
                return 0
            problem = results.content[0]
    except TimeoutError:
        if settings.consolePrintFlag == 1:
            print("Job timed out.")
        return 0
    #Writing results to file
    file = open(filename, 'wb')
    file.write(results.content)
    file.close()
    if settings.consolePrintFlag == 1:
        print("File written to ", filename)
    return filename
    return 0

"""
NAME:       letsGo
RETURNS:    A Photometry object with the data gained during 
            this process OR 0 if the magnitude cannot be
            computed
PARAMETERS: The right ascension of the target star in Decimal Degrees (targetStarRA)
            The declination of the target star in Decimal Degrees (targetStarDec)
            The main .fits file with the image data name (mainFile)
            The dark frame .fits file name (darkFrame)
            The flat field .fits file name (flatField)
PURPOSE:    This method combines the methods in this file to perform 
            full differential photometry on one image file. 
"""
def letsGo(targetStarRA, targetStarDec, mainFile, darkFrame, biasFrame, flatField):
    global settings
    if settings.calibrationFlag == 1:
        # Calibrate the image
        hdul = calibrate(mainFile, darkFrame, biasFrame, flatField)
    else:
        # Use the raw image
        hdul = fits.open(mainFile)
    if hdul == 0:
        return 0

    # Calculate magnitude
    # w is the reference of world coordinates for this image
    w = worldCoordinateSystem(hdul)

    # Convert COORDINATE system data into pixel locations for the image
    X, Y = w.all_world2pix(targetStarRA, targetStarDec, 0)

    # Converts pixel values to integers
    Y = int(Y)  # was pixRA
    X = int(X)  # was pixDec

    # Check if the file has WCS data
    try:
        hdul[0].header['CRVAL1']
    except KeyError:
        if settings.astrometryDotNetFlag == 0:
            return 0
        else:
            # If it doesn't have WCS, add it
            mainFile = getWCS(mainFile, ra=targetStarRA, dec=targetStarDec)
            if mainFile == 0:
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

    # Center the star (IN PROGRESS)
    #pixRA, pixDec = findCenter(int(Y), int(X), hdul[0].data)

    if settings.fwhmFlag == 1:
        try:
            # Check if the file has FWHM
            radius = int(3 * hdul[0].header['FWHM'])
        except KeyError:
            # If not, find radius manually
            radius = findRadius(Y, X, hdul[0].data)
    else:
        # Set the radius to the distance from the center
        # of the star to the farthest edge of the star
        radius = findRadius(Y, X, hdul[0].data)
    # Find the photon counts per pixel of blank sky
    settings.universalBlank = findBlank(hdul[0].data, radius)

    # Find the photon counts in the target star
    targetStarPhotons, targetStarError = starCount(Y, X, hdul[0].data, radius)

    # Find reference stars
    readInReferenceFilename = "0"
    if not(settings.readInReferenceFlag == 0):
        if settings.readInReferenceFlag[len(settings.readInReferenceFlag)-1] == 'v' and settings.readInReferenceFlag[len(settings.readInReferenceFlag)-2] == 's' and settings.readInReferenceFlag[len(settings.readInReferenceFlag)-3] == 'c':
            readInReferenceFilename = settings.readInReferenceFlag
        else:
            for filename in glob.glob(os.path.join(settings.readInReferenceFlag, '*.csv')):
                with open(os.path.join(os.getcwd(), filename), 'r') as f:
                    readInReferenceFilename = filename
        # Read in reference stars from file
        if readInReferenceFilename == "0":
            stars = findOtherStars(Y, X, hdul[0].data, radius, w)
        else:
            stars = readFromFile(readInReferenceFilename, radius, hdul[0].data, w)
    else:
        # Finding new stars automatically
        stars = findOtherStars(Y, X, hdul[0].data, radius, w)

    # Calculate magnitudes, average, and error
    a = radius * radius * math.pi
    sigmaBkg = math.sqrt(a * settings.universalBlank)
    ave, error, stars = calculateMagnitudeAndError(targetStarPhotons, stars, targetStarError, sigmaBkg)
    if ave == 0 and error == 0:
        print("Problem with first calculation.")
        return 0

    # Remove outliers
    if error >= (ave/20):
        print("I calculate there may be some outliers in the data. Review this list "
              + "below for outliers: ")
        for i in range(len(stars)):
            print(str(i+1) + ") " + str(stars[i].targetMagnitude))
        print("The calculated average magnitude is " + str(ave) +" and the calculated error is "
              + str(error) + ".")
        print("\nPlease input the number next to the magnitudes you want to remove from the calculations: "
              + "(enter 100 if there are no more outliers to remove) ")
        x = input()
        x = int(x)
        toRemoveStars = []
        while x != 100:
            toRemoveStars.append(stars[x-1])
            x = input()
            x = int(x)
        for i in range(len(toRemoveStars)):
            stars.remove(toRemoveStars[i])

        # Recalculate average magnitude and standard deviation without outliers
        ave, error, stars = calculateMagnitudeAndError(targetStarPhotons, stars, targetStarError, sigmaBkg)
        if ave == 0 and error == 0:
            print("Problem with second calculation.")
            return 0

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
    ans = Photometry(mainFile, hdul[0].header['JD'], ave, error, stars)
    return ans

"""
NAME:       getFiles
RETURNS:    Array of photometry objects containing names of files and their Julian date
            (other parameters left blank)
PARAMETERS: Directory to search (dirName)
PURPOSE:    Creates an array with the names of files and their Julian date
"""
def getFiles(dirName):
    files = []
    # Check all files
    for filename in glob.glob(os.path.join(dirName, '*.fits')):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            # Open the files
            hdul = fits.open(filename)  # hdul is the computer version of
                                        # the file data
            files.append(Photometry(filename, hdul[0].header['JD'], hdul[0].header['EXPTIME'], 0, []))
            hdul.close()
    for filename in glob.glob(os.path.join(dirName, '*.fts')):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            # Open the files
            hdul = fits.open(filename)  # hdul is the computer version of
                                        # the file data
            files.append(Photometry(filename, hdul[0].header['JD'], hdul[0].header['EXPTIME'], 0, []))
            hdul.close()
    return files

"""
NAME:       matchCal
RETURNS:    File name of matched calibration file
PARAMETERS: Name of the file to match a calibration file to (filename)
            Array of photometry objects containing calibration file names and dates (cals)
PURPOSE:    Match a data file to its closest calibration file by date
"""
def matchCal(filename, cals):
    calName = ""
    possCals = []
    minTimeDiff = 100000000000
    hdul = fits.open(filename)
    currDate = hdul[0].header['JD']
    currExp = hdul[0].header['EXPTIME']
    hdul.close()
    # Find calibration file with smallest time difference from image file
    for i in range(len(cals)):
        if abs(currDate - cals[i].JD) < minTimeDiff:
            minTimeDiff = abs(currDate - cals[i].JD)
            calName = cals[i].fileName
        if abs(currDate - cals[i].JD) < 1.5:
            # If there are multiple files on same day as the target image
            # collect them to measure by exposure
            temp = []
            temp.append(cals[i].fileName)
            temp.append(cals[i].magnitude)
            possCals.append(temp)
    # If only one calibration was taken on the same day, just return that one
    if (len(possCals)) <= 1:
        return calName
    else:
        # Otherwise find the one with the closest exposure time
        minExpDiff = 100000
        calName = ""
        for i in range(len(possCals)):
            if abs(currExp - possCals[i][1]) < minExpDiff:
                minExpDiff = abs(currExp - possCals[i][1])
                calName = possCals[i][0]
        return calName


"""
NAME:       runFiles
RETURNS:    Array of photometry objects containing results for photometry on every file
            in the directory
PARAMETERS: The right ascension of the target star in Decimal Degrees (targetStarRA)
            The declination of the target star in Decimal Degrees (targetStarDec)
            The apikey for use with astrometry.net (key)
            Directory containing files to process (dirName)
            Directory containing dark files (darkDirName)
            Directory containing flat files (flatDirName)
PURPOSE:    Calculates the average magnitude of the target star given the 
            reference stars, and the error in that measurement by calculating
            standard deviation of the calculated magnitudes.
"""
def runFiles(targetStarRA, targetStarDec,
            dirName, darkDirName, biasDirName, flatDirName):
    global settings
    path = dirName
    # Make arrays of files and julian dates
    darkArray = getFiles(darkDirName)
    flatArray = getFiles(flatDirName)
    biasArray = getFiles(biasDirName)

    results = []

    # Check and run all .fits files
    for filename in glob.glob(os.path.join(path, '*.fits')):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            # Get calibration files:
            dark = matchCal(filename, darkArray)
            flat = matchCal(filename, flatArray)
            bias = matchCal(filename, biasArray)

            print("Filename:", filename)
            print("Chosen Calibration files:", dark, flat)
            x = letsGo(targetStarRA, targetStarDec, filename, dark)
            if x != 0:
                results.append(x)

    #Check and run all .fts files
    for filename in glob.glob(os.path.join(path, '*.fts')):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            # Get calibration files:
            dark = matchCal(filename, darkArray)
            flat = matchCal(filename, flatArray)
            #bias = matchCal(filename, biasArray)

            x = letsGo(targetStarRA, targetStarDec, filename, dark, flat)
            if x != 0:
                results.append(x)
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
                    arr[k][1].append(Photometry(info[i].fileName, info[i].JD, info[i].referenceStars[j].targetMagnitude, 0, []))
    for i in range(len(arr)):
        ans = []
        for j in range(len(arr[i][1])):
            ans.append(arr[i][1][j])
        plotResults(ans, chartname="July 14 "+targetStarName+" Option A with "+str(arr[i][0])+".pdf", chartTitle="July 14 "+targetStarName+" Option A with "+str(arr[i][0]))
        printResultsToFile(ans, targetStarName+" Option A with "+str(arr[i][0])+".csv")