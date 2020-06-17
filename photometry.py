from __future__ import division
from __future__ import print_function
from astroquery.simbad import Simbad
from scipy.interpolate import make_interp_spline, BSpline

from client import Client

from photutils import CircularAperture
from photutils import Background2D
from astropy.io import fits
from astropy import wcs
from photutils import aperture_photometry
import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import glob
import os
import json
import requests

import sys
import time
import base64
import optparse
import threading
import io
from email.generator import Generator

from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.application  import MIMEApplication
from email.encoders import encode_noop
import json

"""
Star objects contain everything that makes a star a star, including
the name of the star, its location in decimal degrees, the radius 
being used for the star, the counts calculated over that radius,
the magnitude of the star and the magnitude calculated for the 
target star using this star.
"""
class Star:
    def __init__(self, id, ra, dec, radius, counts, magnitude, targetMagnitude):
        self.id = id
        self.ra = ra #Pixel location
        self.dec = dec #Pixel location
        self.radius = radius
        self.counts = counts
        self.magnitude = magnitude
        self.targetMagnitude = targetMagnitude

"""
Photometry objects contain the results of differential photometry.
It includes the name of the file being operated on, the julian
date of the observation, the magnitude of the target star, and the
error of the magnitude calculation. 
"""
class Photometry:
    def __init__(self, fileName, JD, magnitude, error):
        self.fileName = fileName
        self.JD = JD
        self.magnitude = magnitude
        self.error = error

"""
NAME:       calibrate
RETURNS:    calibrated data
PARAMETERS: the file to be calibrated (filename)
            the dark frame (dark)
            the flat frame (flat)
            the bias frame (bias)
                NOTE: If bias not provided, the bias is not used
PURPOSE:    to calibrate raw .fits files into a form that can 
            be used to calculate accurate magnitude data. 
            This currently only uses dark and flat but can 
            easily be modified to use bias as well
"""
def calibrate(filename, dark, flat, bias=0):
    # Open the files
    hdul = fits.open(filename)  # hdul is the computer version of
                                # the file data
    #the following hdul files are computer readable versions of
    #the dark, and flat .fit files
    hdulDark = fits.open(dark)
    hdulFlat = fits.open(flat)

    # Extract data from the files
    data = hdul[0].data
    dataDark = hdulDark[0].data
    dataFlat = hdulFlat[0].data
    x, y = data.shape

    # Calibrate the files, and then reassign the calibrated data to the
    # original data
    data = data - dataDark #Subtract the dark frame from the data

    #Bias calculations if bias is specified
    if not(bias == 0):
        hdulBias = fits.open(bias)
        dataBias = hdulBias[0].data
        data = data - dataBias

    dataFlatNorm = dataFlat / np.mean(dataFlat) #Normalize the flat frame information
    data = data // dataFlatNorm #Divide out the normalized flat frame data
    hdul[0].data = data #This sets the calibrated data to be a part of the
                        #hdul object
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

    x, y = data.shape
    currY = Y
    currX = X
    #Blank is the number of counts per pixel in empty sky
    blank = findBlankNoRad(data) * 1.5
    r1 = 0
    r2 = 0
    r3 = 0
    r4 = 0
    #Find the distance it takes each direction to get to blank sky
    while data[currX][currY] > blank:
        currY = currY + 1
        if currY >= y:
            currY = Y
            break
    r1 = currY - Y
    currY = Y
    while data[currX][currY] > blank:
        currY = currY - 1
        if currY <= 1:
            currY = Y
            break
    r2 = -1 * (currY - Y)
    currY = Y
    while data[currX][currY] > blank:
        currX = currX - 1
        if currX <= 1:
            currX = X * -1
            break
    r3 = -1 * (currX - X)
    currX = X
    while data[currX][currY] > blank:
        currX = currX + 1
        if currX >= x:
            currX = X
            break
    r4 = currX - X

    # Take the maximum out of the four found radius
    max = 100
    if r1 >= r2 and r1 >= r3 and r1 >= r4:
        max = r1
    elif r2 >= r1 and r2 >= r3 and r2 >= r4:
        max = r2
    elif r3 >= r1 and r3 >= r2 and r3 >= r4:
        max = r3
    elif r4 >= r1 and r4 >= r2 and r4 >= r3:
        max = r4

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
    background = Background2D(data, ((r * 4), (r * 4)))
    return background.background_median

"""
NAME:       starCount
RETURNS:    The amount of photons in a star's image
PARAMETERS: Center of star pixel Y location (Y)
            Center of star pixel X location (X)
            Image data (data)
            Radius of the star (r)
            Average count per blank sky pixel (blank)
PURPOSE:    Taking in information about any star, this method counts
            the photons in a circular area around the center of the star.
"""
def starCount(Y, X, data, r, blank):
    # Creates an aperture centered around the target star of radius r
    targetAperture = CircularAperture((X, Y), r=r)
    targetStarTable = aperture_photometry(data, targetAperture)
    # Counts the sum in that aperture
    targetStarPhotons = targetStarTable['aperture_sum'][0]
    # Subtracts blank counts per every pixel in the star
    targetStarPhotons = targetStarPhotons - ((math.pi * r * r) * blank)
    return targetStarPhotons

"""
NAME:       findOtherStars
RETURNS:    Array of other star objects
PARAMETERS: Center of target star RA (Y)
            Center of target star DEC (X)
            Image data (data)
            Radius to examine (rad)
            Counts in blank sky (blank)
            Parameters for conversion from world coordinates to
            pixel location (w)
PURPOSE:    This method locates stars nearby the target stars
            and then calls a helper method to construct a star
            object for the located stars.
"""
def findOtherStars(Y, X, data, rad, blank, w):
    stars = []
    skyC = SkyCoord.from_pixel(X, Y, w)
    # Query all stars in the image
    customSimbad = Simbad()
    # Sets the flux being looked for to V magnitudes
    customSimbad.add_votable_fields('flux(V)')
    # Query in an area about the size of the image from GBO
    result = customSimbad.query_region(skyC, radius='0d13m0s')  # This line creates warnings
    # If there are any results:
    if not (result is None):
        name = result['MAIN_ID']
        ra = result['RA']
        dec = result['DEC']
        mag = result['FLUX_V']
        x = 0
        for i in range(len(mag)):
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
                # Addition of this star to the array of stars objects
                stars.append(Star(name[i], degreeRa, degreeDec, rad, starCount(starY, starX, data, rad, blank), mag[i], 0))
                x = x+1
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
    file.write("Name,Right Ascension (Decimal Degrees),Declination (Decimal Degrees),Radius (pixels),Photons,Magnitude,\n")
    # For each reference star, print all categories
    for i in range(len(stars)):
        trRA = truncate(stars[i].ra, 6)
        trDec = truncate(stars[i].dec, 6)
        trC = truncate(stars[i].counts)
        file.write(str(stars[i].id) + "," + str(trRA) + "," + str(trDec)
        + "," + str(stars[i].radius) + "," + str(trC) +
                   "," + str(stars[i].magnitude) + ", \n")
    file.close()

"""
NAME:       readFromFile
RETURNS:    An array of star objects
PARAMETERS: The name of the file to read from (fileName)
            The radius to use in calculating star counts (radius)
            Image array of the star (data)
            Average blank counts in sky (blank)
            World coordinate system translation constant (w)
PURPOSE:    This method reads in data on stars from a
            file to use to calculate the magnitude.
"""
def readFromFile(fileName, radius, data, blank, w):
    file = open(fileName, "r")
    stars = []
    for line in file:
        array = line.split(",")
        if not(array[0] == 'Name') and not(array[0] == '\n') and not(array[0] == ''):
            # Pull the star location from file
            X, Y = w.all_world2pix(float(array[1]), float(array[2]), 0)
            # Calculate counts for this image
            starPhotons = starCount(Y, X, data, radius, blank)
            # Create a star object
            stars.append(Star(array[0], float(array[1]), float(array[2]), radius, starPhotons, float(array[5]), array[6]))
    return stars

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
        trErr = truncate(info[i].error)
        file.write(info[i].fileName + "," + str(trJD) + "," + str(trMag) + ","
                   + str(trErr) + ", \n")
    file.close()

"""
NAME:       plotResultsFile
RETURNS:    nothing
PARAMETERS: An output filename (filename)
            An output chartname (chartname)
            An output chart title (chartTitle)
            lineFlag:
                Set to 1 to plot a smooth line 
                over the points
PURPOSE:    Creates a light curve with error bars
"""
def plotResultsFile(filename, chartname="chart.pdf", chartTitle="Light Curve", lineFlag=0):
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
        if not(array[0] == 'File Name') and not(array[0] == '\n') and not(array[0] == ''):
            jd.append(float(array[1]))
            mag.append(float(array[2]))
            err.append(float(array[3]))
    # Smooth line print
    if lineFlag == 1:
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
    plt.show()  # to print to screen
    chartname = "Output/" + chartname
    plt.savefig(chartname)  # to save to file

"""
NAME:       calculateMagnitudeAndError
RETURNS:    Average magnitude, error of the magnitude measurement, and array of stars
            with valid magnitude calculations (i.e. visible stars)
PARAMETERS: A number values representing the counts in the target star (targetStarCounts)
            An array of Star objects representing reference stars (stars)
PURPOSE:    Calculates the average magnitude of the target star given the 
            reference stars, and the error in that measurement by calculating
            standard deviation of the calculated magnitudes.
"""
def calculateMagnitudeAndError(targetStarPhotons, stars):
    # Find the magnitude relative to each other star, and then averages them
    ave = 0
    std = 0
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
        return 0
    ave = ave / len(stars)

    # Calculate standard deviation of magnitudes for error
    for i in range(len(stars)):
        std = std + ((stars[i].targetMagnitude - ave) * (stars[i].targetMagnitude - ave))
    std = std / len(stars)
    std = math.sqrt(std)
    return ave, std, stars

"""
NAME:       worldCoordinateSystem
RETURNS:    World coordinate system constant
PARAMETERS: An unwrapped .fits container (hdul)
PURPOSE:    Uses the .fits data to return the world coordinate system.
"""
def worldCoordinateSystem(hdul):
    w = wcs.WCS(hdul[0].header)
    return w

"""
NAME:       getWCS
RETURNS:    New .fits file with WCS included
PARAMETERS: API key for user login to astrometry.net (key)
            Name of the file to be submitted (file)
PURPOSE:    Add WCS header to .fits file without one using astrometry.net
"""
def getWCS(key, file):
    astrometry = Client()
    astrometry.login(key)
    response = astrometry.upload(fn=file)
    jobid = astrometry.myjobs()[0]
    #print(subid)
    joburl = 'http://nova.astrometry.net/new_fits_file/' + str(jobid)
    #print(joburl)
    results = requests.get(joburl, allow_redirects=True)
    #print(results)
    filename = "CORRECTED_" + file
    file = open(filename, 'wb')
    file.write(results.content)
    file.close()
    check = True
    while check:
        try:
            hdul = fits.open(filename)
            check = False
        except OSError:
            check = True



"""
NAME:       letsGo
RETURNS:    A Photometry object with the data gained during 
            this process OR 0 if the magnitude cannot be
            computed
PARAMETERS: The right ascension of the target star in
            Decimal Degrees (targetStarRA)
            The declination of the target star in
            Decimal Degrees (targetStarDec)
            The apikey for use with astrometry.net (key)
            The main .fits file with the image data name (mainFile)
            The dark frame .fits file name (darkFrame)
                NOTE: Not used when calibrationFlag = 0
            The bias frame .fits file name (biasFrame)
                NOTE: Not used when calibrationFlag = 0
            The flat field .fits file name (flatField)
                NOTE: Not used when calibrationFlag = 0
            calibrationFlag:
                Set 0 to skip calibration, 1 to calibrate
            calibrationOutputFlag:
                Set 0 to not create an output file with the
                calibrated image, set to 1 to output calibrated
                image data.
            calibrationFile:
                Set to name of file for calibrated output, default is
                "output.fits".
            readFlag
                Set 0 to find stars, set 1 to load stars from a
                file
            magnitudeFlag:
                Set 0 to skip magnitude calculation, 1 to 
                calculate magnitude
            fwhmFlag:
                Set 0 to use a calculated radius based on
                the target star gradient, set to 1 to use 
                three times the full width half mass of the 
                target star as the radius
            printFlag:
                Set 0 to not print stars to a file, and set
                1 to print stars to a file
            readInReferenceFilename:
                The name of the file to read in from, if not
                set it will be stars.csv
            consoleFlag:
                Set to 1 to output answer to console
PURPOSE:    This method combines the methods in this file to perform 
            full differential photometry on one image file. 
"""
def letsGo(targetStarRA, targetStarDec, key,
           mainFile, darkFrame, flatField, biasFrame=0,
           calibrationFlag=1, calibrationOutputFlag=0, calibrationFile="output.fits",
           readFlag=0, magnitudeFlag=1, fwhmFlag=1,
           printFlag=1, readInReferenceFilename="stars.csv", consoleFlag=0):
    if calibrationFlag == 1:
        # Calibrate the image
        if biasFrame == 0:
            # No bias
            hdul = calibrate(mainFile, darkFrame, flatField)
        else:
            # Bias
            hdul = calibrate(mainFile, darkFrame, flatField, bias=biasFrame)
        # Print calibrated image to output
        # Create new output directory if none exists
        if not os.path.exists("Output"):
            os.mkdir("Output")
        if calibrationOutputFlag == 1:
            calibrationFile = "Output/" + calibrationFile
            hdul.writeto(calibrationFile, overwrite=True)
    else:
        # Use the raw image
        hdul = fits.open(mainFile)
    if magnitudeFlag == 0:
        # Leave the function if not calculating magnitude
        return 0

    # Check if the file has WCS data
    try:
        hdul[0].header['CRVAL1']
    except KeyError:
        getWCS(key, mainFile)
        mainFile = "CORRECTED_" + mainFile
        if calibrationFlag == 1:
            # Calibrate the image
            if biasFrame == 0:
                # No bias
                hdul = calibrate(mainFile, darkFrame, flatField)
            else:
                # Bias
                hdul = calibrate(mainFile, darkFrame, flatField, bias=biasFrame)
            # Print calibrated image to output
            # Create new output directory if none exists
            if not os.path.exists("Output"):
                os.mkdir("Output")
            if calibrationOutputFlag == 1:
                calibrationFile = "Output/" + calibrationFile
                hdul.writeto(calibrationFile, overwrite=True)
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

    # Center the star
    #pixRA, pixDec = findCenter(int(Y), int(X), hdul[0].data)

    if fwhmFlag == 1:
        try:
            radius = int(3 * hdul[0].header['FWHM'])
        except KeyError:
            radius = findRadius(Y, X, hdul[0].data)
    else:
        # Set the radius to the distance from the center
        # of the star to the farthest edge of the star
        radius = findRadius(Y, X, hdul[0].data)

    # Find the photon counts per pixel of blank sky
    blankSkyPhotons = findBlank(hdul[0].data, radius)

    # Find the photon counts in the target star
    targetStarPhotons = starCount(Y, X, hdul[0].data, radius, blankSkyPhotons)

    # Find reference stars
    if readFlag == 1:
        # Read in reference stars from file
        stars = readFromFile(readInReferenceFilename, radius, hdul[0].data, blankSkyPhotons, w)
    else:
        # Finding new stars automatically
        stars = findOtherStars(Y, X, hdul[0].data, radius, blankSkyPhotons, w)

    # Calculate magnitudes, average, and error
    ave, std, stars = calculateMagnitudeAndError(targetStarPhotons, stars)

    # Remove outliers
    if std >= (ave/20):
        print("I calculate there may be some outliers in the data. Review this list "
              + "below for outliers: ")
        for i in range(len(stars)):
            print(str(i+1) + ") " + str(stars[i].targetMagnitude))
        print("The calculated average magnitude is " + str(ave) +" and the calculated error is "
              + str(std) + ".")
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
        ave, std, stars = calculateMagnitudeAndError(targetStarPhotons, stars)

    # Console output
    if consoleFlag == 1:
        print("The magnitude of the star is ", ave )
        print("The error of this calculation is ", std)

    # Printing reference stars to files
    if printFlag == 1:
        printReferenceToFile(stars)

    # Create and return the results of the photometry
    ans = Photometry(mainFile, hdul[0].header['JD'], ave, std)
    return ans

"""
NAME:       runFiles
RETURNS:    Array of result objects containing results for photometry on every file
            in the directory
PARAMETERS: The right ascension of the target star in
            Decimal Degrees (targetStarRA)
            The declination of the target star in
            Decimal Degrees (targetStarDec)
            The apikey for use with astrometry.net (key)
            Directory containing files to process (dirName)
            The main .fits file with the image data name (mainFile)
            The dark frame .fits file name (darkFrame)
                NOTE: Not used when calibrationFlag = 0
            The bias frame .fits file name (biasFrame)
                NOTE: Not used when calibrationFlag = 0
            The flat field .fits file name (flatField)
                NOTE: Not used when calibrationFlag = 0
            calibrationFlag:
                Set 0 to skip calibration, 1 to calibrate
            calibrationOutputFlag:
                Set 0 to not create an output file with the
                calibrated image, set to 1 to output calibrated
                image data to "output.fits"
            calibrationFile:
                Set to name of file for calibrated output, default is
                "output.fits".
            readFlag
                Set 0 to find stars, set 1 to load stars from a
                file
            magnitudeFlag:
                Set 0 to skip magnitude calculation, 1 to 
                calculate magnitude
            fwhmFlag:
                Set 0 to use a calculated radius based on
                the target star gradient, set to 1 to use 
                three times the full width half mass of the 
                target star as the radius
            printFlag:
                Set 0 to not print stars to a file, and set
                1 to print stars to a file
            consoleFlag:
                Set to 1 to output answer to console
PURPOSE:    Calculates the average magnitude of the target star given the 
            reference stars, and the error in that measurement by calculating
            standard deviation of the calculated magnitudes.
"""
def runFiles(targetStarRA, targetStarDec, key,
            dirName, dark, flat, bias=0,
            calibrationFlag=1, calibrationOutputFlag=0, calibrationOutputFile="output.fits",
            readFlag=0, magnitudeFlag=1, fwhmFlag=1,
            printFlag=1, consoleFlag=0):
    path = dirName
    # Make full file names for the directory
    darkFull = dirName + "/" + dark
    flatFull = dirName + "/" + flat
    if bias != 0:
        biasFull = dirName + "/" + bias
    else:
        biasFull = 0
    results = []
    readFile = "stars.csv"
    # Check for a read in file
    if readFlag == 1:
        for filename in glob.glob(os.path.join(path, '*.csv')):
            with open(os.path.join(os.getcwd(), filename), 'r') as f:
                readFile = filename

    # Check and run all .fits files
    for filename in glob.glob(os.path.join(path, '*.fits')):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            # Check if file is data
            if filename != darkFull and filename != flatFull and filename != biasFull:
                x = letsGo(targetStarRA, targetStarDec, key, filename, darkFull, flatFull, biasFrame=biasFull,
                calibrationFlag=calibrationFlag, calibrationOutputFlag=calibrationOutputFlag,
                calibrationFile=calibrationOutputFile,
                readFlag=readFlag, magnitudeFlag=magnitudeFlag, fwhmFlag=fwhmFlag,
                printFlag=printFlag, readInReferenceFilename=readFile, consoleFlag=consoleFlag)
                if x != 0:
                    results.append(x)
    for filename in glob.glob(os.path.join(path, '*.fts')):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            if filename != darkFull and filename != flatFull and filename != biasFull:
                x = letsGo(targetStarRA, targetStarDec, key, filename, darkFull, flatFull, biasFrame=biasFull,
                           calibrationFlag=calibrationFlag, calibrationOutputFlag=calibrationOutputFlag,
                           calibrationFile=calibrationOutputFile,
                           readFlag=readFlag, magnitudeFlag=magnitudeFlag, fwhmFlag=fwhmFlag,
                           printFlag=printFlag, readInReferenceFilename=readFile, consoleFlag=consoleFlag)
                if x != 0:
                    results.append(x)
    return results