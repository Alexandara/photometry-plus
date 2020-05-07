from __future__ import division
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from scipy.interpolate import make_interp_spline, BSpline

import matplotlib.pyplot as plt
from photutils import CircularAperture
from photutils import Background2D, MedianBackground
from astropy.io import fits
from astropy import wcs
from photutils import aperture_photometry
import astropy.coordinates as coord
import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from matplotlib.colors import LogNorm
from skimage import data, img_as_float
from skimage import exposure


class Star:
    def __init__(self, id, r, d, rad, c, m):
        self.id = id
        self.ra = r #Pixel location
        self.dec = d #Pixel location
        self.radius = rad
        self.counts = c
        self.magnitude = m

class Answer:
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
            the bias frame (bias)
            the flat frame (flat)
PURPOSE:    to calibrate raw .fits files into a form that can 
            be used to calculate accurate magnitude data. 
            This currently only uses dark and flat but can 
            easily be modified to use bias as well
"""
def calibrate(filename, dark, bias, flat):
    # Open the files
    hdul = fits.open(filename)  # hdul is the computer version of
                                # the file data
    #the following hdul files are computer readable versions of
    #the dark, bias, and flat .fit files
    hdulDark = fits.open(dark)
    hdulBias = fits.open(bias)
    hdulFlat = fits.open(flat)

    # Extract data from the files
    data = hdul[0].data
    dataDark = hdulDark[0].data
    #dataBias = hdulBias[0].data # Use this line for a bias frame.
                                 # Leave it out if the bias is included
                                 # in the dark frame
    dataFlat = hdulFlat[0].data
    x, y = data.shape

    # Calibrate the files, and then reassign the calibrated data to the
    # original data
    data = data - dataDark #Subtract the dark frame from the data
    #data = data - dataBias #Use this line to subtract bias separately
    dataFlatNorm = dataFlat / np.mean(dataFlat) #Normalize the flat frame information
    data = data // dataFlatNorm #Divide out the normalized flat frame data
    hdul[0].data = data #This sets the calibrated data to be a part of the
                        #hdul object
    return hdul

"""
NAME:       findCenter
RETURNS:    RA and DEC of the magnitude center of the star
PARAMETERS: Center of star pixel Y location (ra)
            Center of star pixel X location (dec)
            Sky data (data)
PURPOSE:    To take the recorded RA and DEC of a star and 
            center it for the image itself at the point of 
            highest photon counts
"""
def findCenter(ra, dec, data):
    high = 0
    highRA = ra
    highDec = dec
    #Search an area 100 pixels around the reported center
    for i in range(100):
        for j in range(100):
            #If the pixel being looked at has the highest
            #photon count so far, high should be set to that
            #amount and the RA and DEC should be set to that
            #location
            if data[dec + 50 - j][ra + 50 - i] > high:
                high = data[dec + 50 - j][ra + 50 - i]
                highRA = ra + 50 - i
                highDec = dec + 50 - j
    #Return the true center of the star
    return highRA, highDec

"""
NAME:       findBlankNoRad
RETURNS:    Average counts in blank sky
PARAMETERS: Center of star pixel Y location (ra)
            Center of star pixel X location (dec)
            Array of image data for the sky (data1)
PURPOSE:    In the case the FWHM cannot be used as the radius
            of the star, this function can be used to find the 
            counts in a blank portion of the sky without the 
            radius of the main star.
"""
def findBlankNoRad(ra, dec, data1):
    #FindBlankNoRad creates a less correct blank sky count that
    # does not use the radius of the star in any way
    #This allows the radius to be calculated and used to
    # calculate a better blank sky count
    flag = 0
    low = data1[dec][ra]
    y = ra
    x = dec
    currY = ra
    currX = dec
    nextY = ra
    nextX = dec
    # Compute Gradient Descent to find an empty pixel spot
    # This traverses all pixels towards a local minimum, which
    # is then used to find the blank pixel average
    while flag == 0:
        # This checks in every direction for a lower spot,
        # and then steps in that direction
        for i in range(3):
            for j in range(3):
                if data1[currX + 1 - i][currY + 1 - j] < low:
                    low = data1[currX + 1 - i][currY + 1 - j]
                    nextX = currX + 1 - i
                    nextY = currY + 1 - j
        # If there was no movement in the previous step,
        # this sets the end flag to one, causing the while loop
        # to end
        if currX == nextX and currY == nextY:
            flag = 1
        currX = nextX
        currY = nextY
    x = currX
    y = currY
    avBlank = 0
    r = 100 # This is an assumed radius, this can be modified as needed
    #This adds together all pixels in an area 4 times the radius
    # of the target star
    if x <= dec and y <= ra:
        for i in range(r):
            for j in range(r):
                avBlank = avBlank + data1[x - i][y - j]
    if x >= dec and y <= ra:
        for i in range(r):
            for j in range(r):
                avBlank = avBlank + data1[x + i][y - j]
    if x <= dec and y >= ra:
        for i in range(r):
            for j in range(r):
                avBlank = avBlank + data1[x - i][y + j]
    if x >= dec and y >= ra:
        for i in range(r):
            for j in range(r):
                avBlank = avBlank + data1[x + i][y + j]
    #This averages the blank pixels over the area and then
    # returns it
    avBlank = avBlank / (r * r)
    return avBlank

"""
NAME:       findRadius
RETURNS:    Detected max radius
PARAMETERS: Center of star pixel Y location (ra)
            Center of star pixel X location (dec)
            Image array (data)
PURPOSE:    This finds the radius of a target star by detecting
            the transition to blank sky. This is only for use in 
            cases where using the FWHM is not possible.
"""
def findRadius(ra, dec, data):

    x, y = data.shape
    prev = 10000
    curr = 10000
    currY = ra
    currX = dec
    #Blank is the number of counts per pixel in empty sky
    blank = findBlankNoRad(currY, currX, data)
    r1 = 0
    r2 = 0
    r3 = 0
    r4 = 0
    while data[currX][currY] > blank:
        currY = currY + 1
        if currY >= y:
            currY = ra
            break
    r1 = currY - ra
    currY = ra
    while data[currX][currY] > blank:
        currY = currY - 1
        if currY <= 1:
            currY = ra
            break
    r2 = -1 * (currY - ra)
    currY = ra
    while data[currX][currY] > blank:
        currX = currX - 1
        if currX <= 1:
            currX = dec * -1
            break
    r3 = -1 * (currX - dec)
    currX = dec
    while data[currX][currY] > blank:
        currX = currX + 1
        if currX >= x:
            currX = dec
            break
    r4 = currX - dec

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
PARAMETERS: Center of star pixel Y location (ra)
            Center of star pixel X location (dec)
            Image array (data1)
            Radius of primary target star (rad)
PURPOSE:    Using the radius of the primary target star, this 
            method averages out the blank counts over an area 4 times 
            the radius of the primary target star and returns it.
NOTE:       See findBlankNoRad for detailed comments on the procedure.
"""
def findBlank(ra, dec, data1, rad):
    # Data Assignment
    flag = 0
    low = data1[dec][ra]
    y = ra
    x = dec
    currY = ra
    currX = dec
    nextY = ra
    nextX = dec
    options = []
    # Compute Gradient Descent to find an empty pixel spot
    while flag == 0:
        for i in range(3):
            for j in range(3):
                if data1[currX + 1 - i][currY + 1 - j] < low:
                    low = data1[currX + 1 - i][currY + 1 - j]
                    nextX = currX + 1 - i
                    nextY = currY + 1 - j
        if currX == nextX and currY == nextY:
            flag = 1
        currX = nextX
        currY = nextY
    x = currX
    y = currY
    avBlank = 0
    r = 4 * rad
    if x <= dec and y <= ra:
        for i in range(r):
            for j in range(r):
                options.append(data1[x - i][y - j])
                avBlank = avBlank + data1[x - i][y - j]
    if x >= dec and y <= ra:
        for i in range(r):
            for j in range(r):
                options.append(data1[x - i][y - j])
                avBlank = avBlank + data1[x + i][y - j]
    if x <= dec and y >= ra:
        for i in range(r):
            for j in range(r):
                options.append(data1[x - i][y - j])
                avBlank = avBlank + data1[x - i][y + j]
    if x >= dec and y >= ra:
        for i in range(r):
            for j in range(r):
                options.append(data1[x - i][y - j])
                avBlank = avBlank + data1[x + i][y + j]
    avBlank = avBlank / (16 * r * r)
    return options[int(len(options)/2)] # SORT THIS ARRAY

"""
NAME:       starCount
RETURNS:    The amount of photons in a star's image
PARAMETERS: Center of star pixel Y location (ra)
            Center of star pixel X location (dec)
            Image data (data)
            Radius of the star (r)
            Average count per blank sky pixel (blank)
PURPOSE:    Taking in information about any star, this method counts
            the photons in a circular area around the center of the star.
"""
def starCount(ra, dec, data, r, blank):
    A = math.pi * r * r
    pix = 0
    sum = 0
    for i in range(r * 2):
        for j in range(r * 2):
            # Checks if the pixel being checked is within the radius of the
            # star. If so, that pixel's count is added to the total.
            if (math.sqrt((ra - (ra + r - i)) * (ra - (ra + r - i))
                          + (dec - (dec + r - j)) * (dec - (dec
                            + r - j)))) <= r:
                pix = pix + 1
                sum = sum + data[int(dec + r - j)][int(ra + r - i)]
    #print(r, sum, pix, blank)
    tot = sum - (pix * blank)
    return tot

"""
NAME:       starMake
RETURNS:    An array of stars including any star found in this direction
PARAMETERS: Starting DEC for searching (X)
            Starting RA for searching (Y)
            Image data (data)
            Radius to use (rad)
            Counts in blank sky (blank)
            Coordinate to pixel parameter (w)
            Array containing previously searched for stars (stars)
            Length of the image array (length)
            Width of the image array (width)
            Integer indicating what direction to search in (dir)
PURPOSE:    This searches in a certain direction to find a star that
            has a V magnitude to add to the star array. If it doesn't 
            find anything, nothing is added to the array and it is 
            returned as is. The star is only added if it is properly 
            centered and no other stars are nearby. 
"""
def starMake(X, Y, data, rad, blank, w, stars, length, width, dir):
    # If the flag is changed to 1, a star with a magnitude is found.
    flag = 0
    # While no star has been found and while we are not searching
    # beyond the edge of the image array
    while flag == 0 and (X * X) < ((length - 100) * (length - 100)) \
            and (Y * Y) < ((width - 100) * (width - 100)):
        #Check if this pixel has more counts than a blank pixel
        if data[X][Y] > (blank+10):
            #Find the center of the possible star
            tRA = Y
            tDEC = X
            #tRA, tDEC = findCenter(Y, X, data)
            #Set the radius of the star to the radius used
            tR = rad
            #Get the counts in the possible star
            tC = starCount(tRA, tDEC, data, tR, blank)
            #Get the actual coordinates of the star
            # lon, lat = w.all_pix2world(tRA, tDEC, 0) #Line used for manual magnitude
            skyC = SkyCoord.from_pixel(tRA, tDEC, w)
            #Get the magnitude of the star and check it
            customSimbad = Simbad()
            customSimbad.add_votable_fields('flux(V)')
            result = customSimbad.query_region(skyC) #This line creates errors
            if not(result is None):
                Vmag = result['FLUX_V']
                if len(Vmag) == 1:
                    stars.append(Star(result['MAIN_ID'][0], tRA, tDEC, tR, tC, Vmag[0]))
                    flag = 1
            """
            # This code segment is used for manual input of magnitudes
            print(dir, ": What is the magnitude of the star at ", lon, lat,
                  ", if that star has no magnitude, please enter 100.")
            mag = input()
            # Check if the star has a magnitude, set flag to 1 if yes and add
            # to stars array
            if float(mag) != float(100):
                stars.append(Star(tRA, tDEC, tR, tC, mag))
                flag = 1
            """
            if dir == 1:
                X = X + 100
                Y = Y + 100
            elif dir == 2:
                X = X + 100
            elif dir == 3:
                X = X + 100
                Y = Y - 100
            elif dir == 4:
                Y = Y - 100
            elif dir == 5:
                X = X - 100
                Y = Y - 100
            elif dir == 6:
                X = X - 100
            elif dir == 7:
                X = X - 100
                Y = Y + 100
            elif dir == 8:
                Y = Y + 100
        if dir == 1:
            X = X + 1
            Y = Y + 1
        elif dir == 2:
            X = X + 1
        elif dir == 3:
            X = X + 1
            Y = Y - 1
        elif dir == 4:
            Y = Y - 1
        elif dir == 5:
            X = X - 1
            Y = Y - 1
        elif dir == 6:
            X = X - 1
        elif dir == 7:
            X = X - 1
            Y = Y + 1
        elif dir == 8:
            Y = Y + 1
    return stars

"""
NAME:       findOtherStars
RETURNS:    Array of other star objects
PARAMETERS: Center of target star RA (ra)
            Center of target star DEC (dec)
            Image data (data)
            Radius to examine (rad)
            Counts in blank sky (blank)
            Parameters for conversion from world coordinates to
            pixel location (w)
PURPOSE:    This method locates stars nearby the target stars
            and then calls a helper method to construct a star
            object for the located stars.
"""
def findOtherStars(ra, dec, data, rad, blank, w):
    length = len(data)
    width = len(data[0])
    #Create an empty stars array to fill
    stars = []

    # TOP RIGHT
    X = dec + rad + 10
    Y = ra + rad + 10
    stars = starMake(X, Y, data, rad, blank, w, stars, length, width, 1)

    # RIGHT
    X = dec + rad + 10
    Y = ra
    stars = starMake(X, Y, data, rad, blank, w, stars, length, width, 2)

    # BOTTOM RIGHT
    X = dec + rad + 10
    Y = ra - rad - 10
    stars = starMake(X, Y, data, rad, blank, w, stars, length, width, 3)

    # BOTTOM
    X = dec
    Y = ra - rad - 10
    stars = starMake(X, Y, data, rad, blank, w, stars, length, width, 4)

    # BOTTOM LEFT
    X = dec - rad - 10
    Y = ra - rad - 10
    stars = starMake(X, Y, data, rad, blank, w, stars, length, width, 5)

    # LEFT
    X = dec - rad - 10
    Y = ra
    stars = starMake(X, Y, data, rad, blank, w, stars, length, width, 6)

    # TOP LEFT
    X = dec - rad - 10
    Y = ra + rad + 10
    stars = starMake(X, Y, data, rad, blank, w, stars, length, width, 7)

    # TOP
    X = dec
    Y = ra + rad + 10
    stars = starMake(X, Y, data, rad, blank, w, stars, length, width, 8)

    return stars

"""
NAME:       printToFile
RETURNS:    Nothing
PARAMETERS: Array of star objects (stars)
PURPOSE:    This method prints out a file containing
            the stars used for calculating magnitude
"""
def printToFile(stars, w):
    file = open("stars.csv", "w")
    file.write("Name,Right Ascension (Decimal Degrees),Declination (Decimal Degrees),Radius (pixels),Photons,Magnitude,\n")
    for i in range(len(stars)):
        realRA, realDec = w.all_pix2world(stars[i].ra, stars[i].dec, 0)
        file.write(str(stars[i].id) + "," + str(realRA) + "," + str(realDec)
        + "," + str(stars[i].radius) + "," + str(stars[i].counts) +
                   "," + str(stars[i].magnitude) + ", \n")
    file.close()

"""
NAME:       readFromFile
RETURNS:    An array of star objects
PARAMETERS: The name of the file to read from (fileName)
PURPOSE:    This method reads in data on stars from a
            file to use to calculate the magnitude.
"""
def readFromFile(fileName, radius, data, blank, w):
    file = open(fileName, "r")
    stars = []
    for line in file:
        array = line.split(",")
        if not(array[0] == 'Name') and not(array[0] == '\n') and not(array[0] == ''):
            pixRA, pixDec = w.all_world2pix(float(array[1]), float(array[2]), 0)
            aperture = CircularAperture((pixRA, pixDec), r=radius)
            starTable = aperture_photometry(data, aperture)
            starPhotons = starTable['aperture_sum'][0]
            starPhotons = starPhotons - ((math.pi * radius * radius) * blank)
            stars.append(Star(array[0], pixRA, pixDec, radius, starPhotons, float(array[5])))
    return stars

"""
NAME:       output
RETURNS:    nothing
PARAMETERS: An array of Answer objects (info)
PURPOSE:    Output a file with answers including the file name,
            the JD (Julian Date) of the observation, the magnitude 
            reported, and the error of that reported magnitude
"""
def output(info):
    file = open("output.csv", "w")
    file.write("File Name,JD,Magnitude,Error, \n")
    for i in range(len(info)):
        file.write(info[i].fileName + "," + str(info[i].JD) + "," + str(info[i].magnitude) + ","
                   + str(info[i].error) + ", \n")
    file.close()

def plot(filename):
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
    xnew = np.linspace(min(jd), max(jd), 300)
    spl = make_interp_spline(jd, mag, k=3)  # type: BSpline
    power_smooth = spl(xnew)
    plt.title('Light Curve of DW Cnc from February 15th to March 29th')
    plt.errorbar(jd, mag, err, fmt='ko')
    #Smooth line:
    #plt.plot(xnew, power_smooth)
    plt.xlabel('Julian Day')
    plt.ylabel('Magnitude')
    plt.gca().invert_yaxis()
    plt.show()  # to print to screen
    plt.savefig("chart.pdf")  # to save to file

"""
NAME:       letsGo
RETURNS:    An Answer object with the data gained during 
            this process
PARAMETERS: The right ascension of the target star in
            Decimal Degrees (targetStarRA)
            The declination of the target star in
            Decimal Degrees (targetStarDec)
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
PURPOSE:    This is where it all begins.
"""
def letsGo(targetStarRA, targetStarDec,
           mainFile, darkFrame, biasFrame, flatField,
           calibrationFlag=1, calibrationOutputFlag=0,
           readFlag=0, magnitudeFlag=1, fwhmFlag=1,
           printFlag=1):
    if calibrationFlag == 1:
        hdul = calibrate(mainFile, darkFrame, biasFrame, flatField)
        if calibrationOutputFlag == 1:
            hdul.writeto("output.fits", overwrite=True)
    else:
        hdul = fits.open(mainFile)
    if magnitudeFlag == 0:
        return 0
    # w is the reference of world coordinates for this image
    w = wcs.WCS(hdul[0].header)
    #print(w.all_pix2world(2046,2046,0))
    # Convert COORDINATE system data into pixel locations
    X, Y = w.all_world2pix(targetStarRA, targetStarDec, 0)
    # Centers the target star location for this specific image
    pixRA = int(Y)
    pixDec = int(X)
    #pixRA, pixDec = findCenter(int(Y), int(X), hdul[0].data)
    if fwhmFlag == 1:
        radius = int(3 * hdul[0].header['FWHM'])
    else:
        # Set the radius to the distance from the center
        # of the star to the farthest edge of the star
        radius = findRadius(pixRA, pixDec, hdul[0].data)
    # Find the photon counts per pixel of blank sky
    #blankSkyPhotons = findBlank(pixRA, pixDec, hdul[0].data, radius)
    background = Background2D(hdul[0].data, ((radius*4),(radius*4)))
    blankSkyPhotons = background.background_median
    # Find the photon counts in the target star
    #targetStarPhotons = starCount(pixRA, pixDec, hdul[0].data, radius, blankSkyPhotons)
    targetAperture = CircularAperture((pixDec,pixRA), r=radius)
    targetStarTable = aperture_photometry(hdul[0].data, targetAperture)
    targetStarPhotons = targetStarTable['aperture_sum'][0]
    targetStarPhotons = targetStarPhotons - ((math.pi*radius*radius)*blankSkyPhotons)

    # Find an array with up to 8 other stars with magnitudes
    if readFlag == 0:
        stars = findOtherStars(pixRA, pixDec, hdul[0].data, radius, blankSkyPhotons, w)
    elif readFlag == 1:
        stars = readFromFile("stars.csv", radius, hdul[0].data, blankSkyPhotons, w)
    # Find the magnitude relative to each other star, and then averages them
    ave = 0
    magnitudes = []
    std = 0
    #Calculate average magnitude
    for i in range(len(stars)):
        #print(targetStarPhotons, stars[i].counts, stars[i].magnitude)
        #print(targetStarPhotons, stars[i].counts, stars[i].magnitude, blankSkyPhotons)
        if(targetStarPhotons > 0) and (stars[i].counts > 0):
            magnitudes.append((-2.5) * math.log10(targetStarPhotons / stars[i].counts)
                                 + float(stars[i].magnitude))
            ave = ave + ((-2.5) * math.log10(targetStarPhotons / stars[i].counts)
                                 + float(stars[i].magnitude))
    if len(magnitudes) == 0:
        return 0;
    ave = ave / len(magnitudes)
    #Calculate std of magnitudes
    for i in range(len(magnitudes)):
        std = std + ((magnitudes[i] - ave) * (magnitudes[i] - ave))
    std = std / len(magnitudes)
    std = math.sqrt(std)
    #Remove outliers
    if std >= (ave/20):
        print("I calculate there may be some outliers in the data. Review this list "
              + "below for outliers: ")
        for i in range(len(magnitudes)):
            print(str(i+1) + ") " + str(magnitudes[i]))
        print("The calculated average magnitude is " + str(ave) +" and the calculated error is "
              + str(std) + ".")
        print("\nPlease input the number next to the magnitudes you want to remove from the calculations: "
              + "(enter 100 if there are no more outliers to remove) ")
        x = input()
        x = int(x)
        toRemove = []
        toRemoveStars = []
        while x != 100:
            toRemove.append(magnitudes[x-1])
            toRemoveStars.append(stars[x-1])
            x = input()
            x = int(x)
        for i in range(len(toRemove)):
            magnitudes.remove(toRemove[i])
            stars.remove(toRemoveStars[i])
    ave = 0
    std = 0
    for i in range(len(magnitudes)):
        #print(targetStarPhotons, stars[i].counts, stars[i].magnitude)
        ave = ave + magnitudes[i]
    ave = ave / len(magnitudes)
    #Calculate std of magnitudes
    for i in range(len(magnitudes)):
        std = std + ((magnitudes[i] - ave) * (magnitudes[i] - ave))
    std = std / len(magnitudes)
    std = math.sqrt(std)
    #print("The magnitude of the star is ", ave )
    #print("The error of this calculation is ", std)
    if printFlag == 1:
        printToFile(stars, w)
    ans = Answer(mainFile, hdul[0].header['JD'], ave, std)
    return ans

