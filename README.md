# photometry-plus
Photometry+ is an autonomous photometry program designed for  astronomers to perform fast photometry with as much control  as they desire. This program is specifically tailored to the Great Basin Observatory telescope data, but can be adjusted to be used for other .fits files.

# Motivation
This project was created to speed up the lengthy process that is photometry. Doing photometry by hand is slow and inefficient, and it is not easy to do for beginners. With this program it becomes much simpler to do the actual process of photometry as well as, for beginners, to look inside of the program and understand what is happening.

# Build Status
Build successful

# Tech used
 - astropy
 - photutils
 - astroquery 
 - SIMBAD
 - numpy
 - matplotlib
 
# Features
 - Calibration of .fits data and output of calibrated .fits file
 - Autonomous detection of reference stars
 - Complete pipeline for photometry, from calibration to magnitude output
 - Autonomous magnitude and identifier for reference stars from the SIMBAD database
 
# Code Output Example
The output is a .csv file containing calculated magnitudes per day: 
File Name,JD,Magnitude,Error, 
DW Cnc V-20200215at080749_-25-1X1-300-V.fts,2458894.8388310187,15.65760283818426,0.37986752525236084, 

# Installation
This code was run using PyCharm which can be downloaded here: https://www.jetbrains.com/pycharm/

# API Reference
SIMBAD API: https://astroquery.readthedocs.io/en/latest/simbad/simbad.html

# Example Function Calls
letsGo(119.721, 16.279214, "DW Cnc V-20200215at080749_-25-1X1-300-V.fts", dark, bias, flat, fwhmFlag=1, readFlag=1, printFlag=1, calibrationOutputFlag=1)

letsGo(119.721, 16.279214, "DW Cnc V-20200215at080749_-25-1X1-300-V.fts", dark, bias, flat)

# How to start:
Just download photometry.py and call the function 'letsGo' from the main method with the following parameters (note all flags are optional parameters): 

NAME:       letsGo
RETURNS:    An Answer object with the data gained during 
            this process
PARAMETERS: The right ascension of the target star inDecimal Degrees (targetStarRA)
            The declination of the target star in Decimal Degrees (targetStarDec)
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
                
# Contribute
Feel free to download this project and change it to suit your needs!

# Credit 
This project was inspired by using the Great Basin Observatory telescope in the Great Basin National Park. Great thanks to the astropy library.

# License 
GNU @ Alexis Tudor








