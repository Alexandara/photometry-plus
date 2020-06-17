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
 - astrometry.net
 
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
All of the functions in the program have documentation above them describing what they do.  The functions can be used separately or all at once. To do photometry on a single image, call letsGo. The parameters for letsGo are described above the function, but at its simplest it can be run with only the location of the target star in decimal degrees, the image file, and the calibration files. To run multiple files in a folder, use runFiles. runFiles takes in most of the same parameters as letsGo with the addition of dirName, which is the relative file path to the folder you wish to comb. If set to read in from a file the default stars.csv can be used or the program can find a .csv file in the folder it is looking in and use that. For more detailed instructions email alexisrenee1@gmail.com with questions.
                
# Contribute
Feel free to download this project and change it to suit your needs!

# Credit 
This project was inspired by using the Great Basin Observatory telescope in the Great Basin National Park. Great thanks to the astropy library. This project uses a modified version of the astrometry.net client.py file.

# License 
GNU @ Alexis Tudor








