# Photometry+ by Alexis Tudor 
Photometry+ is an autonomous photometry program designed for  astronomers to perform fast photometry with as much control  as they desire. This program is specifically tailored to the Great Basin Observatory telescope data, but can be adjusted to be used for other .fits files.

# Motivation
This project was created to speed up the lengthy process that is differential photometry. Doing differential photometry by hand is slow and inefficient, and it is not easy to do for beginners. With this program it becomes much simpler to do the actual process of photometry and to look inside of the program and understand what is happening.

# Build Status
Build successful

# Full Guide and Tutorials
https://docs.google.com/document/d/156hhJvwQ5JsuQsuCONkSsc4z2qmjbsZ71uLW4XOGpw8/edit?usp=sharing

# Tech used
 - astropy [![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
 - photutils
 - astroquery 
 - SIMBAD
 - numpy
 - matplotlib
 - Astrometry.net
 - VizieR
 - DAOFind
 - PyQT5 
 
# Features
 - Calibration of .fits data and output of calibrated .fits file
 - Autonomous detection of reference stars
 - Complete pipeline for photometry, from calibration to magnitude output
 - Autonomous magnitude and identifier for reference stars from the SIMBAD database
 - Use of astrometry.net to identify WCS for images without WCS attached
 - Plot a light curve based on photometry data
 - Expansive error calculation and correction options
 - User Interface that can be used to teach differential photometry
 
# Code Output Example
**Results file output:**   
File Name,JD,Magnitude,Error,    
DW Cnc V-20200215at080749_-25-1X1-300-V.fts,2458894.8388310187,15.65760283818426,0.37986752525236084,    
X,Reduced Chi Square: 1.20,

# Installation
This code was developed using PyCharm which can be downloaded here: https://www.jetbrains.com/pycharm/

# API Reference
SIMBAD API: https://astroquery.readthedocs.io/en/latest/simbad/simbad.html   
Astrometry.net API: http://nova.astrometry.net/api_help   
VizieR API: https://astroquery.readthedocs.io/en/latest/vizier/vizier.html   

# Example Function Calls
runPhotometry(95.685417, -0.34563889, "MainFile.fts")

runFiles(52.800196, 43.90441, "MainFolder", "DarkFolder", "BiasFolder", "FlatFolder")

# How to use the user interface:
Run the following code to start the user interface:  

import photometryuserinterface  
from PyQt5 import QtWidgets 
import sys

if __name__ == '__main__':  
    app = QtWidgets.QApplication(sys.argv) 
    controller =photometryuserinterface.Controller()  
    controller.show_home()  
    sys.exit(app.exec_())  

# How to start:
All of the functions in the program have documentation above them describing what they do.  The functions can be used separately or all at once. To do photometry on a single image, call runPhotometry. The parameters for runPhotometry are described above the function, but at its simplest it can be run with only the location of the target star in decimal degrees and the image file.   
  
To run multiple files for the same target star in a folder, use runFiles. runFiles takes in most of the same parameters as runPhotometry with the addition of dirName, which is the relative file path to the folder you wish to comb. For more detailed instructions email alexisrenee1@gmail.com with questions.
                
# Contribute
Feel free to download this project and change it to suit your needs. Credit Photometry+ if you use any part of the code.

# How to Cite Photometry+
**In Publications:**

If you use Photometry+ for work/research presented in a publication (whether directly, or as a dependency to another package), we ask that you please use the following citation:
     
     Tudor, A. 2020, Photometry+, v3.0, Github at https://github.com/Alexandara/photometry-plus

If there is no place to cite the papers, please use this acknowledgement:

     This research made use of Photometry+ at https://github.com/Alexandara/photometry-plus.
     
**In Projects:**

If you are using Photometry+ as part of a code project (e.g., affiliated packages), a useful way to acknowledge your use of Photometry+ is in your README. Use this badge: [![Generic badge](https://img.shields.io/badge/powered%20by-Photometry+-blue.svg)](https://github.com/Alexandara/photometry-plus)

**In Presentations:**

When using Photometry+ in a presentation, please include the Photometry+ logo in your presentation:
![Photometry+ Logo (2)](https://user-images.githubusercontent.com/6069321/86058691-23bf2e00-ba16-11ea-8f97-5ef990d68a4c.png)

# Credit 
This project was inspired by using the Great Basin Observatory telescope in the Great Basin National Park. Great thanks to the astropy library. This project uses a modified version of the astrometry.net client.py file.

# License 
GNU @ Alexis Tudor
[![Open Source Love svg2](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)](https://github.com/ellerbrock/open-source-badges/)








