"""
Program written by Alexis Tudor at the University of Nevada, Reno
Email at alexisrenee1@gmail.com
Copyright and Licensing: GNU @ Alexis Tudor
"""
import photometry
import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk
from tkinter import messagebox
from tkinter import font as tkfont
import os
from django.template.defaultfilters import slugify
import tkinter.filedialog

#root.option_add('*tearoff', False)
#filename = filedialog.askopenfile()

class Project:
    def __init__(self, name, ra, dec, file, dark, bias, flat):
        self.projectName = name
        self.targetStarRA = ra
        self.targetStarDec = dec
        self.mainFile = file
        self.darkFile = dark
        self.biasFile = bias
        self.flatFile = flat
        self.projectSettings = photometry.Settings()

class Start(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        self.title_font = tkfont.Font(family='Helvetica', size=18, weight="bold", slant="italic")

        # the container is where we'll stack a bunch of frames
        # on top of each other, then the one we want visible
        # will be raised above the others
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        # Menu
        self.logoimg = Image.open("assets/Photometry+ Logo (gold).jpg")
        self.logoimg = self.logoimg.resize((400, 250), Image.ANTIALIAS)
        self.logo = ImageTk.PhotoImage(self.logoimg)

        menubar = tk.Menu(self)
        self.config(menu=menubar)
        file = tk.Menu(menubar)
        help_ = tk.Menu(menubar)

        menubar.add_cascade(menu=file, label='File')
        menubar.add_cascade(menu=help_, label='Help')

        file.add_command(label='New Project +', command=lambda: self.show_frame("ProjectPage"))

        self.smalllogoimg = self.logoimg.resize((10, 10), Image.ANTIALIAS)
        self.smallLogo = ImageTk.PhotoImage(self.smalllogoimg)
        file.entryconfig('New Project +', accelerator='Ctrl+N', image=self.smallLogo, compound='left')

        """
        # How to add menus in menus
        file.delete('Save')
        save = Menu(file)
        file.add_cascade(menu=save, label='Save')
        save.add_command(label='Save As', command=lambda: print('Saving As...'))
        save.add_command(label='Save All', command=lambda: print('Saving All...'))
        """

        self.frames = {}
        for F in (HomePage, ProjectPage):
            page_name = F.__name__
            frame = F(parent=container, controller=self)
            self.frames[page_name] = frame

            # put all of the pages in the same location;
            # the one on the top of the stacking order
            # will be the one that is visible.
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame("HomePage")

    def show_frame(self, page_name):
        '''Show a frame for the given page name'''
        frame = self.frames[page_name]
        frame.tkraise()

class HomePage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        # Blue = #100197

        # Get size of computer screen
        screenWidth = self.winfo_screenwidth()
        screenHeight = self.winfo_screenheight()

        # Configure screen
        self.columnconfigure(0, weight=1, minsize=screenWidth*.2)
        self.columnconfigure(1, weight=1, minsize=screenWidth*.6)
        self.rowconfigure(0, weight=1, minsize=screenHeight*.6)
        self.rowconfigure(1, weight=1, minsize=screenHeight*.2)

        # Set style

        self.style = ttk.Style()
        self.style.configure('TFrame', background='#100197')
        self.style.configure('Big.TButton', background='#100197', font=('Space Mono', 20, 'bold'))
        self.style.configure('TLabel', background='#100197')#, font=('Arial', 11))
        self.style.configure('Header.TLabel', font=('Arial', 18, 'bold'))

        # Photos
        self.spaceimg = Image.open("assets/Space.gif")
        self.spaceimg = self.spaceimg.resize((450, 350), Image.ANTIALIAS)
        self.space = ImageTk.PhotoImage(self.spaceimg)
        self.logoimg = Image.open("assets/Photometry+ Logo (gold).jpg")
        self.logoimg = self.logoimg.resize((400, 250), Image.ANTIALIAS)
        self.logo = ImageTk.PhotoImage(self.logoimg)

        # Start Screen Logo
        self.logoLabel = ttk.Label(self)
        self.logoLabel.image = self.logo
        self.logoLabel.config(image=self.logo)
        self.logoLabel.grid(row=0, column=1)

        # New Project Button
        self.newProjectButton = ttk.Button(self, text="New Project + ", style='Big.TButton', command=lambda: controller.show_frame("ProjectPage"))
        self.newProjectButton.grid(row=1, column=1, sticky=tk.N)

        # Treeview
        treeview = ttk.Treeview(self)
        treeview.grid(row=0,column=0,rowspan=2, sticky=tk.NSEW)
        treeview.insert('', '0', 'home', text='Home')
        treeview.insert('', '1', 'about', text='About')
        treeview.insert('', '2', 'settings', text='Settings')
        treeview.insert('', '3', 'myprojects', text='My Projects')

class ProjectPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        # Get size of computer screen
        screenWidth = self.winfo_screenwidth()
        screenHeight = self.winfo_screenheight()

        # Images
        self.logoimg = Image.open("assets/Photometry+ Logo (gold).jpg")
        self.logoimg = self.logoimg.resize((120, 75), Image.ANTIALIAS)
        self.logo = ImageTk.PhotoImage(self.logoimg)

        # Configure screen
        self.columnconfigure(0, weight=1, minsize=screenWidth * .4)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)
        self.columnconfigure(3, weight=1)
        self.columnconfigure(4, weight=1)

        # Logo
        self.logoLabel = ttk.Label(self)
        self.logoLabel.image = self.logo
        self.logoLabel.config(image=self.logo)
        self.logoLabel.grid(row=0, column=0, sticky=tk.NW)
        self.logoLabel.bind('<Button-1>', self.gohome)

        # Settings Block
        # Name of Project
        ttk.Label(self, text="Project Name").grid(row=0, column=1)
        self.projectName = ttk.Entry(self)
        self.projectName.grid(row=0, column=2, columnspan=3, sticky=tk.E+tk.W)

        # Coordinates
        ttk.Label(self, text="Target Star Coordinates").grid(row=1,column=1, columnspan=4)
        ttk.Label(self, text="Right Ascension").grid(row=2,column=1)
        ttk.Label(self, text="Declination").grid(row=3,column=1)
        self.decimalDegreeRA = ttk.Entry(self)
        self.decimalDegreeDec = ttk.Entry(self)
        self.degreeRA = ttk.Entry(self)
        self.minuteRA = ttk.Entry(self)
        self.secondRA = ttk.Entry(self)
        self.degreeDec = ttk.Entry(self)
        self.minuteDec = ttk.Entry(self)
        self.secondDec = ttk.Entry(self)

        self.coordinateChoice = tk.StringVar()
        self.coordinateComboBox = ttk.Combobox(self, textvariable=self.coordinateChoice, state="readonly")
        self.coordinateComboBox.grid(row=4, column=1, columnspan=4)
        self.coordinateComboBox.config(values=('Decimal Degrees', 'Degrees'))

        if photometry.getcoordinateChoiceFlag() == "DEG":
            self.degreeRA.grid(row=2, column=2, sticky=tk.E+tk.W)
            self.minuteRA.grid(row=2, column=3, sticky=tk.E+tk.W)
            self.secondRA.grid(row=2, column=4, sticky=tk.E+tk.W)
            self.degreeDec.grid(row=3, column=2, sticky=tk.E+tk.W)
            self.minuteDec.grid(row=3, column=3, sticky=tk.E+tk.W)
            self.secondDec.grid(row=3, column=4, sticky=tk.E+tk.W)
            self.coordinateComboBox.current([1])
        else:
            self.decimalDegreeRA.grid(row=2, column=2, columnspan=3, sticky=tk.E+tk.W)
            self.decimalDegreeDec.grid(row=3, column=2, columnspan=3, sticky=tk.E+tk.W)
            self.coordinateComboBox.current([0])

        # Calibration File Uploads
        if photometry.getcalibrationFlag() == 1:
            ttk.Label(self, text="Calibration File Uploads").grid(row=5, column=1, columnspan=4)

            ttk.Label(self, text="Dark Frame").grid(row=6, column=1, columnspan=2)
            self.biasButton = ttk.Button(self, text="Select File or Directory", command=self.getDark)
            self.biasButton.grid(row=6, column=3, columnspan=2)

            ttk.Label(self, text="Bias Frame").grid(row=7, column=1, columnspan=2)
            self.biasButton = ttk.Button(self, text="Select File or Directory", command=self.getBias)
            self.biasButton.grid(row=7, column=3, columnspan=2)

            ttk.Label(self, text="Flat Field").grid(row=8, column=1, columnspan=2)
            self.biasButton = ttk.Button(self, text="Select File or Directory", command=self.getFlat)
            self.biasButton.grid(row=8, column=3, columnspan=2)

        # Settings Start
        ttk.Label(self, text="Settings").grid(row=9, column=1, columnspan=4)

        # API key
        ttk.Label(self, text="Astrometry.net API Key").grid(row=10, column=1, columnspan=2)
        self.apiKeyEntry = ttk.Entry(self)
        self.apiKeyEntry.grid(row=10, column=3, columnspan=2, sticky=tk.E + tk.W)
        if not(photometry.getastrometryDotNetFlag() == 0):
            self.apiKeyEntry.insert(0, str(photometry.getastrometryDotNetFlag()))

        # calibrationFlag setting
        ttk.Label(self, text="Calibrate Image?").grid(row=11, column=1, columnspan=2)
        self.calibrationChoice = tk.StringVar()
        self.calibrationComboBox = ttk.Combobox(self, textvariable=self.calibrationChoice, values=('No', 'Yes'), state="readonly")
        self.calibrationComboBox.grid(row=11, column=3, columnspan=2, sticky=tk.E + tk.W)
        self.calibrationComboBox.current([photometry.getcalibrationFlag()])

        # useDarkFlag setting
        ttk.Label(self, text="Use Dark Frame?").grid(row=12, column=1, columnspan=2)
        self.useDarkChoice = tk.StringVar()
        self.useDarkComboBox = ttk.Combobox(self, textvariable=self.useDarkChoice, values=('No', 'Yes'), state="readonly")
        self.useDarkComboBox.grid(row=12, column=3, columnspan=2, sticky=tk.E + tk.W)
        self.useDarkComboBox.current([photometry.getuseDarkFlag()])

        # subtractBiasFromDarkFlag setting
        ttk.Label(self, text="Subtract Bias from Dark Frame?").grid(row=13, column=1, columnspan=2)
        self.subBiasDarkChoice = tk.StringVar()
        self.subBiasDarkComboBox = ttk.Combobox(self, textvariable=self.subBiasDarkChoice, values=('No', 'Yes'), state="readonly")
        self.subBiasDarkComboBox.grid(row=13, column=3, columnspan=2, sticky=tk.E + tk.W)
        self.subBiasDarkComboBox.current([photometry.getsubtractBiasFromDarkFlag()])

        # blankPerStarFlag setting
        ttk.Label(self, text="Background Subtraction Method:").grid(row=14, column=1, columnspan=2)
        self.blankChoice = tk.StringVar()
        self.blankComboBox = ttk.Combobox(self, textvariable=self.blankChoice, values=('Background from Entire Image', 'Background from Individual Stars'), state="readonly")
        self.blankComboBox.grid(row=14, column=3, columnspan=2, sticky=tk.E + tk.W)
        self.blankComboBox.current([photometry.getblankPerStarFlag()])

        # catalogChoice setting
        ttk.Label(self, text="Catalog to Search:").grid(row=15, column=1, columnspan=2)
        self.catalogEntry = ttk.Entry(self)
        self.catalogEntry.grid(row=15, column=3, columnspan=2, sticky=tk.E + tk.W)
        if not(photometry.getcatalogChoice() == 0):
            self.catalogEntry.insert(0,photometry.getcatalogChoice())
        else:
            self.catalogEntry.insert(0,"SIMBAD")

        # filterChoice setting
        ttk.Label(self, text="Filter to Search:").grid(row=16, column=1, columnspan=2)
        self.filterChoice = tk.StringVar()
        self.filterComboBox = ttk.Combobox(self, textvariable=self.filterChoice,
                                          values=('V', 'B', 'g', 'r','i', 'Other: Please type other filter'))
        self.filterComboBox.grid(row=16, column=3, columnspan=2, sticky=tk.E + tk.W)
        if photometry.getfilterChoice() == 'V':
            self.filterComboBox.current([0])
        if photometry.getfilterChoice() == 'B':
            self.filterComboBox.current([1])
        if photometry.getfilterChoice() == 'g':
            self.filterComboBox.current([2])
        if photometry.getfilterChoice() == 'r':
            self.filterComboBox.current([3])
        if photometry.getfilterChoice() == 'i':
            self.filterComboBox.current([4])

        # lightCurveLineFlag setting
        ttk.Label(self, text="Plot Light Curve with Line?").grid(row=17, column=1, columnspan=2)
        self.lineChoice = tk.StringVar()
        self.lineComboBox = ttk.Combobox(self, textvariable=self.lineChoice, values=('No', 'Yes'),
                                                state="readonly")
        self.lineComboBox.grid(row=17, column=3, columnspan=2, sticky=tk.E + tk.W)
        self.lineComboBox.current([photometry.getlightCurveLineFlag()])

        # errorChoice setting
        ttk.Label(self, text="How should error be calculated?").grid(row=18, column=1, columnspan=2)
        self.errorChoice = tk.StringVar()
        self.errorComboBox = ttk.Combobox(self, textvariable=self.errorChoice, values=('Standard Deviation', 'Weighted Magnitude', 'Jack Knife'),
                                                state="readonly")
        self.errorComboBox.grid(row=18, column=3, columnspan=2, sticky=tk.E + tk.W)
        if photometry.geterrorChoice() == 'STD':
            self.errorComboBox.current([0])
        if photometry.geterrorChoice() == 'WMG':
            self.errorComboBox.current([1])
        if photometry.geterrorChoice() == 'JKF':
            self.errorComboBox.current([2])

        # readInReferenceFlag
        ttk.Label(self, text="Use Reference Stars from File?").grid(row=19, column=1, columnspan=2)
        self.referenceChoice = tk.StringVar()
        self.referenceComboBox = ttk.Combobox(self, textvariable=self.referenceChoice, values=('No', 'Yes'),
                                         state="readonly")
        self.referenceComboBox.grid(row=19, column=3, columnspan=2, sticky=tk.E + tk.W)
        self.referenceComboBox.current([photometry.getreadInReferenceFlag()])
        if photometry.getreadInReferenceFlag() == 1:
            ttk.Label(self, text="Reference Stars").grid(row=20, column=1, columnspan=2)
            self.biasButton = ttk.Button(self, text="Select .csv File")
            self.biasButton.grid(row=20, column=3, columnspan=2)

        # readInRadiusFlag setting
        ttk.Label(self, text="Use Target Star Radius for All Stars?").grid(row=21, column=1, columnspan=2)
        self.radiusChoice = tk.StringVar()
        self.radiusComboBox = ttk.Combobox(self, textvariable=self.radiusChoice,
                                          values=('No', 'Yes'),
                                          state="readonly")
        self.radiusComboBox.grid(row=21, column=3, columnspan=2, sticky=tk.E + tk.W)
        if photometry.getreadInRadiusFlag() == 0:
            self.radiusComboBox.current([1])
        if photometry.getreadInRadiusFlag() == 1:
            self.radiusComboBox.current([0])

        # fwhmFlag setting
        ttk.Label(self, text="Target Star Radius:").grid(row=22, column=1, columnspan=2)
        self.fwhmChoice = tk.StringVar()
        self.fwhmComboBox = ttk.Combobox(self, textvariable=self.fwhmChoice,
                                           values=('Find Radius Manually', 'Use Full-Width Half-Maximum'),
                                           state="readonly")
        self.fwhmComboBox.grid(row=22, column=3, columnspan=2, sticky=tk.E + tk.W)
        self.fwhmComboBox.current([photometry.getfwhmFlag()])

        # astrometryTimeOutFlag setting
        ttk.Label(self, text="Iterations to wait for Astrometry.net:").grid(row=23, column=1, columnspan=2)
        self.timeOutEntry = ttk.Entry(self)
        self.timeOutEntry.grid(row=23, column=3, columnspan=2, sticky=tk.E + tk.W)
        if not(photometry.getastrometryTimeOutFlag() == 0):
            self.timeOutEntry.insert(0, str(photometry.getastrometryTimeOutFlag()))
        else:
            self.timeOutEntry.insert(0, "Until it Completes")

        # removeReferenceOutliersFlag setting
        ttk.Label(self, text="Remove Reference Star Outliers Over Z-score:").grid(row=24, column=1, columnspan=2)
        self.timeOutEntry = ttk.Entry(self)
        self.timeOutEntry.grid(row=24, column=3, columnspan=2, sticky=tk.E + tk.W)
        if not (photometry.getremoveReferenceOutliersFlag() == 0):
            self.timeOutEntry.insert(0, str(photometry.getremoveReferenceOutliersFlag()))
        else:
            self.timeOutEntry.insert(0, "Do Not Automatically Remove Outliers")

        # Submit
        self.submitButton = ttk.Button(self, text="Submit Project", command=self.submit)
        self.submitButton.grid(row=25, column=1, columnspan=4)

    def submit(self):
        self.check()
        continueAns = messagebox.askokcancel(title='Are you sure?', message='These settings cannot be changed later. Submit anyways?',
                                             icon="warning")
        if continueAns == False:
            return 0

    def check(self):
        errorMessage = ""
        # Check project name input
        name = self.projectName.get()
        if name == "" or len(name) < 1:
            errorMessage = errorMessage + "\nProject Name cannot be blank"
        else:
            path = os.path.isdir("PhotPSaveData/" + slugify(name)) #If path is true, need a different file name
            if path == True:
                errorMessage = errorMessage + "\nProject Name already exists, choose a unique project name"
        # Check coordinate input
        if self.coordinateChoice == "Decimal Degrees":
            if self.decimalDegreeRA.get() == "" or len(self.decimalDegreeRA.get()) < 1:
                errorMessage = errorMessage + "\nRight Ascension cannot be blank"
            if self.decimalDegreeDec.get() == "" or len(self.decimalDegreeDec.get()) < 1:
                errorMessage = errorMessage + "\nDeclination cannot be blank"
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

        # Check

    def gohome(self, event):
        self.controller.show_frame("HomePage")

    def getBias(self):
        self.bias = tk.filedialog.askopenfile()

    def getDark(self):
        self.dark = tk.filedialog.askopenfile()

    def getFlat(self):
        self.flat = tk.filedialog.askopenfile()