# spectral_emission_finder
This code finds faint emission signatures situated near a quasar continuum in a 2D spectra

GUI code to find faint emission signatures situated either on top of, or near the quasar continuum of a 2D spectra. The code will also compare the significance of the emission as compared to the quasar emission and find the spatial distance of the centroid of the emission to that of the quasar.  

This is a GUI code that takes in a spectra from a 2D spectroscopic source file and find faint emission signatures in a quasar spectra. The code will only work if the quasar continuum has previously been removed from the vicinity of the emission. For the purpose of this code, we will assume that the mentioned continuum subtraction has already been done on the source file.

The GUI features are currently available only for one rest frame emission transition. The code further lets you choose three regions - The region will possible emission inside and two left and right regions where the quasar flux is present. The code will further compare these three regions to check the significance of emission and finding the spatial impact parameter. The code can be run by simply downloading all the files in the folder and compiling the 'main_GUI.py' through a linux terminal (if all required modules are available, see below). The code is linux terminal friendly and has not been checked for other operating systems.

Required dependancies - This is a simple python code (version - 2.7). Due to issues in running matplotlib widget in version 3, unfortuantely the code might not work in python 3. However, experts are welcome to change the code and try. Here is the list of python modules required to run this code -

astropy

numpy

scipy

matplotlib

PyAstronomy

sys

pathlib

os

pyfits

csv

itertools

tabulate



Other essential files include -

initial_parameters.dat - This is the primary parameter file and is required to give initial values for many paramters used by the GUI. The file is self descriptive. Please have a look to change the default.

combined_J0025+1145_nir_OIII_5007.fits - This is an example of an acceptable format of the source file. It has to be a FITS file (preferably with ESO telescope compatible headers).

custom_functions.py - This is a file with all relevant custom functions defined which is utilised in the main GUI code ('main_GUI.py').



Here is a short description of interactive functions that one can perform and all the buttons that you see in the GUI -

Custom functions to play around with the 2D data - There are three buttons - 'Norm', 'Clip' and 'Smoothing', that can be used to customise the 2-D data until you make the emission clear enough to be seen by eye. 'Norm' is a scaling of the data, such that, flux_new = flux_original/(10^Norm). 'Clip' performs the numpy clip function on the data, such that, flux_new = np.clip(flux_original, -clip, +clip) and 'Smoothing' just performs a 2D gaussian smoothing of the data. Once a combination of these functions make the faint emission visible, one can further save all these parameters by cliking on the 'Save Params' button so that they dont have to go through this entire hassle again. Once these parameters are saved, it can be reloaded back instantly using the 'Get Params' button. To analyse this emission further, one has to declare windows in which the emission and the residual quasar continuum exist.      


Selecting 2-D regions for analysis - The quasar spectra might be noisy with emission faint features not immediately identifiable. For the purpose of analysing the faint emission, I have created a selection method to identify a large region where the emission might possibly lie. To do the same, one has to press the buttons 't' -> 'y' -> 'h' in that sequence. This activates the selection function and now one can draw a rectangle around the 2D data. Once the rectangle is drawn, one can specify using buttons whether the regions is a possible emission feature (button - 'Set Emission Region') or one of the quasar continuum regions (left or right, that can be declared by buttons - 'Set Left Region' and 'Set Right Region' consecutively). Once any of the three buttons are clicked, an intermediate figure will pop out of the GUI to show you the selected part. If this is not correct, one can repeat the process iteratively until one is satisfied with the selection regions. Once the regions are finalised, press the 'Save NIR Data' button to finally save the selections.

Once the regions and parameters are set, one has to click on 'Analyse Data' button. This will perform the comparison of flux, and fit the gaussian for identifying both the centroid of emission and that of the quasar to give an estimate of the spatial impact parameter in one direction (i.e., the spatial distance of the quasar line of sight from the center of emission). Further if the emission is resolved, one can find the width of the emission line. All these information will be printed out in the terminal after the data analysis is complete. The button also produces a nice figure showing the 1D and 2D emission. 

The code is free to use for all. Thank you for using the code and please feel free to contact me at - 'ranjan_adarsh@yahoo.com' for any comments, suggestions and intimation of any bugs.
