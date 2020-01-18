from __future__ import print_function, division


#astropy
from astropy.io import fits
from astropy.io.fits import getdata
from astropy.visualization import PercentileInterval, SqrtStretch
from astropy.wcs import WCS
from astropy.wcs import find_all_wcs
#from astropy import wcs
from astropy.cosmology import Planck15 as LCDM
import astropy.units as u
from astropy.io import fits as pf
from astropy.cosmology import FlatLambdaCDM
from astropy.modeling import models, fitting
import astropy.units as u


#numpy
import numpy
import numpy as np
import numpy.polynomial.chebyshev as cheb
import numpy.polynomial.polynomial as poly
import numpy.ma as ma
from numpy.polynomial import polynomial as P


#scipy
import scipy
import scipy as sp
import scipy.optimize as opt
import scipy.ndimage as ndimage
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy.optimize import fmin
from scipy import integrate
from scipy.integrate import quad
from scipy import signal
import scipy.ndimage
from scipy.stats import norm


#matplotlib
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons, RectangleSelector
#from matplotlib.widgets import TextBox
import matplotlib.patches as mpatches
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec

plt.rcParams["font.family"] = "Times New Roman"


#pyastronomy
from PyAstronomy import pyaGui
from PyAstronomy import pyasl
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy import funcFit as fuf


#others
import sys
from pathlib import Path
import os
import os.path
from os.path import exists
import pyfits
import csv as csv
import itertools
from tabulate import tabulate


#Cutom import
from custom_functions import *






#Getting source filenames and initial guess of the required paramters 
initial_guesses = np.genfromtxt('initial_parameters.dat', dtype=[('mystring','S50')], comments='#')
wl_type = str(initial_guesses[0][0])
size_of_font = int(initial_guesses[1][0])
file_name1 = str(initial_guesses[2][0])
redshift = float(initial_guesses[3][0])
fits_image_filename1 = file_name1
file_emission_data = str(initial_guesses[4][0])
arcsec_lim_lower = float(initial_guesses[5][0])
arcsec_lim_upper = float(initial_guesses[6][0])
wavelength_lim_lower = float(initial_guesses[7][0])
wavelength_lim_upper = float(initial_guesses[8][0])





#Initialising some constants and formatting variables
# Speed in Light in Km/s
c = 299792.458 

#Setting up cosmology
cosmo = FlatLambdaCDM(H0=67.8 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.308)

#Some random formatting code
underline = '\033[4m'
end_formatting = end_format = reset = '\033[0m'

#Initialising normalisation parameter
normalisation_constant = 1

#Setting GUI screen text font type and size 
font = {'family' : 'times new roman',
        #'weight' : 'bold',
        'size'   : size_of_font}

matplotlib.rc('font', **font)




######################CORRECTING FOR WAVELENGTH AND SPECIFYING POSITION ANGLE OF OBSERVATION##################################

def wavelength_solution_function(fits_file_name, redshift):
	fname = fits_file_name
	print ("")
	print (underline + "Working on file: "+fname + reset)
	image = pf.getdata(fname)*1.e17
	image = image[5:-5,:]
	header = pf.getheader(fname)
	posang = header['HIERARCH ESO ADA POSANG']
	wl = np.arange(image.shape[1])*header['CD1_1']+header['CRVAL1']
	vrad = header['ESO QC VRAD HELICOR']
	wl *= (1. + vrad/299792.)
	header['CD1_1'] *= 10.
	wl *= 10.

	### Convert air wavelength to vacuum:
	if wl_type == 'air':
		l_vac = air2vac(wl)
	elif wl_type == 'vac':
		l_vac = wl

	return(l_vac, posang)
	

######################CORRECTING FOR WAVELENGTH AND SPECIFYING POSITION ANGLE OF OBSERVATION##################################





# -- Dictionary of emission line and rest-frame wavelengths:
l_cen = {
    'OIII_5007': 5008.24, # [OIII] 5007 
}
line_list = sorted(l_cen.items(), key=lambda x: x[1])







#Loading data
data1, hdr1 = getdata(fits_image_filename1, 0, header=True)
data_err1, hdr_err1 = getdata(fits_image_filename1, 1, header=True)

#Initialising normalisation paramters
trish1 = data1
data1 = data1/normalisation_constant


#Defining wavelength array and other parameters obtained from the fits file
wavelength_solution_1, posang1 = wavelength_solution_function(fits_image_filename1, redshift)
posang11 = "PA = " + str(int(np.round(posang1, 0))) + str(r"$^{\circ}$")
arcsec_solution_1 = np.arange(data1.shape[0])*hdr1['CD2_2']+hdr1['CRVAL2']


#Defining limit for wavelength and arcsec axes to create a custom box for analysis
idx_arcsec_lim_lower_1 = np.searchsorted(arcsec_solution_1, float(arcsec_lim_lower))
idx_wavelength_lim_lower_1 = np.searchsorted(wavelength_solution_1, float(wavelength_lim_lower))
idx_arcsec_lim_upper_1 = np.searchsorted(arcsec_solution_1, float(arcsec_lim_upper))
idx_wavelength_lim_upper_1 = np.searchsorted(wavelength_solution_1, float(wavelength_lim_upper))




#Defining full range of x and y axes
y_axis = np.arange(0, data1.shape[0])
x_axis = np.arange(0, data1.shape[1])


start_value_x_axis_physical = (float(hdr1['CRVAL1'])) - ((int(hdr1['CRPIX1'])) * (float(hdr1['CDELT1'])))
start_value_y_axis_physical = (float(hdr1['CRVAL2'])) - ((int(hdr1['CRPIX2'])) * (float(hdr1['CDELT2'])))

x_axis_real = start_value_x_axis_physical
y_axis_real = start_value_y_axis_physical

for i in range(len(x_axis)-1):
	x_axis_real = np.append(x_axis_real, (start_value_x_axis_physical + ((i+1)*(float(hdr1['CDELT1'])))))
	
for i in range(len(y_axis)-1):
	y_axis_real = np.append(y_axis_real, (start_value_y_axis_physical + ((i+1)*(float(hdr1['CDELT2'])))))






#Plotting

fig, (ax1) = plt.subplots(1, sharex=True, sharey=True)
line1 = ax1.imshow(data1, cmap='Blues_r', aspect='equal', origin='lower')
ax_norm = plt.axes([0.10, 0.94, 0.3, 0.02])
snorm = Slider(ax_norm, 'Norm', (-20), (+10), valinit=-0.1)
ax_clipping = plt.axes([0.60, 0.94, 0.3, 0.02])
sclipping = Slider(ax_clipping, 'Clip', (0.0), (10.0), valinit=0.0)


#Defining smoothening parameter
global data_smooth
data_smooth = 1.0
axdata_smooth = plt.axes([0.6, 0.04, 0.3, 0.02])
sdata_smooth = Slider(axdata_smooth, 'Smoothing', (data_smooth+0), (data_smooth+5), valinit=data_smooth)


#updating values in GUI in realtime
def update(val):
	clip = sclipping.val
	normal = snorm.val
	smooth_sigma = sdata_smooth.val
	data_new1 = data1 / 10**(normal)
	data_new1 = np.clip(data_new1, -clip, +clip)
	sigma = [smooth_sigma, smooth_sigma]
	data_new1 = sp.ndimage.filters.gaussian_filter(data_new1, sigma, mode='constant')
	line1.set_data(data_new1)

fig.canvas.draw_idle()
sclipping.on_changed(update)
snorm.on_changed(update)
sdata_smooth.on_changed(update)


################################LOAD_AND_SAVE_SELECTION_POINTS################################

#Initialising some strings
text = ax1.text(0, 0, "")
str1 = ''



#Defining the selection function
def line_select_callback(eclick, erelease):
	global str1
	global x1, y1, x2, y2
	x1, y1 = eclick.xdata, eclick.ydata
	x2, y2 = erelease.xdata, erelease.ydata

	idx_yaxis_min = int(np.searchsorted(x_axis, min(x1, x2)))
	idx_yaxis_max = int(np.searchsorted(x_axis, max(x1, x2)))
	idx_xaxis_min = int(np.searchsorted(y_axis, min(y1, y2)))
	idx_xaxis_max = int(np.searchsorted(y_axis, max(y1, y2)))

	data_selected_region1 = []

	data_selected_emission1 = []
	
	rect_region = []
	rect_emission = []


	if (str1 == 'add'):
		data_selected_region1	= np.append(data_selected_region1, data1[idx_xaxis_min:idx_xaxis_max,idx_yaxis_min:idx_yaxis_max])
	

	elif (str1 == 'rem'):
		data_selected_emission1	= np.append(data_selected_emission1, data1[idx_xaxis_min:idx_xaxis_max,idx_yaxis_min:idx_yaxis_max])


	else:
        	print ("Function Inactive.....")


	# Create a Rectangle patch
	heigth_rect_region = idx_yaxis_max - idx_yaxis_min
	width_rect_region = idx_xaxis_max - idx_xaxis_min
	rect = patches.Rectangle((idx_xaxis_min,idx_yaxis_min),width_rect_region,heigth_rect_region,linewidth=1,edgecolor='y',facecolor='none')

	# Add the patch to the Axes
	ax1.add_patch(rect)

	fig.canvas.draw()



def toggle_selector(event):
	global str1
	print(' Key pressed.')
	if event.key in ['T', 't'] and toggle_selector.RS.active:
        	print(' RectangleSelector deactivated.')
        	toggle_selector.RS.set_active(False)
	if event.key in ['Y', 'y'] and not toggle_selector.RS.active:
        	print(' RectangleSelector activated.')
        	toggle_selector.RS.set_active(True)
	if event.key in ['H', 'h'] and toggle_selector.RS.active:
        	print('Add function activated')
        	str1 = 'add'
        	toggle_selector.RS.set_active(True)
	if event.key in ['J', 'j'] and toggle_selector.RS.active:
        	print('Remove function activated')
        	str1 = 'rem'
        	toggle_selector.RS.set_active(True)


toggle_selector.RS = RectangleSelector(ax1, line_select_callback, drawtype='box', useblit=False, button=[1], minspanx=5, minspany=5, spancoords='pixels', interactive=True)

plt.connect('key_press_event', toggle_selector)


################################LOAD_AND_SAVE_SELECTION_POINTS################################









#######################################SET_REGIONS_FOR_ANALYSIS##########################################


set_region_left_ax = plt.axes([0.001, 0.61, 0.08, 0.02])
button12 = Button(set_region_left_ax, 'Set Left Region')
def set_region_left_ax(event):



	x11, y11, x22, y22 = x1, y1, x2, y2


	
	nir_lr_idx_yaxis_min = int(np.searchsorted(x_axis, min(x11, x22)))
	nir_lr_idx_yaxis_max = int(np.searchsorted(x_axis, max(x11, x22)))
	nir_lr_idx_xaxis_min = int(np.searchsorted(y_axis, min(y11, y22)))
	nir_lr_idx_xaxis_max = int(np.searchsorted(y_axis, max(y11, y22)))

	global nir_lr_region_x_min, nir_lr_region_x_max, nir_lr_region_y_min, nir_lr_region_y_max
	nir_lr_region_x_min, nir_lr_region_x_max, nir_lr_region_y_min, nir_lr_region_y_max = nir_lr_idx_xaxis_min, nir_lr_idx_xaxis_max, nir_lr_idx_yaxis_min, nir_lr_idx_yaxis_max


	global nir_lr_data_selected_region_final1, nir_lr_data_selected_region_final2, nir_lr_data_selected_region_final3, nir_lr_err_selected_region_final1, nir_lr_err_selected_region_final2, nir_lr_err_selected_region_final3

	nir_lr_data_selected_region_final1 = []
	nir_lr_err_selected_region_final1 = []
	nir_lr_data_selected_region_final1 = data1[nir_lr_idx_xaxis_min:nir_lr_idx_xaxis_max,nir_lr_idx_yaxis_min:nir_lr_idx_yaxis_max]
	nir_lr_err_selected_region_final1 = data_err1[nir_lr_idx_xaxis_min:nir_lr_idx_xaxis_max,nir_lr_idx_yaxis_min:nir_lr_idx_yaxis_max]



	fig99, (ax_test91) = plt.subplots(1, 1, sharex=True)
	clip = sclipping.val
	normal = snorm.val
	smooth_sigma = sdata_smooth.val
	nir_lr_data_selected_region_final11 = nir_lr_data_selected_region_final1 / 10**(normal)
	nir_lr_data_selected_region_final11 = np.clip(nir_lr_data_selected_region_final11, -clip, +clip)
	sigma = [smooth_sigma, smooth_sigma]
	nir_lr_data_selected_region_final11 = sp.ndimage.filters.gaussian_filter(nir_lr_data_selected_region_final11, sigma, mode='constant')
	ax_test91.imshow(nir_lr_data_selected_region_final11)
	plt.show()

	a1, a2, b1, b2 = y11, y22, x11, x22

	main_rect_height = b2-b1
	main_rect_width = a2-a1

	rect_emission1 = patches.Rectangle(([b1,a1]),main_rect_height,main_rect_width,linewidth=1,edgecolor='y',facecolor='red', alpha=0.2)
	ax1.add_patch(rect_emission1)

	fig.canvas.draw()
	
	print ('Left Region is set!')
   
button12.on_clicked(set_region_left_ax)



set_emission_nir_ax = plt.axes([0.001, 0.56, 0.08, 0.02])
button13 = Button(set_emission_nir_ax, 'Set Emission Region')
def set_emission_nir_ax(event):
	x11, y11, x22, y22 = x1, y1, x2, y2
	
	nir_em_idx_yaxis_min = int(np.searchsorted(x_axis, min(x11, x22)))
	nir_em_idx_yaxis_max = int(np.searchsorted(x_axis, max(x11, x22)))
	nir_em_idx_xaxis_min = int(np.searchsorted(y_axis, min(y11, y22)))
	nir_em_idx_xaxis_max = int(np.searchsorted(y_axis, max(y11, y22)))

	global nir_em_region_x_min, nir_em_region_x_max, nir_em_region_y_min, nir_em_region_y_max
	nir_em_region_x_min, nir_em_region_x_max, nir_em_region_y_min, nir_em_region_y_max = nir_em_idx_xaxis_min, nir_em_idx_xaxis_max, nir_em_idx_yaxis_min, nir_em_idx_yaxis_max


	global nir_em_data_selected_region_final1, nir_em_data_selected_region_final2, nir_em_data_selected_region_final3, nir_em_err_selected_region_final1, nir_em_err_selected_region_final2, nir_em_err_selected_region_final3

	nir_em_data_selected_region_final1 = []
	nir_em_err_selected_region_final1 = []
	nir_em_data_selected_region_final1 = data1[nir_em_idx_xaxis_min:nir_em_idx_xaxis_max,nir_em_idx_yaxis_min:nir_em_idx_yaxis_max]
	nir_em_err_selected_region_final1 = data_err1[nir_em_idx_xaxis_min:nir_em_idx_xaxis_max,nir_em_idx_yaxis_min:nir_em_idx_yaxis_max]
	
	fig91, (ax_test92) = plt.subplots(1, 1, sharex=True)
	clip = sclipping.val
	normal = snorm.val
	smooth_sigma = sdata_smooth.val
	nir_em_data_selected_region_final11 = nir_em_data_selected_region_final1 / 10**(normal)
	nir_em_data_selected_region_final11 = np.clip(nir_em_data_selected_region_final11, -clip, +clip)
	sigma = [smooth_sigma, smooth_sigma]
	nir_em_data_selected_region_final11 = sp.ndimage.filters.gaussian_filter(nir_em_data_selected_region_final11, sigma, mode='constant')
	ax_test92.imshow(nir_em_data_selected_region_final11)
	plt.show()


	a1, a2, b1, b2 = y11, y22, x11, x22


	main_rect_height = b2-b1
	main_rect_width = a2-a1

	rect_emission1 = patches.Rectangle(([b1,a1]),main_rect_height,main_rect_width,linewidth=1,edgecolor='y',facecolor='red', alpha=0.4)
	ax1.add_patch(rect_emission1)

	fig.canvas.draw()

	print ('Emission Region is set!')
   
button13.on_clicked(set_emission_nir_ax)









set_region_right_ax = plt.axes([0.001, 0.51, 0.08, 0.02])
button14 = Button(set_region_right_ax, 'Set Right Region')
def set_region_right_ax(event):
	x11, y11, x22, y22 = x1, y1, x2, y2
	
	nir_rr_idx_yaxis_min = int(np.searchsorted(x_axis, min(x11, x22)))
	nir_rr_idx_yaxis_max = int(np.searchsorted(x_axis, max(x11, x22)))
	nir_rr_idx_xaxis_min = int(np.searchsorted(y_axis, min(y11, y22)))
	nir_rr_idx_xaxis_max = int(np.searchsorted(y_axis, max(y11, y22)))

	global nir_rr_region_x_min, nir_rr_region_x_max, nir_rr_region_y_min, nir_rr_region_y_max
	nir_rr_region_x_min, nir_rr_region_x_max, nir_rr_region_y_min, nir_rr_region_y_max = nir_rr_idx_xaxis_min, nir_rr_idx_xaxis_max, nir_rr_idx_yaxis_min, nir_rr_idx_yaxis_max



	global nir_rr_data_selected_region_final1, nir_rr_data_selected_region_final2, nir_rr_data_selected_region_final3, nir_rr_err_selected_region_final1, nir_rr_err_selected_region_final2, nir_rr_err_selected_region_final3

	nir_rr_data_selected_region_final1 = []
	nir_rr_err_selected_region_final1 = []
	nir_rr_data_selected_region_final1 = data1[nir_rr_idx_xaxis_min:nir_rr_idx_xaxis_max,nir_rr_idx_yaxis_min:nir_rr_idx_yaxis_max]
	nir_rr_err_selected_region_final1 = data_err1[nir_rr_idx_xaxis_min:nir_rr_idx_xaxis_max,nir_rr_idx_yaxis_min:nir_rr_idx_yaxis_max]
	
	fig93, (ax_test93) = plt.subplots(1, 1, sharex=True)
	clip = sclipping.val
	normal = snorm.val
	smooth_sigma = sdata_smooth.val
	nir_rr_data_selected_region_final11 = nir_rr_data_selected_region_final1 / 10**(normal)
	nir_rr_data_selected_region_final11 = np.clip(nir_rr_data_selected_region_final11, -clip, +clip)
	sigma = [smooth_sigma, smooth_sigma]
	nir_rr_data_selected_region_final11 = sp.ndimage.filters.gaussian_filter(nir_rr_data_selected_region_final11, sigma, mode='constant')
	ax_test93.imshow(nir_rr_data_selected_region_final11)
	plt.show()

	a1, a2, b1, b2 = y11, y22, x11, x22

	main_rect_height = b2-b1
	main_rect_width = a2-a1

	rect_emission1 = patches.Rectangle(([b1,a1]),main_rect_height,main_rect_width,linewidth=1,edgecolor='r',facecolor='red', alpha=0.2)
	ax1.add_patch(rect_emission1)

	fig.canvas.draw()

	print ('Right Region is set!')
   
button14.on_clicked(set_region_right_ax)

	










#Save and load parameters
save_nir_ax = plt.axes([0.001, 0.46, 0.08, 0.02])
button16 = Button(save_nir_ax, 'Save NIR Data')
def save_nir_ax(event):


	global nir_lr_region_x_min2, nir_lr_region_x_max2, nir_lr_region_y_min2, nir_lr_region_y_max2
	global nir_em_region_x_min2, nir_em_region_x_max2, nir_em_region_y_min2, nir_em_region_y_max2
	global nir_rr_region_x_min2, nir_rr_region_x_max2, nir_rr_region_y_min2, nir_rr_region_y_max2

	nir_lr_region_x_min2, nir_lr_region_x_max2, nir_lr_region_y_min2, nir_lr_region_y_max2 = nir_lr_region_x_min, nir_lr_region_x_max, nir_lr_region_y_min, nir_lr_region_y_max

	nir_em_region_x_min2, nir_em_region_x_max2, nir_em_region_y_min2, nir_em_region_y_max2 = nir_em_region_x_min, nir_em_region_x_max, nir_em_region_y_min, nir_em_region_y_max

	nir_rr_region_x_min2, nir_rr_region_x_max2, nir_rr_region_y_min2, nir_rr_region_y_max2 = nir_rr_region_x_min, nir_rr_region_x_max, nir_rr_region_y_min, nir_rr_region_y_max

	array_for_nir_data = np.array([nir_lr_region_x_min2, nir_lr_region_x_max2, nir_lr_region_y_min2, nir_lr_region_y_max2, nir_em_region_x_min2, nir_em_region_x_max2, nir_em_region_y_min2, nir_em_region_y_max2, nir_rr_region_x_min2, nir_rr_region_x_max2, nir_rr_region_y_min2, nir_rr_region_y_max2])
	
	np.savetxt('emission_data_boxes.txt', array_for_nir_data)
		
	print ('Saved!')
   
button16.on_clicked(save_nir_ax)



get_nir_ax = plt.axes([0.001, 0.41, 0.08, 0.02])
button17 = Button(get_nir_ax, 'Get NIR Data')
def get_nir_ax(event):


	global nir_lr_region_x_min2, nir_lr_region_x_max2, nir_lr_region_y_min2, nir_lr_region_y_max2
	global nir_em_region_x_min2, nir_em_region_x_max2, nir_em_region_y_min2, nir_em_region_y_max2
	global nir_rr_region_x_min2, nir_rr_region_x_max2, nir_rr_region_y_min2, nir_rr_region_y_max2

	nir_lr_region_x_min2, nir_lr_region_x_max2, nir_lr_region_y_min2, nir_lr_region_y_max2, nir_em_region_x_min2, nir_em_region_x_max2, nir_em_region_y_min2, nir_em_region_y_max2, nir_rr_region_x_min2, nir_rr_region_x_max2, nir_rr_region_y_min2, nir_rr_region_y_max2 = np.loadtxt('emission_data_boxes.txt', unpack=True)

	

	global nir_lr_region_x_min, nir_lr_region_x_max, nir_lr_region_y_min, nir_lr_region_y_max
	nir_lr_region_x_min, nir_lr_region_x_max, nir_lr_region_y_min, nir_lr_region_y_max = nir_lr_region_x_min2, nir_lr_region_x_max2, nir_lr_region_y_min2, nir_lr_region_y_max2
	

	nir_lr_idx_xaxis_min, nir_lr_idx_xaxis_max, nir_lr_idx_yaxis_min, nir_lr_idx_yaxis_max = int(nir_lr_region_x_min), int(nir_lr_region_x_max), int(nir_lr_region_y_min), int(nir_lr_region_y_max)


	global nir_lr_data_selected_region_final1, nir_lr_data_selected_region_final2, nir_lr_data_selected_region_final3, nir_lr_err_selected_region_final1, nir_lr_err_selected_region_final2, nir_lr_err_selected_region_final3

	nir_lr_data_selected_region_final1 = []
	nir_lr_err_selected_region_final1 = []
	nir_lr_data_selected_region_final1 = data1[nir_lr_idx_xaxis_min:nir_lr_idx_xaxis_max,nir_lr_idx_yaxis_min:nir_lr_idx_yaxis_max]
	nir_lr_err_selected_region_final1 = data_err1[nir_lr_idx_xaxis_min:nir_lr_idx_xaxis_max,nir_lr_idx_yaxis_min:nir_lr_idx_yaxis_max]

	a1, a2, b1, b2 = nir_lr_idx_xaxis_min, nir_lr_idx_xaxis_max, nir_lr_idx_yaxis_min, nir_lr_idx_yaxis_max

	main_rect_height = b2-b1
	main_rect_width = a2-a1

	rect_emission1 = patches.Rectangle(([b1,a1]),main_rect_height,main_rect_width,linewidth=1,edgecolor='y',facecolor='red', alpha=0.2)
	ax1.add_patch(rect_emission1)

	fig.canvas.draw()
	
	print ('Left Region is set!')


	global nir_em_region_x_min, nir_em_region_x_max, nir_em_region_y_min, nir_em_region_y_max
	nir_em_region_x_min, nir_em_region_x_max, nir_em_region_y_min, nir_em_region_y_max = nir_em_region_x_min2, nir_em_region_x_max2, nir_em_region_y_min2, nir_em_region_y_max2
	

	nir_em_idx_xaxis_min, nir_em_idx_xaxis_max, nir_em_idx_yaxis_min, nir_em_idx_yaxis_max = int(nir_em_region_x_min), int(nir_em_region_x_max), int(nir_em_region_y_min), int(nir_em_region_y_max)

	global nir_em_data_selected_region_final1, nir_em_data_selected_region_final2, nir_em_data_selected_region_final3, nir_em_err_selected_region_final1, nir_em_err_selected_region_final2, nir_em_err_selected_region_final3

	nir_em_data_selected_region_final1 = []
	nir_em_err_selected_region_final1 = []
	nir_em_data_selected_region_final1 = data1[nir_em_idx_xaxis_min:nir_em_idx_xaxis_max,nir_em_idx_yaxis_min:nir_em_idx_yaxis_max]
	nir_em_err_selected_region_final1 = data_err1[nir_em_idx_xaxis_min:nir_em_idx_xaxis_max,nir_em_idx_yaxis_min:nir_em_idx_yaxis_max]

	a1, a2, b1, b2 = nir_em_idx_xaxis_min, nir_em_idx_xaxis_max, nir_em_idx_yaxis_min, nir_em_idx_yaxis_max

	main_rect_height = b2-b1
	main_rect_width = a2-a1

	rect_emission1 = patches.Rectangle(([b1,a1]),main_rect_height,main_rect_width,linewidth=1,edgecolor='y',facecolor='red', alpha=0.4)
	ax1.add_patch(rect_emission1)

	fig.canvas.draw()


	print ('Emission Region is set!')








	global nir_rr_region_x_min, nir_rr_region_x_max, nir_rr_region_y_min, nir_rr_region_y_max
	nir_rr_region_x_min, nir_rr_region_x_max, nir_rr_region_y_min, nir_rr_region_y_max = nir_rr_region_x_min2, nir_rr_region_x_max2, nir_rr_region_y_min2, nir_rr_region_y_max2
	

	nir_rr_idx_xaxis_min, nir_rr_idx_xaxis_max, nir_rr_idx_yaxis_min, nir_rr_idx_yaxis_max = int(nir_rr_region_x_min), int(nir_rr_region_x_max), int(nir_rr_region_y_min), int(nir_rr_region_y_max)

	global nir_rr_data_selected_region_final1, nir_rr_data_selected_region_final2, nir_rr_data_selected_region_final3, nir_rr_err_selected_region_final1, nir_rr_err_selected_region_final2, nir_rr_err_selected_region_final3


	
	nir_rr_region_x_min, nir_rr_region_x_max, nir_rr_region_y_min, nir_rr_region_y_max = nir_rr_idx_xaxis_min, nir_rr_idx_xaxis_max, nir_rr_idx_yaxis_min, nir_rr_idx_yaxis_max



	global nir_rr_data_selected_region_final1, nir_rr_data_selected_region_final2, nir_rr_data_selected_region_final3, nir_rr_err_selected_region_final1, nir_rr_err_selected_region_final2, nir_rr_err_selected_region_final3

	nir_rr_data_selected_region_final1 = []
	nir_rr_err_selected_region_final1 = []
	nir_rr_data_selected_region_final1 = data1[nir_rr_idx_xaxis_min:nir_rr_idx_xaxis_max,nir_rr_idx_yaxis_min:nir_rr_idx_yaxis_max]
	nir_rr_err_selected_region_final1 = data_err1[nir_rr_idx_xaxis_min:nir_rr_idx_xaxis_max,nir_rr_idx_yaxis_min:nir_rr_idx_yaxis_max]

	a1, a2, b1, b2 = nir_rr_idx_xaxis_min, nir_rr_idx_xaxis_max, nir_rr_idx_yaxis_min, nir_rr_idx_yaxis_max

	main_rect_height = b2-b1
	main_rect_width = a2-a1

	rect_emission1 = patches.Rectangle(([b1,a1]),main_rect_height,main_rect_width,linewidth=1,edgecolor='r',facecolor='red', alpha=0.2)
	ax1.add_patch(rect_emission1)

	fig.canvas.draw()


	print ('Right Region is set!')

	print ('Loaded!')
   
button17.on_clicked(get_nir_ax)










#Analyse the emission data. After analysis, the function will print values of flux and impact parameter on the terminal and save figure 
analyse_nir_ax = plt.axes([0.001, 0.36, 0.08, 0.02])
button15 = Button(analyse_nir_ax, 'Analyse Data')
def analyse_nir_ax(event):

	global nir_lr_region_x_min2, nir_lr_region_x_max2, nir_lr_region_y_min2, nir_lr_region_y_max2
	global nir_em_region_x_min2, nir_em_region_x_max2, nir_em_region_y_min2, nir_em_region_y_max2
	global nir_rr_region_x_min2, nir_rr_region_x_max2, nir_rr_region_y_min2, nir_rr_region_y_max2

	nir_lr_region_x_min2, nir_lr_region_x_max2, nir_lr_region_y_min2, nir_lr_region_y_max2, nir_em_region_x_min2, nir_em_region_x_max2, nir_em_region_y_min2, nir_em_region_y_max2, nir_rr_region_x_min2, nir_rr_region_x_max2, nir_rr_region_y_min2, nir_rr_region_y_max2 = np.loadtxt('emission_data_boxes.txt', unpack=True)



	######################LEFT######################

	global nir_lr_region_x_min, nir_lr_region_x_max, nir_lr_region_y_min, nir_lr_region_y_max
	nir_lr_region_x_min, nir_lr_region_x_max, nir_lr_region_y_min, nir_lr_region_y_max = nir_lr_region_x_min2, nir_lr_region_x_max2, nir_lr_region_y_min2, nir_lr_region_y_max2
	

	nir_lr_idx_xaxis_min, nir_lr_idx_xaxis_max, nir_lr_idx_yaxis_min, nir_lr_idx_yaxis_max = int(nir_lr_region_x_min), int(nir_lr_region_x_max), int(nir_lr_region_y_min), int(nir_lr_region_y_max)


	global nir_lr_data_selected_region_final1, nir_lr_data_selected_region_final2, nir_lr_data_selected_region_final3, nir_lr_err_selected_region_final1, nir_lr_err_selected_region_final2, nir_lr_err_selected_region_final3

	nir_lr_data_selected_region_final1 = []
	nir_lr_err_selected_region_final1 = []

	######################LEFT######################






	######################EMISSION######################

	global nir_em_region_x_min, nir_em_region_x_max, nir_em_region_y_min, nir_em_region_y_max
	nir_em_region_x_min, nir_em_region_x_max, nir_em_region_y_min, nir_em_region_y_max = nir_em_region_x_min2, nir_em_region_x_max2, nir_em_region_y_min2, nir_em_region_y_max2
	

	nir_em_idx_xaxis_min, nir_em_idx_xaxis_max, nir_em_idx_yaxis_min, nir_em_idx_yaxis_max = int(nir_em_region_x_min), int(nir_em_region_x_max), int(nir_em_region_y_min), int(nir_em_region_y_max)

	global nir_em_data_selected_region_final1, nir_em_data_selected_region_final2, nir_em_data_selected_region_final3, nir_em_err_selected_region_final1, nir_em_err_selected_region_final2, nir_em_err_selected_region_final3

	nir_em_data_selected_region_final1 = []
	nir_em_err_selected_region_final1 = []

	######################EMISSION######################




	######################RIGHT######################

	global nir_rr_region_x_min, nir_rr_region_x_max, nir_rr_region_y_min, nir_rr_region_y_max
	nir_rr_region_x_min, nir_rr_region_x_max, nir_rr_region_y_min, nir_rr_region_y_max = nir_rr_region_x_min2, nir_rr_region_x_max2, nir_rr_region_y_min2, nir_rr_region_y_max2
	

	nir_rr_idx_xaxis_min, nir_rr_idx_xaxis_max, nir_rr_idx_yaxis_min, nir_rr_idx_yaxis_max = int(nir_rr_region_x_min), int(nir_rr_region_x_max), int(nir_rr_region_y_min), int(nir_rr_region_y_max)

	global nir_rr_data_selected_region_final1, nir_rr_data_selected_region_final2, nir_rr_data_selected_region_final3, nir_rr_err_selected_region_final1, nir_rr_err_selected_region_final2, nir_rr_err_selected_region_final3


	
	nir_rr_region_x_min, nir_rr_region_x_max, nir_rr_region_y_min, nir_rr_region_y_max = nir_rr_idx_xaxis_min, nir_rr_idx_xaxis_max, nir_rr_idx_yaxis_min, nir_rr_idx_yaxis_max



	global nir_rr_data_selected_region_final1, nir_rr_data_selected_region_final2, nir_rr_data_selected_region_final3, nir_rr_err_selected_region_final1, nir_rr_err_selected_region_final2, nir_rr_err_selected_region_final3

	nir_rr_data_selected_region_final1 = []
	nir_rr_err_selected_region_final1 = []

	######################RIGHT######################




	idx_arcsec_lim_lower_1_local_qso = min(nir_lr_idx_xaxis_min, nir_rr_idx_xaxis_min)
	idx_arcsec_lim_upper_1_local_qso = min(nir_lr_idx_xaxis_max, nir_rr_idx_xaxis_max)


	nir_lr_data_selected_region_final1 = data1[idx_arcsec_lim_lower_1_local_qso:idx_arcsec_lim_upper_1_local_qso,nir_lr_idx_yaxis_min:nir_lr_idx_yaxis_max]

	nir_lr_err_selected_region_final1 = data_err1[idx_arcsec_lim_lower_1_local_qso:idx_arcsec_lim_upper_1_local_qso,nir_lr_idx_yaxis_min:nir_lr_idx_yaxis_max]


	idx_arcsec_lim_lower_1_local_em = nir_em_idx_xaxis_min
	idx_arcsec_lim_upper_1_local_em = nir_em_idx_xaxis_max


	nir_em_data_selected_region_final1 = data1[idx_arcsec_lim_lower_1_local_em:idx_arcsec_lim_upper_1_local_em,nir_em_idx_yaxis_min:nir_em_idx_yaxis_max]

	nir_em_err_selected_region_final1 = data_err1[idx_arcsec_lim_lower_1_local_em:idx_arcsec_lim_upper_1_local_em,nir_em_idx_yaxis_min:nir_em_idx_yaxis_max]




	nir_rr_data_selected_region_final1 = data1[idx_arcsec_lim_lower_1_local_qso:idx_arcsec_lim_upper_1_local_qso,nir_rr_idx_yaxis_min:nir_rr_idx_yaxis_max]

	nir_rr_err_selected_region_final1 = data_err1[idx_arcsec_lim_lower_1_local_qso:idx_arcsec_lim_upper_1_local_qso,nir_rr_idx_yaxis_min:nir_rr_idx_yaxis_max]
	

	################################################FIT_SPECTRAL_SIDE################################################


	data4 = data1
	data_spectral_4 = np.zeros([data1.shape[1]])
	data_spectral_err_4 = np.zeros([data1.shape[1]]) 
	for i in range(data4.shape[1]):
		data_spectral_4[i] = np.sum(data1[idx_arcsec_lim_lower_1_local_em:idx_arcsec_lim_upper_1_local_em, i])
		data_spectral_err_4[i] = np.sqrt(np.sum(data_err1[idx_arcsec_lim_lower_1_local_em:idx_arcsec_lim_upper_1_local_em, i]**2))

	

	x = wavelength_solution_1[nir_em_idx_yaxis_min:nir_em_idx_yaxis_max]
	y = data_spectral_4[nir_em_idx_yaxis_min:nir_em_idx_yaxis_max]
	err = data_spectral_err_4[nir_em_idx_yaxis_min:nir_em_idx_yaxis_max]

	popt_emm, perr_emm = gaus_fitting(x, y, err)
	amp_emm, mean_emm, std_emm = popt_emm[:]
	std_emm = np.abs(std_emm)
	amp_err_emm, mean_err_emm, std_err_emm = perr_emm[:]
	emitter_redshift = ((float(mean_emm)) / (line_list[0][1])) - 1

	fwhm_emm = 2 * np.sqrt(2 * np.log(2)) * np.abs(std_emm) 
	fwhm_err_emm = np.abs(2 * np.sqrt(2 * np.log(2))) * std_err_emm 
	H_par_emm = (0.3989 * amp_emm)/std_emm
	H_par_err_emm = np.abs(H_par_emm) * np.sqrt ( (amp_err_emm/amp_emm)**2 + (std_err_emm/std_emm)**2 )
	fitted_area_emm = (float(H_par_emm) * float(fwhm_emm))/(2.35*0.3989)
	fitted_area_err_emm = np.abs(fitted_area_emm) * np.sqrt( (float(H_par_err_emm)/float(H_par_emm))**2 + (float(fwhm_err_emm)/float(fwhm_emm))**2 )

	print ('\n\n')
	print ('Emission - ')

	print ('Area from fit emission spectral -')
	fitted_area_emm, fitted_area_err_emm = error_on_gaussian_fit(x, amp_emm, amp_err_emm, mean_emm, mean_err_emm, std_emm, std_err_emm)
	fitted_area_emm = np.sum(gaus(x, *popt_emm)) * (x[1]-x[0])
	individual_sigma_area_emm = np.std(y - (gaus(x, *popt_emm)))
	fitted_area_err_emm = np.sqrt((2*fwhm_emm)/(x[1]-x[0]))*individual_sigma_area_emm 
	print ('Area by residual method emission  -', fitted_area_emm, fitted_area_err_emm)

	print ('Amplitude - ', (float(amp_emm)), "+/-", (float(amp_err_emm)))
	print ('FHWM - ', (float(fwhm_emm)), "+/-", (float(fwhm_err_emm)))
	print ('Center - ', (float(mean_emm)), "+/-", (float(mean_err_emm)))
	print ('Flux (1e-17) - ', (float(fitted_area_emm)), "+/-", (float(fitted_area_err_emm)))

	dL = ((cosmo.luminosity_distance(emitter_redshift)).value)*3.086e+24
	luminosity = ((fitted_area_emm*(4*np.pi*(dL**2)))*1e-17)/1e40	
	luminosity_err = ((fitted_area_err_emm*(4*np.pi*(dL**2)))*1e-17)/1e40
	print ('Luminosity (1e40) - ', (float(luminosity)), "+/-", (float(luminosity_err)))



	if ((line_list[0][0])=='Halpha'):
		dL = ((cosmo.luminosity_distance(emitter_redshift)).value)*3.086e+24
		luminosity = fitted_area_emm*4*np.pi*dL**2
		luminosity_err = fitted_area_err_emm*4*np.pi*dL**2
		SFR = luminosity*1.4e-41
		SFR_err = luminosity_err*1.4e-41
		SFR = SFR/1.59
		SFR_err = SFR_err/1.59
		print ('SFR - ', (float(SFR)), "+/-", (float(SFR_err)))
        elif ((line_list[0][0])=='OII_3727'):
		dL = ((cosmo.luminosity_distance(emitter_redshift)).value)*3.086e+24
		luminosity = fitted_area_emm*4*np.pi*dL**2
		luminosity_err = fitted_area_err_emm*4*np.pi*dL**2
		SFR = luminosity*1.4e-41
		SFR_err = luminosity_err*1.4e-41
		SFR = SFR/1.59
		SFR_err = SFR_err/1.59
		print ('SFR - ', (np.round(float(SFR)), 3), "+/-", (np.round(float(SFR_err)), 3))

	
	print ('Emitter Redshift - ', emitter_redshift)

	################################################FIT_SPECTRAL_SIDE################################################






	################################################PLOTTING_2D_IMAGE_DATA################################################


	clip = sclipping.val
	normal = snorm.val
	smooth_sigma = sdata_smooth.val
	data_new1 = data1 / 10**(normal)
	data_new1 = np.clip(data_new1, -clip, +clip)
	sigma = [smooth_sigma, smooth_sigma]
	data_new1 = sp.ndimage.filters.gaussian_filter(data_new1, sigma, mode='constant')


	################################################PLOTTING_2D_IMAGE_DATA################################################





	################################################FIRST_WINDOW################################################


	data_spatial_qso_left_1 = np.zeros([data1.shape[0]])
	data_spatial_qso_err_left_1 = np.zeros([data1.shape[0]]) 
	for i in range(data1.shape[0]):
		data_spatial_qso_left_1[i] = np.sum(data1[i, nir_lr_idx_yaxis_min:nir_lr_idx_yaxis_max])
		data_spatial_qso_err_left_1[i] = np.sqrt(np.sum(data_err1[i, nir_lr_idx_yaxis_min:nir_lr_idx_yaxis_max]**2))

	x1 = arcsec_solution_1[idx_arcsec_lim_lower_1_local_qso:idx_arcsec_lim_upper_1_local_qso]
	y1 = data_spatial_qso_left_1[idx_arcsec_lim_lower_1_local_qso:idx_arcsec_lim_upper_1_local_qso]  
	err1 = data_spatial_qso_err_left_1[idx_arcsec_lim_lower_1_local_qso:idx_arcsec_lim_upper_1_local_qso]
	popt_qso_left_1, perr_qso_left_1 = gaus_fitting(x1, y1, err1)
	amp_qso_left_1, mean_qso_left_1, std_qso_left_1 = popt_qso_left_1[:]
	std_qso_left_1 = np.abs(std_qso_left_1)
	amp_err_qso_left_1, mean_err_qso_left_1, std_err_qso_left_1 = perr_qso_left_1[:]
	center_qso_left_1 = (float(mean_qso_left_1))



	data_spatial_em_left_1 = np.zeros([data1.shape[0]])
	data_spatial_em_err_left_1 = np.zeros([data1.shape[0]]) 
	for i in range(data1.shape[0]):
		data_spatial_em_left_1[i] = np.sum(data1[i, nir_em_idx_yaxis_min:nir_em_idx_yaxis_max])
		data_spatial_em_err_left_1[i] = np.sqrt(np.sum(data_err1[i, nir_em_idx_yaxis_min:nir_em_idx_yaxis_max]**2))

	x_em_1 = arcsec_solution_1[idx_arcsec_lim_lower_1_local_em:idx_arcsec_lim_upper_1_local_em]
	y_em_1 = data_spatial_em_left_1[idx_arcsec_lim_lower_1_local_em:idx_arcsec_lim_upper_1_local_em]  
	err_em_1 = data_spatial_em_err_left_1[idx_arcsec_lim_lower_1_local_em:idx_arcsec_lim_upper_1_local_em]
	popt_em_1, perr_em_1 = gaus_fitting(x_em_1, y_em_1, err_em_1)
	amp_em_1, mean_em_1, std_em_1 = popt_em_1[:]
	std_em_1 = np.abs(std_em_1)
	amp_err_em_1, mean_err_em_1, std_err_em_1 = perr_em_1[:]
	center_em_1 = (float(mean_em_1))


	################################################FIRST_WINDOW################################################





	######################FIRST_WINDOW_IMPACT_PARAMETER_AND_GALAXY_FWHM########################



	fwhm_qso_1 = 2 * np.sqrt(2 * np.log(2)) * np.abs(std_qso_left_1) 
	fwhm_err_qso_1 = np.abs(2 * np.sqrt(2 * np.log(2))) * np.abs(std_err_qso_left_1) 
	H_par_qso_1 = (0.3989 * amp_qso_left_1)/std_qso_left_1
	H_par_qso_err_1 = np.abs(H_par_qso_1) * np.sqrt ( (amp_err_qso_left_1/amp_qso_left_1)**2 + (std_err_qso_left_1/std_qso_left_1)**2 )
	fitted_area_qso_1 = (float(H_par_qso_1) * float(fwhm_qso_1))/(2.35*0.3989)
	fitted_area_err_qso_1 = np.abs(fitted_area_qso_1) * np.sqrt( (float(H_par_qso_err_1)/float(H_par_qso_1))**2 + (float(fwhm_err_qso_1)/float(fwhm_qso_1))**2 )

	
	fwhm_em_1 = 2 * np.sqrt(2 * np.log(2)) * np.abs(std_em_1) 
	fwhm_err_em_1 = np.abs(2 * np.sqrt(2 * np.log(2))) * np.abs(std_err_em_1) 
	H_par_em_1 = (0.3989 * amp_em_1)/std_em_1
	H_par_em_err_1 = np.abs(H_par_em_1) * np.sqrt ( (amp_err_em_1/amp_em_1)**2 + (std_err_em_1/std_em_1)**2 )
	fitted_area_em_1 = (float(H_par_em_1) * float(fwhm_em_1))/(2.35*0.3989)
	fitted_area_err_em_1 = np.abs(fitted_area_em_1) * np.sqrt( (float(H_par_em_err_1)/float(H_par_em_1))**2 + (float(fwhm_err_em_1)/float(fwhm_em_1))**2 )



	print ('\n')
	print ('Emission 1 - ')

	print ('Area from fit emission 1 -')
	fitted_area_em_1, fitted_area_err_em_1 = error_on_gaussian_fit(x_em_1, amp_em_1, amp_err_em_1, mean_em_1, mean_err_em_1, std_em_1, std_err_em_1)
	fitted_area_em_1 = np.sum(gaus(x_em_1, *popt_em_1)) * (x_em_1[1]-x_em_1[0])
	individual_sigma_area_em_1 = np.std(y_em_1 - (gaus(x_em_1, *popt_em_1)))
	fitted_area_err_em_1 = np.sqrt((2*fwhm_em_1)/(x_em_1[1]-x_em_1[0]))*individual_sigma_area_em_1 
	print ('Area by residual method for emission 1 -', fitted_area_em_1, fitted_area_err_em_1)

	print ('Area from fit qso 1 -')
	fitted_area_qso_1, fitted_area_err_qso_1 = error_on_gaussian_fit(x1, amp_qso_left_1, amp_err_qso_left_1, mean_qso_left_1, mean_err_qso_left_1, std_qso_left_1, std_err_qso_left_1)
	fitted_area_qso_1 = np.sum(gaus(x1, *popt_qso_left_1)) * (x1[1]-x1[0])
	individual_sigma_area_qso_1 = np.std(y1 - (gaus(x1, *popt_qso_left_1)))
	fitted_area_err_qso_1 = np.sqrt((2*fwhm_qso_1)/(x1[1]-x1[0]))*individual_sigma_area_qso_1 
	print ('Area by residual method qso 1 -', fitted_area_qso_1, fitted_area_err_qso_1)

	print ('Position Angle - ', str(posang11))
	print ('Amplitude - ', (float(amp_em_1)), "+/-", (float(amp_err_em_1)))
	print ('FHWM qso (arcsec) - ', (float(fwhm_qso_1)), "+/-", (float(fwhm_err_qso_1)))
	print ('FHWM emission (arcsec) - ', (float(fwhm_em_1)), "+/-", (float(fwhm_err_em_1)))
	print ('Center - ', (float(mean_em_1)), "+/-", (float(mean_err_em_1)))
	print ('Area - ', (float(fitted_area_em_1)), "+/-", (float(fitted_area_err_em_1)))



	if ((float(fwhm_qso_1+fwhm_err_qso_1)) < (float(fwhm_em_1-fwhm_err_em_1))):
		galaxy_width1 = np.sqrt( (float(fwhm_em_1))**2 - (float(fwhm_qso_1))**2 )
		galaxy_width1_err = ( np.sqrt((4 * (float(fwhm_em_1))**2 * (float(fwhm_err_em_1))**2 ) + (4 * (float(fwhm_qso_1))**2 * (float(fwhm_err_qso_1))**2 ))  * np.sqrt( galaxy_width1 ) )/2.
		
		print ('Emission seems to be resolved in window-1 FWHM (arcsec) - ', galaxy_width1, '+/-', galaxy_width1_err)
		scale_distance_gal_1 = ((cosmo.angular_diameter_distance(emitter_redshift).value)/206.264806) * galaxy_width1
		scale_distance_gal_err_1 = ((cosmo.angular_diameter_distance(emitter_redshift).value)/206.264806) * galaxy_width1_err
		print ('Emission seems to be resolved in window-1. FWHM (kpc) - ', scale_distance_gal_1, '+/-', scale_distance_gal_err_1)



	impact_par_arcsec_1 = (np.abs((float(mean_em_1)) - (float(mean_qso_left_1))))
	impact_par_arcsec_err_1 = np.sqrt( (mean_err_em_1)**2 + (mean_err_qso_left_1)**2 )
	print ('Impact Parameter (arcsec) 1 - ', impact_par_arcsec_1, '+/-', impact_par_arcsec_err_1)
	scale_distance_1 = ((cosmo.angular_diameter_distance(emitter_redshift).value)/206.264806) * impact_par_arcsec_1
	scale_distance_err_1 = ((cosmo.angular_diameter_distance(emitter_redshift).value)/206.264806) * impact_par_arcsec_err_1
	print ('Impact Parameter (kpc) 1 - ', scale_distance_1, '+/-', scale_distance_err_1)




	if ((line_list[0][0])=='Halpha'):
		#dL = LCDM.luminosity_distance(emitter_redshift).to('cm')
		dL = ((cosmo.luminosity_distance(emitter_redshift)).value)*3.086e+24
		luminosity = (fitted_area_em_1*4*np.pi*dL**2)*1e-17
		luminosity_err = (fitted_area_err_em_1*4*np.pi*dL**2)*1e-17
		SFR = luminosity*1.4e-41
		SFR_err = luminosity_err*1.4e-41
		SFR = SFR/1.59
		SFR_err = SFR_err/1.59
		print ('SFR at 1 - ', (float(SFR)), "+/-", (float(SFR_err)))
        elif ((line_list[0][0])=='OII_3727'):
		#dL = LCDM.luminosity_distance(emitter_redshift).to('cm')
		dL = ((cosmo.luminosity_distance(emitter_redshift)).value)*3.086e+24
		luminosity = fitted_area_em_1*4*np.pi*dL**2
		luminosity_err = fitted_area_err_em_1*4*np.pi*dL**2
		SFR = luminosity*1.4e-41
		SFR_err = luminosity_err*1.4e-41
		SFR = SFR/1.59
		SFR_err = SFR_err/1.59
		print ('SFR at 1- ', (np.round(float(SFR)), 3), "+/-", (np.round(float(SFR_err)), 3))


											     


	########################FIRST_WINDOW_IMPACT_PARAMETER_AND_GALAXY_FWHM#############################








	################################################PLOT################################################

	fig88 = plt.figure(figsize=[10,5])
	grid = gridspec.GridSpec(6, 6, wspace=0.0, hspace=0.0)

	line12 = fig88.add_subplot(grid[1:3, :-1])
	line11 = fig88.add_subplot(grid[1:3, 5])
	line52 = fig88.add_subplot(grid[4:6, :-1])
	cmap = 'rainbow'

	line12.imshow(data_new1, origin='lower', interpolation='none', aspect='auto', extent=[wavelength_solution_1[0], wavelength_solution_1[-1], arcsec_solution_1[0], arcsec_solution_1[-1]], cmap=cmap)
	line12.set_xticks([])
	line12.text(16640, -4.5, str(posang11), ha='left', color='black', bbox=dict(facecolor='white', alpha=0.8))
	line12.tick_params(axis='both', labelsize=size_of_font)
	line12.set_ylim(arcsec_solution_1[idx_arcsec_lim_lower_1],arcsec_solution_1[idx_arcsec_lim_upper_1])
	line11.errorbar((y1 )/(0.5*y1.max()), x1, xerr=0., color='black', lw=1, alpha=1.0, linestyle='steps-mid')
	line11.set_yticks([])
	line11.set_xticks([])
	line11.errorbar((y_em_1 )/y_em_1.max(), x_em_1, xerr=0., color='red', lw=1, alpha=1.0, linestyle='steps-mid')
	line11.set_yticks([])
	line11.set_xticks([])
	line12.axhline(float(center_qso_left_1), color='black', linestyle='dashed', lw=1)
	line11.axhline(float(center_qso_left_1), color='black', linestyle='dashed', lw=1)
	line12.axhline(float(center_em_1), color='red', linestyle='dashed', lw=1)
	line11.axhline(float(center_em_1), color='red', linestyle='dashed', lw=1)
	line11.tick_params(axis='both', labelsize=size_of_font)

	wavelength_solution_2_rest = wavelength_solution_1/(1+emitter_redshift)
	wavelength_solution_1_rest = wavelength_solution_1/(1+emitter_redshift)

	upward_limit_plotting_1d_trace = np.searchsorted(arcsec_solution_1, (center_em_1 + 1.5*fwhm_em_1))
	downward_limit_plotting_1d_trace = np.searchsorted(arcsec_solution_1, (center_em_1 - 1.5*fwhm_em_1))

	data_spectral_1_plot = np.zeros([data1.shape[1]])
	data_spectral_err_1_plot = np.zeros([data1.shape[1]])
	for i in range(data1.shape[1]):
		data_spectral_1_plot[i] = np.sum(data1[downward_limit_plotting_1d_trace:upward_limit_plotting_1d_trace, i])
		data_spectral_err_1_plot[i] = np.sqrt(np.sum(data_err1[downward_limit_plotting_1d_trace:upward_limit_plotting_1d_trace, i]**2))


	center_of_emission_real = (1+emitter_redshift)*float(line_list[0][1])
	xlimit_absorption_lower = vel_prof(wavelength_solution_1[0], center_of_emission_real)
	xlimit_absorption_upper = vel_prof(wavelength_solution_1[-1], center_of_emission_real)
	velocity_array_emission = vel_prof(wavelength_solution_1, center_of_emission_real)

	line52.plot(velocity_array_emission, data_spectral_1_plot)
	line52.set_xlim(xlimit_absorption_lower, xlimit_absorption_upper)
	line52.set_ylim(-0.3, 0.7)
	line52.axhline(0, color='green', linestyle='dashed')
	line52.axvline(((float(mean_emm)/(1+emitter_redshift))), color='green', linestyle='dashed')
	line52.plot(velocity_array_emission, gaus(wavelength_solution_1, *popt_emm), 'r-')
	line52.tick_params(axis='both', labelsize=size_of_font)
	line12.set_xlabel(r'Observed Wavelength ($\lambda$)', fontsize=1.0*size_of_font)
	line52.set_xlabel(r'Relative Velocity (km s$\rm ^{-1}$) at z=%2.3f' % (emitter_redshift), fontsize=1.0*size_of_font)
	line52.set_ylabel(r'NF', fontsize=1.0*size_of_font)


	################################################PLOT_SPECTRAL_SIDE################################################






	line12.yaxis.set_major_locator(ticker.MultipleLocator(1))
	line12.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
	line12.xaxis.set_major_locator(ticker.MultipleLocator(500))
	line12.xaxis.set_minor_locator(ticker.MultipleLocator(100))
	line12.tick_params(axis = 'both', which = 'major', direction='in', length=size_of_font, width=2, colors='k')	
	line12.tick_params(axis = 'both', which = 'minor', direction='in', length=size_of_font/2, width=1, colors='k')
	line12.axes.get_xaxis().set_visible(True)


	left, width = .25, .5
	bottom, height = .25, .5
	right = left + width
	top = bottom + height

	fig88.text(0.40, 0.85, r'J0025+1145', zorder=13)
	fig88.text(0.08, 0.60, r'Distance along slit (arcsec)', va='center', rotation='vertical', zorder=12)
	plt.savefig('figure_emission_1.pdf')
	print ('Success!')
	#plt.show()
   
button15.on_clicked(analyse_nir_ax)








#Save the parameters
save_data_ax = plt.axes([0.10, 0.03, 0.1, 0.02])
button6 = Button(save_data_ax, 'Save Params')
def save_data_ax(event):
	global clip_val, normal_val, smooth_sigma_val, array_new
	clip_val = sclipping.val
	normal_val = snorm.val
	smooth_sigma_val = sdata_smooth.val
	array_new = []
	array_new = np.append(array_new, [float(clip_val), float(normal_val)])
	array_new = np.append(array_new, [float(smooth_sigma_val)])
	np.savetxt(file_emission_data, array_new)
	print ("Data Saved")
button6.on_clicked(save_data_ax)



#Load saved parameters
get_data_ax = plt.axes([0.20, 0.03, 0.1, 0.02])
button7 = Button(get_data_ax, 'Get Params')
def get_data_ax(event):

	global emission_x_min, emission_x_max, emission_y_min, emission_y_max, region_x_min, region_x_max, region_y_min, region_y_max
	#global array_new in globals
	saved_array = np.loadtxt(file_emission_data)	
	clip_val, normal_val, smooth_sigma_val = saved_array
		
	

	sclipping.set_val(float(clip_val))
	snorm.set_val(float(normal_val))
	sdata_smooth.set_val(float(smooth_sigma_val))


	clip = sclipping.val
	normal = snorm.val
	smooth_sigma = sdata_smooth.val

	

	data_new1 = data1 / 10**(normal)
	data_new1 = np.clip(data_new1, -clip, +clip)
	sigma = [smooth_sigma, smooth_sigma]
	data_new1 = sp.ndimage.filters.gaussian_filter(data_new1, sigma, mode='constant')
	line1.set_data(data_new1)

	fig.canvas.draw_idle()


button7.on_clicked(get_data_ax)





plt.show()


