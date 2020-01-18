

from scipy.optimize import curve_fit
import numpy as np
c = 299792.458 # Speed in Light in Km/s

######################GAUSSIAN RELATED FUNCTIONS##################################

def error_on_gaussian_fit(x_em_4, amp_em_4, amp_err_em_4, mean_em_4, mean_err_em_4, std_em_4, std_err_em_4):
	amp_em_4_for_error = np.random.normal(amp_em_4, amp_err_em_4, 1000) 
	std_em_4_for_error = np.random.normal(std_em_4, std_err_em_4, 1000) 
	mean_em_4_for_error = np.random.normal(mean_em_4, mean_err_em_4, 1000) 
	fitted_area_em_4_err_array = np.array([])
	for amp_em_4_for_error2, std_em_4_for_error2, mean_em_4_for_error2 in zip(amp_em_4_for_error, std_em_4_for_error, mean_em_4_for_error):
		profile_em_4 = np.sum(gaus(x_em_4, amp_em_4_for_error2, mean_em_4_for_error2, std_em_4_for_error2)) * (x_em_4[1]-x_em_4[0])
		fitted_area_em_4_err_array = np.append(fitted_area_em_4_err_array, profile_em_4)

	print (np.mean(fitted_area_em_4_err_array))
	print (np.std(fitted_area_em_4_err_array))

	return((np.mean(fitted_area_em_4_err_array)), (np.std(fitted_area_em_4_err_array)))




def gaus(x,a,x0,sigma):
	return a*np.exp(-(x-x0)**2/(2*sigma**2))


def gaus_fitting(x, y, err):
	amp = y.max()
	sigma = len(x)/4.
	mu = x.max() - ((x.max() - x.min())/2.)
	popt,pcov = curve_fit(gaus,x,y, sigma=err, p0=[5,mu,sigma])
	perr = np.sqrt(np.diag(pcov))	
	return (popt, perr)



def gaus_profiling(x, popt, perr, random_factor):
	amp_dist = np.random.normal(popt[0], 3*perr[0], random_factor)
	mean_dist = np.random.normal(popt[1], 3*perr[1], random_factor)
	std_dist = np.random.normal(popt[2], 3*perr[2], random_factor)
	actual_dist = gaus(x, *popt)	
	test_dist = np.array([])

	for i in range(len(amp_dist)):
		for j in range(len(mean_dist)):
			for k in range(len(std_dist)):
				for l in range(len(x)):
					test_dist = np.append(test_dist, gaus(x[l],amp_dist[i],mean_dist[j],std_dist[k]))

	test_dist = np.reshape(test_dist, (random_factor**3, len(x)))
	test_dist_lower = np.zeros([len(x)])
	test_dist_upper = np.zeros([len(x)])

	for i in range(len(x)):
		test_dist_lower[i] =  np.mean(test_dist[:,i]) - 3*(np.std(test_dist[:,i]))
		test_dist_upper[i] =  np.mean(test_dist[:,i]) + 3*(np.std(test_dist[:,i]))
	
	return (actual_dist, test_dist_lower, test_dist_upper)


######################GAUSSIAN RELATED FUNCTIONS##################################




	



######################ARRAY_SMOOTHING_FUNCTION##################################

def smooth(y, box_pts):
	box = np.ones(box_pts)/box_pts
	y_smooth = np.convolve(y, box, mode='same')
	return y_smooth

######################ARRAY_SMOOTHING_FUNCTION##################################




######################VELOCITY PROFILE##################################

def vel_prof(x, centre):
	xnew = c * ((x-centre)/x)	
	return (xnew)

######################VELOCITY PROFILE##################################



######################AIR TO VACUUM WAVELENGTH CONVERTER##################################

def air2vac(air):
    """ Air to Vaccuum conversion from D. Morton 1991, ApJS 77, 119 """
    if type(air) == float or type(air) == int:
        air = np.array(air)
    air = np.array(air)
    ij = (np.array(air) >= 2000)
    out = np.array(air).copy()
    sigma2 = (1.e4/air)**2
    fact = 1.0 + 6.4328e-5 + 2.94981e-2/(146.0 - sigma2) + 2.5540e-4/(41.0 - sigma2)
    out[ij] = air[ij]*fact[ij]
    return out

######################AIR TO VACUUM WAVELENGTH CONVERTER##################################













