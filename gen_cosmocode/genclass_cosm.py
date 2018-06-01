import numpy as np 
import matplotlib.pyplot as plt 

from scipy.integrate import quad
from astropy.cosmology import FlatLambdaCDM
from cmb_prior import z_lss 
from precalc_func import mk_icov_cmb_withOb_al as cmb_icov
 
class ezs_for_cosmos:

	"""
	A storage class for different E(z) expressions for the cosmologies
	E(z) = sqrt(H^2/H_0^2)

	"""
        def ez_sr(self, z,theta):
    		om, w0, h0 = theta

    		#equation 20; Planck mod-grav paper
    		mat_term = om * (1+z)**3
    		ode = 1 - om
    		rad_term = pow(h0/100., -2.) * 2.469e-5 * (1+0.2271*3.04) * (1+z)**4.
    		de_term  = ode*pow((1+z)**3/(om*(1+z)**3 + ode) ,w0/ode)    
    		ezsq  = mat_term + de_term + rad_term
    		return 1./np.sqrt(ezsq)

	def ez_hlg(self, z, theta):
		om, beta = theta 
		mat_term = om*(1+z)**3. 	#while defining this for each cosmology is "annoying" I'd rather continue like this for safety
		rad_term = pow(h0/100., -2.) * 2.469e-5 * (1+0.2271*3.04) * (1+z)**4.
		alpha = np.exp(1 - om)		
		mod_term = np.log(alpha + beta*z)
		ezsq = mat_term + de_term + rad_term
		return 1./np.sqrt(ezsq)

	def ez_pngb(self, z, theta):
		om, w0, F = theta
		mat_term = om*(1+z)**3. 	#while defining this for each cosmology is "annoying" I'd rather continue like this for safety
		rad_term = pow(h0/100., -2.) * 2.469e-5 * (1+0.2271*3.04) * (1+z)**4.
		ode = 1 - om
		def eos_intterm(z, w0, F):
		    eos = -1 + (1+w0)*pow(1+z, -F)
		    return (1 + eos)/(1+z)
		de_term = ode*np.exp(3*quad(eos_intterm, 0, z, args=(w0, F))	
		ezsq = mat_term + de_term + rad_term	#same as the case for the individual terms
		return 1./np.sqrt(ezsq)

	def ez_qcdghost(self, z, theta):
		"""
		QCD Ghost Free Dark Energy: Zhang et al. find that this model is 
		"""
		om, gamma = theta
		mat_term = om*(1+z)**3. 	#while defining this for each cosmology is "annoying" I'd rather continue like this for safety
		omega_rad = pow(h0/100., -2.) * 2.469e-5 * (1+0.2271*3.04) 		
		rad_term =  omega_rad * (1+z)**4.
		kappa = (1 - (om + omega_rad)/gamma)/2.
		ezsq = kappa + np.sqrt(kappa**2+(mat_term+rad_term)/gamma)
		return 1./np.sqrt(ezsq)


class distance_for_cosmo:

	def __init__(self, ez, theta):
		"""
		Initiate the class with the expression  E(z) = sqrt(H^2/H_0^2) for the given cosmology 
		(see the ezs_for_cosmos class where these are catalogued)

		The second variable is the set of the parameters for which the cosmological observables should be 
		calculated.
		Note: the format should be (Omega_m, xx, yy, h0) such that the first parameter is om and the last is h0


		Definitions: 
		1. c
		2. Obh2
		3. Omega_gamma
		"""

		self.ez = ez
		self.theta = theta
		self.c = 299792.458
		self.h0 = theta[-1]
		self.om = theta[0]
		self.obh2 = 0.02262
		self.omega_gam = 2.469e-5
		self.zstar = z_lss(theta[0], theta[-1], self.obh2)

	def luminosity_distance(self,z, log=True):
		"""
		Use the expression for the luminosity distance from Carroll et al. 1992 (or any cosmology book)
		d_L = c*(1+z)*int(0, z) dz/(E(z)*H0)
		(assume flatness)
		returns if log=True, distance modulus, else distance in mega parsec

		"""
		out = np.zeros_like(z)
		npar = len(self.theta)

		for i, zval in enumerate(z):
			out[i] = quad(self.ez, 0, zval, args=(self.theta))[0]	
		dl_mpc = c*out*(1+z)/self.h0	

		if log:
			mu = 5*np.log10(dl_mpc)+25	
			return mu
		else:
			return dl_mpc

	def angular_diameter_distance(self, z):
		"""
		D_A = D_L/(1+z)**2.
		return D_A in Mpc
		"""
		dl_mpc = self.luminosity_distance(z, log=False)
		return dl_mpc/(1+z)**2.

	def cmb_shift(self):
		"""
		CMB shift: R = sqrt(Om * h0^2)*Dc(z*)/c

		"""	
		om = self.om
		ezterm = quad(ez, 0., zstar, args=(self.theta))[0]
		return np.sqrt(om)*ezterm

	def rs(self, z):
		"""
		The size of the sound horizon at a redshift
		"""
    		a = 1./(1+z)
    		#define E(z)
    		ha = self.h0*(1./self.ez(z, self.theta))**2.
    		denom = a**2* ha * np.sqrt((1 + (3*a*self.obh2/(4*self.omega_gam))))
    		return 1./denom
	
	def cmb_la(self):
		"""
		Position of the first acoustic peak in the power spectrum
		la = pi*DC/rs --> see Planck 2015
		"""
		zstar = z_lss(self.om, self.h0, self.obh2)
		dc_decoupling = angular_diameter_distance(zstar)*(1+zstar)
		rs_decoup = quad(self.rs, 0., zstar, args=(self.theta))
		return np.pi*dc_decoupling/rs_decoup



class likelihoods_cosmology:
	"""
	Defined the likelihoods for the individual cosmological probes.
	For initialising, it requires the E(z) expression, the value of theta for evaluating and the fiducial cosmology
	"""
	def __init__(self, ez, theta, theta_fid, datasets):

		self.ez = ez
		self.theta = theta
		self.theta_fid = theta_fid

		self.func_fid = distance_for_cosmo(ez, theta_fid)

		#self.obh2 = 0.02262
		#define the CMB fiducial observations
		self.cmb_obs = np.array([func_fid.cmb_shift(), func_fid.cmb_la(), func_fid.obh2])
		self.icov_mat = cmb_icov()

		#load the SN covariance matrix for the DESIRE survey	
		self.zdist = np.arange(0, 1.56, 0.025)
		self.z_ast = zdist[2:]	
		self.sn_obs = func_fid.luminosity_distance(z_ast)*self.h0/(c*(1+z_ast))		
		self.sn_cov_mat = np.loadtxt('../anc_data/piecewise_weight.ascii', skiprows=1)
				
		#H(z) measurements
		self.measure_types = ['chronometers', 'radBAO']
		#self.hz_obs = func_fid.
		
		#local H0
		self.h0val = 73.24
		self.efrac = 0.015		
		
		#which probes to use
		self.addcmb, self.addsn, self.addlocal_h0 = datasets
		
	"""
	To-dos: 	
		1. A dictionary with input "datasets" as an attribute for the class
	"""
	def cmb_likelihood(self):
		"""
		The chisquare for the CMB include the R, la and Obh2 values along with the covariance matrix from the Planck 2015 paper
		"""
		func_init = distance_for_cosmo(ez, theta)
		#define the "theoretical" R and lA for the given cosmology
		rth = func_init.cmb_shift()
		lath = func_init.cmb_la()

		#make the theoretical array
		th_array = np.array(rth, lath, func_init.obh2)
		dif_arr = th_array - self.cmb_obs
		cmb_chisq = np.dot(dif_arr.T, np.dot(self.icov_mat, dif_arr))
		return cmb_chisq
	
	def sn_likelihood(self):
			
		func_init = distance_for_cosmo(ez, theta)
		sn_dist = func_init.luminosity_distance(z_ast)*self.h0/(c*(1+z_ast))

		dif_arr = sn_dist - sn_obs
		sn_chisq = np.dot(dif_arr.T, np.dot(self.sn_cov_mat, dif_arr))
		return sn_chisq

	def local_h0_likelihood(self):
		"""
		Likelihood for the local H0 measurements. Currently just a gaussian centered on the value defined 
		and fractional error
		Possibilities:
		1. Actually derive a likelihood from e.g. Feeney et al. 2018
		"""
		h0_chisq = (self.h0 - self.h0val)**2./((self.h0val*self.efrac)**2)		
		return h0_chisq
	
	def total_llhood(self):
		"""
		Return the log likelihood for the combination of probes requested
		
		"""	
		tchisq = 0.
		if addcmb:
			cmb = cmb_likelihood()
			tchisq += cmb

		if addsn:		
			sn = sn_likelihood()
			tchisq += sn

		if addlocal_h0:		
			local_h0 = local_h0_likelihood()
			tchisq += h0

		return -0.5*chisq

