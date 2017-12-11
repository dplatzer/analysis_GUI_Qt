from scipy import constants
import numpy as np

###### DATABANK #############
HEV = constants.physical_constants['Planck constant in eV s'][0] #Planck constant in eV
ME = constants.electron_mass
QE = constants.elementary_charge
C = constants.speed_of_light

IP_HE = 24.58739 #Helium
IP_NE = 21.56454 #Neon
IP_AR = 15.75961  #Argon
IP_KR = 13.99961 #Krypton
IP_XE = 12.12984 #Xenon
IP_N2_X = 15.58 #Nitrogen X

GASLIST = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'N2_X']
IPLIST = [IP_HE, IP_NE, IP_AR, IP_KR, IP_XE, IP_N2_X]

RFANO_HE = 60.15 # Helium
RFANO_AR = 26.6 # Argon

###### GLOBAL VARIABLES ###########
''' Global means here the variables shared by the different objects in the program'''
#energy calibration
afit = 0.0
t0fit = 0.0
cfit = 0.0

#experimental parameters
cur_Ip = 0.0
cur_nu = 0.0
cur_Vp = 0.0
cur_L = 0
cur_SBi = 0
first_harm = 0
lambda_start = 800
cur_nu = C / (lambda_start * 1e-9)
cur_Vp = 0 # retarding potential
cur_L = 2 # length of the TOF
first_harm = 13

#energy conversion
elow = (first_harm - 1) * HEV * cur_nu
ehigh = float(52 * HEV * cur_nu)
dE = 0.01

#scan steps
scanstep_nm = 25
scanstep_fs = scanstep_nm*2/(C*1e-6)
stepsnb = 0

delay_vect = np.zeros([1,1])
energy_vect = np.zeros([1,1])
rabbit_mat = np.zeros([1,1])
rabbitxuvsub_mat = np.zeros([1,1])
xuvonly_vect = np.zeros([1,1])
bandsnb = 0
bands_vect = np.zeros([1,1])

# after FT
FT_ampl = np.zeros([1,1])
FT_phase = np.zeros([1,1])
fpeak = np.zeros([1,1])
peak = np.zeros([1,1])
peak_phase = np.zeros([1,1])
fpeak_main = 0
rabbitmode = "normal"
contrast = np.zeros([1,1])

### FT and 2 omega parameters ####
FT_padding = False
FT_window = False
FT_zero_order = True
FT_npad = 4096
two_w_average = False
two_w_bfilter = False
two_w_integral = False
two_w_phioffset = 0.0

chirp_as = 0.0
chirp_as_eV = 0.0

# booleans
rabnormalized = False
rabsmoothed = False
xuvsubstracted = False

''' The list of the variables displayed in the "variable explorer" on the left of the main window '''
varlist = [['afit',afit], ['t0fit',t0fit], ['cfit',cfit], ['Ip',cur_Ip], ['nu',cur_nu], ['Vp',cur_Vp], ['L',cur_L],
		   ['1st harm', first_harm], ['elow', elow], ['ehigh', ehigh], ['dE', dE], ['steps_nm', scanstep_nm],
		   ['steps_fs', scanstep_fs], ['stepsnb', stepsnb], ['delay', delay_vect], ['energy', energy_vect],
		   ['XUV', xuvonly_vect], ['rabbit', rabbit_mat], ['rabbitXUVsub', rabbitxuvsub_mat], ['bandsnb', bandsnb],
		   ['bands', bands_vect], ['normalized', rabnormalized], ['smoothed', rabsmoothed], ['subXUV', xuvsubstracted],
		   ['rabbitmode', rabbitmode], ['FT_ampl', FT_ampl], ['FT_phase', FT_phase], ['fpeak', fpeak], ['peak', peak],
		   ['peak_phase', peak_phase], ['fpeak_main', fpeak_main], ['chirp_as', chirp_as], ['chirp_as_eV', chirp_as_eV],
		   ['contrast', contrast]]

def update_varlist():
	global varlist
	varlist = [['afit', afit], ['t0fit', t0fit], ['cfit', cfit], ['Ip', cur_Ip], ['nu', cur_nu], ['Vp', cur_Vp],
			   ['L', cur_L], ['1st harm', first_harm], ['elow', elow], ['ehigh', ehigh], ['dE', dE],
			   ['steps_nm', scanstep_nm],['steps_fs', scanstep_fs], ['stepsnb', stepsnb], ['delay', delay_vect],
			   ['energy', energy_vect], ['XUV', xuvonly_vect], ['rabbit', rabbit_mat], ['rabbitXUVsub', rabbitxuvsub_mat],
			   ['bandsnb', bandsnb], ['bands', bands_vect], ['normalized', rabnormalized], ['smoothed', rabsmoothed],
			   ['subXUV', xuvsubstracted], ['rabbitmode', rabbitmode], ['FT_ampl', FT_ampl], ['FT_phase', FT_phase],
			   ['fpeak', fpeak], ['peak', peak], ['peak_phase', peak_phase], ['fpeak_main', fpeak_main], ['chirp_as', chirp_as],
			   ['chirp_as_eV', chirp_as_eV], ['contrast', contrast]]

def clear_varlist():
	global afit, t0fit, cfit, delay, energy_vect, rabbit_mat, xuvonly_vect, bandsnb, bands_vect
	afit = 0.0
	t0fit = 0.0
	cfit = 0.0
	delay = np.zeros([1, 1])
	energy_vect = np.zeros([1, 1])
	rabbit_mat = np.zeros([1, 1])
	xuvonly_vect = np.zeros([1, 1])
	bandsnb = 5
	bands_vect = np.zeros([1, 1])

def print_cts():
	print("Ip = " + str(cur_Ip) + " eV")
	print("nu = %.2e Hz" % cur_nu)
	print("Vp (retarding pot.) = " + str(cur_Vp) + " V")
	print("L (TOF length) = " + str(cur_L) + " m")
	print("First harmonic: " + str(first_harm))

	#just a synthax memo for exception handling:
	#try:
	#	...
	#except Exception:
	#	print(traceback.format_exception(*sys.exc_info()))