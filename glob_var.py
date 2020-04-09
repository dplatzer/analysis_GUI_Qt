from scipy import constants
import numpy as np

# DATABANK #############
HEV = constants.physical_constants['Planck constant in eV s'][0]  # Planck constant in eV
ME = constants.electron_mass
QE = constants.elementary_charge
C = constants.speed_of_light

IP_HE = 24.58739  # Helium
IP_NE = 21.56454  # Neon
IP_AR = 15.75961  # Argon
IP_KR = 13.99961  # Krypton
IP_XE = 12.12984  # Xenon
IP_N2_X = 15.58  # N2 X

SOFTWARE_LIST = ['Labview 2013', 'Labview 2016', 'PyMoDAQ']
GASLIST = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'N2_X']
IPLIST = [IP_HE, IP_NE, IP_AR, IP_KR, IP_XE, IP_N2_X]
FIRST_HARMLIST = [17, 15, 13, 13, 13, 13]  # to update

RFANO_HE = 60.15  # Helium
RFANO_AR = 26.6  # Argon

# GLOBAL VARIABLES
''' Global means here the variables shared by the different objects in the program'''

path = ""

TOF_resolution = 1e-9  # 1e-9 on SE1, 5e-11 on SE10

software_version = 2  # PyMoDAQ
minus_sign = False  # does the spectrum need to be flipped vertically?
skiplines = 0  # are there any header/do we want to skip some data in the tof spectra?
filenametype = 'Not Needed'
scan_i = 0  # used only in PyMoDAQ to explore the scans in the .h5 file
str_scan_i = '000'
spectrum_i = 0  # used only in PyMoDAQ to explore the spectra in the .h5 file

# background substraction
rmbg_avg_min = 300
rmbg_avg_max = 400

# energy calibration
afit = 0.0
t0fit = 0.0
cfit = 0.0

# experimental parameters
cur_gas_index = 2  # Argon
cur_Ip = IPLIST[cur_gas_index]
cur_SBi = 0
first_harm = FIRST_HARMLIST[cur_gas_index]
lambda_start = 800
cur_nu = C / (lambda_start * 1e-9)
cur_Vp = 0  # retarding potential
cur_L = 2  # length of the TOF

# VMI
vrep = 4.8  # in kV

# energy conversion
# elow = (first_harm - 1) * HEV * cur_nu
elow = cur_Ip + 0.01
ehigh = float(35 * HEV * cur_nu)
dE = 0.01

# scan steps
scanstep_nm = 25
scanstep_fs = scanstep_nm*2/(C*1e-6)
stepsnb = 0

delay_vect = np.zeros([1, 1])
tof_vect = np.zeros([1, 1])
energy_vect = np.zeros([1, 1])
rabbit_mat = np.zeros([1, 1])
rabbitxuvsub_mat = np.zeros([1, 1])
rabbit_TOF_mat = np.zeros([1, 1])
xuvonly_vect = np.zeros([1, 1])
bandsnb = 5
bands_vect = np.zeros([bandsnb, 2])
bands_indices = np.zeros([bandsnb, 2])
SBi = 0

# FT and 2 omega parameters
FT_N = 0  # either Nsteps or npad
FT_Nsteps = 0
FT_padding = True
FT_window = True
FT_zero_order = False
FT_npad = 2048
phase_corr = False
two_w_average = False
two_w_bfilter = False
two_w_integral = False
two_w_phioffset = 0.0

# after FT
freqnorm = np.zeros([1, 1])
FT_ampl = np.zeros([1, 1])
FT_ampl_TOF = np.zeros([1, 1])
FT_phase = np.zeros([1, 1])
FT_phase_TOF = np.zeros([1, 1])
fpeak = np.zeros([1, 1])
fpeak_TOF = np.zeros([1, 1])
peak = np.zeros([1, 1])
peak_TOF = np.zeros([1, 1])
peak_phase = np.zeros([1, 1])
peak_phase_TOF = np.zeros([1, 1])
SNR = np.zeros([1, 1])
phase_err = np.zeros([1, 1])
phase_err_TOF = np.zeros([1, 1])
fit_ampl = np.zeros([1, 1])
fit_phase = np.zeros([1, 1])
fit_phase_err = np.zeros([1, 1])
fit_freq = np.zeros([1, 1])
fpeak_main = 0
fpeak_fit = 0  # used for the cos fit
fpeak_main_TOF = 0
fpeak_lock = 0
rabbitmode = "normal"
contrast = np.zeros([1, 3])

chirp_as = 0.0
chirp_as_TOF = 0.0
chirp_as_eV = 0.0
chirp_as_TOF_eV = 0.0

# booleans
calibloaded = False
rabnormalized = False
rabsmoothed = False
xuvsubstracted = False

''' The list of the variables displayed in the "variable explorer" on the left of the main window '''
varlist = [['*(-1)', minus_sign], ['path', path], ['scan', str_scan_i], ['skip lines', skiplines],
			['Ip', cur_Ip], ['nu', cur_nu],
			['1st harm', first_harm],
			['stepsnb', stepsnb], ['delay', delay_vect], ['energy', energy_vect],
			['XUV', xuvonly_vect], ['rabbit', rabbit_mat], ['rabbitXUVsub', rabbitxuvsub_mat],
			['bands', bands_vect],['bandsindices', bands_indices], ['normalized', rabnormalized], ['smoothed', rabsmoothed], ['subXUV', xuvsubstracted],
			['rabbitmode', rabbitmode], ['FT_ampl', FT_ampl], ['FT_phase', FT_phase],
			['freqnorm', freqnorm], ['fpeak', fpeak], ['peak', peak],
			['peak_phase', peak_phase], ['phase_err', phase_err],
			['SNR', SNR], ['fit_phase', fit_phase],
			['fit_phase_err', fit_phase_err], ['fit_freq', fit_freq],
			['fpeak_main', fpeak_main],
			['fpeak_lock', fpeak_lock],
			['chirp_as', chirp_as], ['chirp_as_eV', chirp_as_eV], ['contrast', contrast]]

def update_varlist():
	global varlist
	varlist = [['*(-1)',minus_sign], ['path',path], ['scan', str_scan_i], ['skip lines', skiplines],
			   ['Ip', cur_Ip],
			   ['nu', cur_nu],
			   ['1st harm', first_harm],
			   ['stepsnb', stepsnb], ['delay', delay_vect],
			   ['energy', energy_vect], ['XUV', xuvonly_vect], ['rabbit', rabbit_mat], ['rabbitXUVsub', rabbitxuvsub_mat],
			   ['bands', bands_vect],['bandsindices', bands_indices], ['normalized', rabnormalized], ['smoothed', rabsmoothed],
			   ['subXUV', xuvsubstracted], ['rabbitmode', rabbitmode], ['FT_ampl', FT_ampl], ['FT_phase', FT_phase],
			   ['freqnorm', freqnorm],['fpeak', fpeak], ['peak', peak],
			   ['peak_phase', peak_phase], ['phase_err', phase_err],
			   ['SNR', SNR],
			   ['fit_phase', fit_phase], ['fit_phase_err', fit_phase_err],
			   ['fit_freq', fit_freq],
			   ['fpeak_main', fpeak_main],
			   ['fpeak_lock', fpeak_lock],
			   ['chirp_as', chirp_as], ['chirp_as_eV', chirp_as_eV], ['contrast', contrast]]

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
	bands_indices = np.zeros([1, 1])

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