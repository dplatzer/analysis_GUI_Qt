# _Rabbit_functions.py

from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtCore import Qt
from glob import glob
import traceback, sys, os
import numpy as np
import scipy.optimize as opt
import scipy
import h5py

import glob_var as cts
import analysis_functions as af


class Rabbit_functions_mixin:
    """ This is not really a standalone class. It's a trick I use to split the RabbitWin class definition
    ("Rabbit-win.py") into two files: one for the GUI, one for the function; but keeping the advantage
    of sharing instance variables (self.something) between the two

    Here are defined the majority of the functions (methods to be correct) used in the RabbitWin object,
    especially analysis functions.
    """

    def importdata_lr(self):
        ''' "[Import] data" listener. Loads all the tof files, one by one, converts them to energy and regroup them into the
            cts.rabbit_mat array'''
        data_f = QFileDialog.getOpenFileName(self, 'Import data')
        data_filename = data_f[0]
        data_taken = False
        if (data_filename):
            if cts.software_version != 2:
                fdir, fname = os.path.split(data_filename)
                fdir = fdir + '/'
                #flist = glob(fdir + '*[0-9][0-9][0-9][0-9].txt')
                flist = glob(fdir + cts.filenametype)
                cts.path = fdir

                if (len(flist) != 0):
                    data_taken = True
                    cts.stepsnb = len(flist)
                    cts.delay_vect = np.linspace(0, cts.scanstep_fs * cts.stepsnb * 1e-15, cts.stepsnb, endpoint=False)

                    data = np.transpose(self.import_tof(data_filename))
                    cts.tof_vect = data[0, :]
                    cts.tof_vectlength = data.shape[1]

                    i = 0
                    for fn in flist:
                        self.window().statusBar().showMessage("Processing " + str(i))
                        self.window().statusBar().repaint()

                        data = np.transpose(self.import_tof(fn))
                        if cts.minus_sign:
                            data[1, :] = (-1) * data[1, :]
                        counts = data[1, :] - \
                                 data[1, cts.rmbg_avg_min:cts.rmbg_avg_max].mean()
                        cts.energy_vect, counts2 = af.jacobian_transform(cts.tof_vect[::-1], counts)

                        self.elength = cts.energy_vect.shape[0]

                        if i == 0:
                            cts.rabbit_mat = np.zeros([cts.stepsnb, cts.energy_vect.shape[0]])
                            self.data_tof = []
                            #print(cts.rabbit_mat.shape)
                        cts.rabbit_mat[i, :] = counts2
                        i += 1

                    cts.rabbit_mat = np.array(cts.rabbit_mat)
                    # cts.rabbit_mat = cts.rabbit_mat[::-1, :]
                    # self.sig_int_custom()
                else:
                    self.window().statusBar().showMessage("Incorrect directory")

            # PyMoDAQ
            else:
                cts.path = data_filename
                data_taken = True

                with h5py.File(data_filename, 'r') as file:
                    data_path = ''.join(['Scan', cts.str_scan_i, '/Detector000/Data1D/Ch000/'])
                    delay_path = ''.join(['Scan', cts.str_scan_i, '/Scan_x_axis'])

                    raw_datas = file['Raw_datas']
                    data_scan = np.array(raw_datas[data_path + 'Data'])[:,cts.skiplines:]  # whole rabbit trace
                    cts.tof_vect = np.array(raw_datas[data_path + 'X_axis'])[cts.skiplines:]

                    cts.delay_vect = np.array(raw_datas[delay_path])

                cts.scanstep_nm = (cts.delay_vect[1:] - cts.delay_vect[:-1]).mean() * 1000
                cts.delay_vect = (cts.delay_vect - cts.delay_vect[0])*2/(cts.C*1e-9)
                cts.stepsnb = data_scan.shape[0]
                self.scanparam_le.setText("{:.2f}".format(cts.scanstep_nm))

                for i in range(len(cts.delay_vect)):
                    if cts.minus_sign:
                        data_scan[i, :] = (-1) * data_scan[i, :]
                    counts = data_scan[i, :] - \
                             data_scan[i, cts.rmbg_avg_min:cts.rmbg_avg_max].mean()
                    cts.energy_vect, counts2 = af.jacobian_transform(cts.tof_vect[::-1], counts)

                    if i == 0:
                        cts.rabbit_mat = np.zeros([cts.stepsnb, cts.energy_vect.shape[0]])
                        # print(cts.rabbit_mat.shape)
                    cts.rabbit_mat[i, :] = counts2

            if data_taken:
                # Updating the graph
                # self.rab_widget.refreshplot_fn(cts.energy_vect, cts.delay_vect, cts.rabbit_mat)
                self.rab_widget.colorauto_cb.setEnabled(True)
                self.rab_widget.logscale_cb.setEnabled(True)

                self.window().statusBar().showMessage("Data loaded successfully")
                self.dataloaded = True
                self.importXUV_btn.setEnabled(True)

                # Updating the "Signal Treatment" Groupbox
                self.smooth_btn.setEnabled(True)
                self.normalize_btn.setEnabled(True)

                # Updating the "Select bands" Groupbox
                self.showbands_cb.setChecked(True)
                self.color = "black"

                # finding the oscillation frequency using rainbow rabbit analysis
                cts.FT_padding = True
                cts.rabbitmode = 'rainbow'
                self.rabbit_fn()
                # self.rabbit_done = False # so I don't display amp and phase
                self.refreshplot()

                # Updating the "RABBIT" Groupbox
                self.rainbowrab_btn.setEnabled(True)
                self.exportrab_btn.setEnabled(True)
                self.lock_2w_cb.setEnabled(True)

                self.int_rab_cb.setEnabled(True)

                self.window().updateglobvar_fn()

    def importXUV_lr(self):
        ''' "[Import] XUV only" button listener. Loads XUV only from file (it is converted from tof to energy)'''
        self.update_envect_fn()
        xuv_fname = QFileDialog.getOpenFileName(self, "Import XUV only")
        xuv_f = xuv_fname[0]
        if (xuv_f):
            data = []
            try:
                data = np.loadtxt(xuv_f, unpack=True)
                cts.tof_vect = data[0, :]
                cts.tof_vect = cts.tof_vect[::-1]
                cts.xuvonly_vect = data[1, :] - data[1, 300:600].mean()
                cts.energy_vect, cts.xuvonly_vect = af.jacobian_transform(cts.tof_vect, cts.xuvonly_vect)
                cts.xuvonly_vect = np.array(cts.xuvonly_vect)
                self.xuvonlyloaded = True
                if self.bandselected:
                    self.subXUV_btn.setEnabled(True)

                self.window().updateglobvar_fn()

            except ValueError:
                self.window().statusBar().showMessage("Incorrect XUV only file")

    def importRabbitVMI_lr(self):
        """ Used to load VMI matrices and to convert radius to energy"""

        fname = QFileDialog.getOpenFileName(self, "Import RABBIT matrix")
        f = fname[0]
        if f:
            cts.path = f
            try:
                rab = np.loadtxt(f)
            except Exception:
                print(traceback.format_exception(*sys.exc_info()))

            rab = rab.T
            dshape = rab.shape
            cts.stepsnb = dshape[0]

            Vrep = cts.vrep  # in kV
            elow = 0.01 + cts.cur_Ip
            ehigh = 3.078 * Vrep + cts.cur_Ip
            self.elow_le.setText("{:.2f}".format(elow))
            self.ehigh_le.setText("{:.2f}".format(ehigh))
            self.update_envect_fn()

            cts.delay_vect = np.linspace(0, cts.scanstep_nm * 2 / cts.C * cts.stepsnb * 1e-9, cts.stepsnb,
                                         endpoint=False)
            nbsteps_E = (cts.ehigh - cts.elow) / cts.dE + 1
            cts.energy_vect = np.linspace(cts.elow, cts.ehigh, nbsteps_E)

            if self.need_calib_cb.isChecked():
                r = np.linspace(1, 256, 256) # pixel vector
                # calib from Guillaume, I added the +Ip and put Vrep as a parameter
                a = 2.239E-5*14.7872/7.6606 * Vrep
                b = 7.665E-9*14.7872/7.6606 * Vrep
                EeV = a*r**2 + b*r**3 + cts.cur_Ip
                Ederiv_over_r = 2*a + 3*b*r # see remark below

                # from now on I mimic the function af.jacobian_transform except that the jacobian transform
                # is slightly different: as a 1/r correction is already performed during the Abel inversion
                # process, we have to multiply by r*dr/dE (ie divide by 1/r * dE/dr)
                countsnew = np.zeros((dshape[0], dshape[1]))
                for i in range(dshape[0]):
                    countsnew[i, :] = rab[i, :]/Ederiv_over_r
                f = scipy.interpolate.interp1d(EeV, countsnew)
                cts.rabbit_mat = f(cts.energy_vect)
            else:
                cts.rabbit_mat = rab


            # Updating the graph
            # self.rab_widget.refreshplot_fn(cts.energy_vect, cts.delay_vect, cts.rabbit_mat)
            self.rab_widget.colorauto_cb.setEnabled(True)
            self.rab_widget.logscale_cb.setEnabled(True)

            self.window().statusBar().showMessage("Data loaded successfully")
            self.dataloaded = True
            self.importXUV_btn.setEnabled(True)

            # Updating the "Signal Treatment" Groupbox
            self.smooth_btn.setEnabled(True)
            self.normalize_btn.setEnabled(True)

            # Updating the "Select bands" Groupbox
            self.showbands_cb.setChecked(True)
            self.color = "black"

            # finding the oscillation frequency using rainbow rabbit analysis
            cts.FT_padding = True
            cts.rabbitmode = 'rainbow'
            self.rabbit_fn()
            # self.rabbit_done = False # so I don't display amp and phase
            self.refreshplot()

            # Updating the "RABBIT" Groupbox
            self.rainbowrab_btn.setEnabled(True)
            self.exportrab_btn.setEnabled(True)
            self.lock_2w_cb.setEnabled(True)

            self.int_rab_cb.setEnabled(True)

            self.window().updateglobvar_fn()

    def smoothrab_lr(self):
        ''' "smooth RABBIT" button listener'''
        try:
            sm = int(self.smooth_le.text())
            for i in range(cts.stepsnb):
                cts.rabbit_mat[i, :] = af.smooth(cts.rabbit_mat[i, :], sm)
            self.rab_widget.colorauto_cb.setCheckState(Qt.Checked)
            self.rab_widget.logscale_cb.setCheckState(Qt.Unchecked)

            self.smooth_btn.setText("smoothed")
            self.smooth_btn.setEnabled(False)
            self.smooth_le.setEnabled(False)
            cts.rabsmoothed = True
            self.window().updateglobvar_fn()

            self.rab_widget.refreshplot_fn(signal=cts.rabbit_mat)
        except ValueError:
            self.window().statusBar().showMessage("Smooth must be an integer")

    def normalizerab_lr(self):
        ''' "normalize RABBIT" button listener'''
        for i in range(cts.rabbit_mat.shape[0]):
            cts.rabbit_mat[i, :] = cts.rabbit_mat[i, :] / cts.rabbit_mat[i, :].sum()
        # self.rab_widget.colorauto_cb.setCheckState(Qt.Checked)
        # self.rab_widget.logscale_cb.setCheckState(Qt.Unchecked)

        self.normalize_btn.setText("normalized")
        self.normalize_btn.setEnabled(False)
        cts.rabnormalized = True
        self.window().updateglobvar_fn()

        self.rab_widget.refreshplot_fn(signal=cts.rabbit_mat, keep_lims=True)

    def sig_int_custom(self):
        # integrates the signal along custom ranges
        hnu = cts.HEV * cts.cur_nu
        harm = [25, 27]
        ranges = [[(h-0.5)*1.55, (h+0.5)*1.55] for h in harm]
        i=0
        for xi, xf in ranges:
            ji = np.argmin(abs(cts.energy_vect - xi))
            jf = np.argmin(abs(cts.energy_vect - xf)) + 1
            x_data = np.trapz(cts.rabbit_mat[:, ji:jf], cts.energy_vect[ji:jf], axis=1)
            x_data = np.sum(x_data)
            print('HH',harm[i],'range ', ranges[i], 'integral', x_data)
            i+=1

    def normal_rabbit_lr(self):
        cts.rabbitmode = 'normal'
        self.rabbit_fn()
        self.refreshplot()

    def rainbow_rabbit_lr(self):
        cts.rabbitmode = 'rainbow'
        self.rabbit_fn()
        if self.clickcount == 2*cts.bandsnb:
            self.done_btn.setEnabled(True)
        self.refreshplot()

    def rabbit_fn(self):
        ''' ANALYSIS - either for normal or Rainbow RABBIT '''
        self.peak_phase = []
        self.phase_err = []
        self.fpeak = []
        self.ampl = []
        self.ang = []
        self.energy = []
        data = []

        hnu = cts.HEV * cts.cur_nu
        self.SBi = cts.first_harm + 1

        if cts.xuvsubstracted:
            signal = cts.rabbitxuvsub_mat
        else:
            signal = cts.rabbit_mat

        jx = []
        # if we don't select bands (rainbow mode), we analyse the whole spectrogram
        if cts.rabbitmode == 'rainbow':
            ji = np.argmin(cts.energy_vect)
            jf = np.argmax(cts.energy_vect)
            jx.append([ji, jf])
        else:
            for xi, xf in cts.bands_vect:
                ji = np.argmin(abs(cts.energy_vect - xi))
                jf = np.argmin(abs(cts.energy_vect - xf)) + 1
                jx.append([ji, jf])
        jx = np.array(jx)
        cts.bands_indices = jx

        if cts.rabbitmode == 'normal':
            self.peak = np.zeros(cts.bandsnb)
            i = 0

            print('integral along bands:')
            for ji, jf in jx[:]:
                if jf - ji > 1:  # we integrate only if we have more than one point
                    x_data = np.trapz(signal[:, ji:jf], cts.energy_vect[ji:jf], axis=1)
                else:
                    x_data = signal[:, ji]
                print(x_data[0])

                self.freqnorm, ampl, ang = af.FFT(x_data, cts.scanstep_fs*1e-15)
                print(ang)

                # if self.lock_2w_cb.isChecked():
                #     pass
                # else:
                #     fpeak, peak, peak_phase, phase_err = af.find_2w(x_data, cts.scanstep_fs*1e-15)

                data.append(x_data)
                self.ampl.append(np.array(ampl))
                # self.peak[i] = np.absolute(peak)
                # self.peak_phase.append(peak_phase)
                # self.phase_err.append(phase_err)
                # self.fpeak.append(fpeak)
                self.ang.append(ang)

                i += 1

        elif cts.rabbitmode == 'rainbow':
            self.peak = []

            for ji, jf in jx[:]:
                for x in range(ji,jf+1):
                    x_data = signal[:, x]
                    # f2om, peak, peak_phase, phase_err = af.find_2w(x_data, cts.scanstep_fs*1e-15)
                    self.freqnorm, ampl, ang =  af.FFT(x_data, cts.scanstep_fs*1e-15)

                    data.append(x_data)
                    self.ang.append(ang)
                    self.ampl.append(ampl)
                    # self.peak.append(np.absolute(peak))
                    self.energy.append(cts.energy_vect[x])
                    # self.peak_phase.append(peak_phase)
                    # self.phase_err.append(phase_err)
                    # self.fpeak.append(f2om)

            self.energy = np.array(self.energy)


        else:
            print('Incorrect cts.rabbitmode value')
            raise ValueError

        cts.FT_ampl = np.array(self.ampl)
        cts.FT_phase = np.array(self.ang)
        # cts.fpeak = np.array(self.fpeak)
        # cts.peak = np.array(self.peak)
        # cts.peak_phase = (-1) * np.array(self.peak_phase)
        # cts.phase_err = np.array(self.phase_err)

        cts.freqnorm = np.array(self.freqnorm)

        if self.lock_2w_cb.isChecked():
            self.fpeak_index = int(self.set_2w_le.text())
            cts.fpeak_main = cts.freqnorm[self.fpeak_index]
            cts.fpeak_lock = cts.fpeak_main
        else:
            start_i = cts.FT_N//10 # we look for a peak that is not the 0 frequency component
            # so we investigate from here to the end of the spectrum

            # The 2w frequency that we choose (and consider the main one) is the one for which
            # the 2w ampl is the biggest

            # see the argmax definition (and the examples)
            fpeak_main_i = np.unravel_index(np.argmax(cts.FT_ampl[:, start_i:], axis=None),
                                            cts.FT_ampl[:, start_i:].shape)[1]
            fpeak_main_i += start_i
            cts.fpeak_main = cts.freqnorm[fpeak_main_i]
            print('cts.fpeak_main', cts.fpeak_main)
            self.fpeak_index = fpeak_main_i
            self.set_2w_le.setText(str(self.fpeak_index))

        # see Q:\LIDyL\Atto\ATTOLAB\SE1\Analyse_RABBIT\DFT_Phase_in_Rabbit\DFT_Phase_in_Rabbit.pdf
        # works only if nu is determined with 0-padding, and for npad >= 2048 the correction is not
        # necessary
        if cts.phase_corr:
            print('with Phase correction')
            N = cts.FT_N  # either Nsteps or npad
            Nsteps = cts.FT_Nsteps
            print(Nsteps, N)
            dt = cts.scanstep_fs*1e-15
            nu = cts.fpeak_main*cts.cur_nu
            Extraphase_DP = np.zeros(N // 2)
            for j in range(N // 2):
                Extraphase_DP[j] = np.pi * (Nsteps - 1) * (nu * dt - j / N)
                cts.FT_phase[:, j] = cts.FT_phase[:, j] - Extraphase_DP[j]
            # cts.FT_phase = af.wrap2pmpi(np.unwrap(cts.FT_phase, axis=1))



        cts.fpeak = cts.freqnorm[self.fpeak_index]
        cts.peak = cts.FT_ampl[:, self.fpeak_index]
        cts.peak_phase = cts.FT_phase[:, self.fpeak_index]
        cts.peak_phase = np.unwrap(cts.peak_phase)

        # calculating SNR and phase error based on SNR ("Thierry's" technique)
        l = cts.FT_ampl.shape[1]

        # for steps 75nm
        # cts.SNR = cts.FT_ampl[:, self.fpeak_index]/\
        #           np.sqrt((cts.FT_ampl[:,200:600]**2).mean()*400)
        # for steps 25nm or similar
        cts.SNR = cts.FT_ampl[:, self.fpeak_index] / \
                  np.sqrt((cts.FT_ampl[:, l//2:] ** 2).mean() * l)

        cts.phase_err = (0.03421*2048+19.26)/(2048+57.47)/cts.SNR

        fit_ampl = []
        fit_phase = []
        fit_phase_err = []
        fit_freq = []
        for i in range(len(cts.peak)):
            """ Not working well yet (I have to lock also the fit frequency)"""
            # doing a 4 parameters cosine fit to measure the phase
            # the fit params are: amplitude 'a', offset 'b', phase 'phi' and frequency 'f0'
            # (see definition of af.cosine)
            max = data[i].max()
            mean = data[i].mean()
            popt = [0,0,0,0]
            pcov = [1,1,1,1]
            # popt, pcov = opt.curve_fit(af.cosine, cts.delay_vect, data[i],
            #                            p0=[max-mean, mean, 1, cts.fpeak_main*cts.cur_nu])
            if popt[0] < 0:
                popt[0] = -1 * popt[0]
                popt[2] += + np.pi
            popt[2] = -1 * popt[2] # flipping the phase, just like when doing FT analysis
            if i == 0:
                popt[2] = np.unwrap([0, popt[2]])[1]
            fit_ampl.append(popt[0])
            fit_phase.append(popt[2])
            fit_phase_err.append(np.sqrt(np.diag(pcov))[2])
            fit_freq.append(popt[3]/cts.cur_nu)

        cts.fit_ampl = np.array(fit_ampl)
        cts.fit_phase = np.array(fit_phase)
        cts.fit_phase_err = np.array(fit_phase_err)
        # cts.fit_phase = af.wrap2pmpi(cts.fit_phase)
        cts.fit_freq = np.array(fit_freq)

        if cts.rabbitmode == 'normal':
            self.SBi = cts.SBi
            self.sborder = np.linspace(self.SBi, self.SBi + 2 * (cts.bandsnb - 1), cts.bandsnb)
            self.pa = np.polyfit(self.sborder, cts.peak_phase, 1)
            self.attochirp = self.pa[0] * np.append(self.sborder, self.sborder[-1]+2) + self.pa[1]
            cts.chirp_as = self.pa[0] / (
                    cts.cur_nu * 2 * np.pi) * 1e18  # there is a factor of 2 due to the way pa is calculated
            cts.chirp_as_eV = float(cts.chirp_as / (2 * hnu))

            self.FTcontrast_btn.setEnabled(True)

        self.export_2w_btn.setEnabled(True)
        self.export_2w_std_btn.setEnabled(True)
        self.plotSBvsdelay_btn.setEnabled(True)
        self.clear_btn.setEnabled(True)

        self.rabbit_done = True
        self.window().updateglobvar_fn()

    def FTselect_fn(self):
        """Used to select energy bands to analyze using FT. The bands are defined here as
        being where the 2w amplitude is above a certain threshold, defined with a number
        between 0 and 1 in self.FT_threshold_le QLineEdit.
        Example: for a threshold of 0.7, the band will be the region where the 2w amplitude
        is greater than 0.7*(max 2w amplitude) for a given sideband.

        :return: None
        """
        amp_threshold = float(self.FTthreshold_le.text())
        hnu = cts.HEV * cts.cur_nu
        self.SBi = cts.SBi

        if cts.rabbitmode == 'normal':
            cts.rabbitmode = 'rainbow'
            self.rabbit_fn()

        if self.integrate_cb.isChecked():
            for i in range(cts.bandsnb):
                iemin = np.argmin(abs(cts.energy_vect - (self.SBi - 0.6 + 2 * i) * hnu))
                iemax = np.argmin(abs(cts.energy_vect - (self.SBi + 0.6 + 2 * i) * hnu))

                ampl = cts.peak[iemin:iemax]
                dampl = np.zeros(ampl.shape[0])

                for k in range(dampl.shape[0] - 1):  # dampl[-1] = 0
                    dampl[k] = ampl[k + 1] - ampl[k]
                    # print(ampl.max())
                xi = []
                for l in range(dampl.shape[0]):
                    if (ampl[l] > amp_threshold * ampl.max() and dampl[l] > 0):
                        xi.append(l)
                for l in range(dampl.shape[0]):
                    if (ampl[l] > amp_threshold * ampl.max() and dampl[l] < 0):
                        xi.append(l)
                # in this way, xi is probably messy but what import are the first and last values

                cts.bands_vect[i, 0] = cts.energy_vect[xi[0] + iemin]
                cts.bands_vect[i, 1] = cts.energy_vect[xi[-1] + iemin]
        else:
            for i in range(cts.bandsnb):
                iemin = np.argmin(abs(cts.energy_vect - (self.SBi - 0.3 + 2 * i) * hnu))
                iemax = np.argmin(abs(cts.energy_vect - (self.SBi + 0.3 + 2 * i) * hnu))

                ampl = cts.peak[iemin:iemax]
                xi = np.argmax(ampl)
                cts.bands_vect[i, 0] = cts.energy_vect[xi + iemin]
                cts.bands_vect[i, 1] = cts.energy_vect[xi + iemin]

        self.clickcount = 2*cts.bandsnb
        self.clearbands_btn.setEnabled(True)
        # self.done_btn.setEnabled(True)

        self.refreshplot(keep_lims=True)

    def theoryselect_fn(self):
        """Used to select energy bands centered around expected values from the calibration
        and with width defined in self.theorythreshold_le QLineEdit.
        Example: with threshold 0.3, the bands will be 2*0.3*(photon energy) wide

        :return: None
        """
        hnu = cts.HEV * cts.cur_nu
        self.SBi = cts.SBi

        if self.integrate_cb.isChecked():
            for i in range(cts.bandsnb):
                emin = np.argmin(
                    abs(cts.energy_vect - (self.SBi - float(self.theorythreshold_le.text()) + 2 * i) * hnu))
                emax = np.argmin(
                    abs(cts.energy_vect - (self.SBi + float(self.theorythreshold_le.text()) + 2 * i) * hnu))

                cts.bands_vect[i, 0] = cts.energy_vect[emin]
                cts.bands_vect[i, 1] = cts.energy_vect[emax]

        else:
            for i in range(cts.bandsnb):
                ecenter = np.argmin(abs(cts.energy_vect - (self.SBi + 2 * i) * hnu))
                cts.bands_vect[i, 0] = cts.energy_vect[ecenter]
                cts.bands_vect[i, 1] = cts.energy_vect[ecenter]

        self.clickcount = 2 * cts.bandsnb

        self.refreshplot(keep_lims=True)
        self.clearbands_btn.setEnabled(True)

    def refreshplot(self, keep_lims=False):
        """This function plots/replots the three graphs (rabbit trace (self.rab_widget), amplitude
         (self.FT_ax) and phase (self.phase_ax) of the Fourier Transformed trace).

        :param keep_lims: if True, the limits (both in x and y) of self.phase_ax and
                                    self.FT_ax are preserved
        :return: None
        """
        hnu = cts.HEV * cts.cur_nu
        self.rab_widget.refreshplot_fn(cts.energy_vect, cts.delay_vect,
                                       cts.rabbit_mat, keep_lims=keep_lims)

        if self.showbands_cb.isChecked() and not cts.np.array_equal(
                cts.bands_vect,np.zeros([cts.bandsnb, 2])):
            for i in range(2 * cts.bandsnb):
                xval = cts.bands_vect[i // 2, i % 2]
                if xval != 0:
                    self.rab_widget.ax.axvline(xval, color=self.color)
                    self.rab_widget.fc.draw()

        if self.rabbit_done:
            xmin_ph, xmax_ph = self.phase_ax.get_xlim()
            ymin_ph, ymax_ph = self.phase_ax.get_ylim()
            xmin_FT, xmax_FT = self.FT_ax.get_xlim()
            ymin_FT, ymax_FT = self.FT_ax.get_ylim()

            self.phase_ax.cla()
            self.phase_ax.set_xlabel("Photon Energy [eV]")
            self.phase_ax.set_ylabel("2w Phase [rad]")
            self.FT_ax.cla()

            if cts.rabbitmode == "normal":
                self.FT_ax.set_xlabel("freq")
                self.FT_ax.set_ylabel("ampl")

                for i in range(cts.bandsnb):
                    self.FT_ax.plot(self.freqnorm, self.ampl[i] / self.ampl[i].max(),
                                    label='SB%d' % (self.SBi + 2 * i), linewidth=1)
                self.phase_ax.errorbar(self.sborder * hnu, cts.peak_phase, yerr=cts.phase_err, color='m',
                                       markersize=5)
                # self.phase_ax.errorbar(self.sborder*hnu, cts.fit_phase, yerr=cts.fit_phase_err)
                self.phase_ax.plot(np.append(self.sborder, self.sborder[-1] + 2) * hnu, self.attochirp, 'k')

            elif cts.rabbitmode == 'rainbow':
                self.FT_ax.set_xlabel("Photon Energy [eV]")
                self.FT_ax.set_ylabel("2w amplitude [arb. u.]")

                self.phase_ax.errorbar(self.energy, cts.peak_phase, yerr=0, color='m',
                                       markersize=4)
                # self.phase_ax.errorbar(self.energy, cts.fit_phase, yerr=cts.fit_phase_err)
                # self.FT_ax.plot(self.energy, cts.FT_ampl[:, self.fpeak_index] / cts.FT_ampl.max())
                self.FT_ax.plot(self.energy, cts.peak, color='m', label='2w ampl.')
                # self.FT_ax.plot(self.energy, cts.fit_ampl)

                if self.showbands_cb.isChecked() and not cts.np.array_equal(
                        cts.bands_vect, np.zeros([cts.bandsnb, 2])):
                    for i in range(2 * cts.bandsnb):
                        xval = cts.bands_vect[i // 2, i % 2]
                        if xval != 0:
                            self.phase_ax.axvline(xval, color='k')
                            self.FT_ax.axvline(xval, color='k')

                if self.int_rab_cb.isChecked():
                    delay_int_rab = cts.rabbit_mat.sum(axis=0)
                    # rescaling to compare with cts.FT_ampl
                    delay_int_rab = delay_int_rab / delay_int_rab.max() * cts.FT_ampl.max()
                    self.FT_ax.plot(self.energy, delay_int_rab, label="delay integrated")

            if keep_lims:
                self.phase_ax.set_xlim(xmin_ph, xmax_ph)
                self.phase_ax.set_ylim(ymin_ph, ymax_ph)
                self.FT_ax.set_xlim(xmin_FT, xmax_FT)
                self.FT_ax.set_ylim(ymin_FT, ymax_FT)

            if cts.calibloaded:
                phase_secax = self.phase_ax.secondary_xaxis(
                    'top', functions=(af.eeV2tof2, af.tof2EeV2))
                phase_secax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
                phase_secax.set_xlabel('TOF [s]')
                phase_secax.tick_params(labelsize=8)

                FT_secax = self.FT_ax.secondary_xaxis(
                    'top', functions=(af.eeV2tof2, af.tof2EeV2))
                FT_secax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
                FT_secax.set_xlabel('TOF [s]')
                FT_secax.tick_params(labelsize=8)

            self.FT_ax.legend(fontsize=10)
            self.phase_fc.draw()
            self.FT_fc.draw()

    def selectbands_lr(self):
        """Connected to self.selectbands_cb (CheckBox in the "Select Bands" QGroupBox)
        """
        if self.selectbands_cb.isChecked():
            self.showbands_cb.setChecked(True)
            self.showbands_cb.setEnabled(False)
        else:
            self.showbands_cb.setEnabled(True)

    def showbands_lr(self):
        """Connected to self.showbands_cb (CheckBox in the "Select Bands" QGroupBox)
        """
        self.refreshplot(keep_lims=True)

    def done_lr(self):
        """Connected to self.done_btn (PushButton in the "Select Bands" QGroupBox)
        """
        self.done_btn.setEnabled(False)
        self.normalrab_btn.setEnabled(True)
        self.bandselected = True

        cts.contrast = np.zeros([cts.bandsnb, 3])
        for i in range(cts.bandsnb):
            cts.contrast[i, 0] = int(self.SBi + 2 * i)
        self.window().updateglobvar_fn()

        if self.xuvonlyloaded:
            self.subXUV_btn.setEnabled(True)

    def clearbands_lr(self):
        """Connected to self.clearbands_btn (PushButton in the "Select Bands" QGroupBox)
        """
        cts.bands_vect = np.zeros([cts.bandsnb, 2])
        self.done_btn.setEnabled(False)
        self.clickcount = 0
        self.clearbands_btn.setEnabled(False)
        self.refreshplot(keep_lims=True)

    def lock_2w_lr(self):
        """Connected to self.lock_2w_cb (PushButton in the "RABBIT" QGroupBox)
        """
        if self.lock_2w_cb.isChecked():
            self.set_2w_le.setEnabled(True)
        else:
            self.set_2w_le.setEnabled(False)

    def set_2w_lr(self):
        """Connected to self.set_2w_cb (PushButton in the "RABBIT" QGroupBox)
        """
        cts.fpeak_lock = cts.freqnorm[int(self.set_2w_le.text())]
        self.window().updateglobvar_fn()

    def clear_lr(self):
        # completely outdated

        self.peak = []
        self.peak_phase = []
        self.fpeak = []
        self.ampl = []
        self.ang = []

        self.phase_ax.cla()
        self.FT_ax.cla()

        self.phase_fc.draw()
        self.FT_fc.draw()