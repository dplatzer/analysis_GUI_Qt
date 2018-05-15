# _calib_functions.py

from PyQt5.QtWidgets import QFileDialog
from matplotlib.ticker import FormatStrFormatter

import traceback, sys

from scipy import optimize as opt
import numpy as np

import glob_var as cts
import analysis_functions as af

class calib_functions_mixin:

    def loadfile_lr(self):
        ''' "load data" button listener'''
        if self.dataloaded:
            self.clear_lr()

        TOFfile = QFileDialog.getOpenFileName(self, 'Load data')
        fname = TOFfile[0]
        if (fname):

            self.counts = self.import_tof(fname)

            if cts.minus_sign:
                self.counts[:, 1] = (-1) * self.counts[:, 1]

            self.thxmin = 0
            self.thy = 0
            self.thxmax = self.counts[-1, 0]

            self.tof_ax.plot(self.counts[:, 0], self.counts[:, 1])

            self.tof_ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
            self.tof_fc.draw()

            self.rmbgnd_btn.setEnabled(True)
            self.findpeaks_btn.setEnabled(True)
            self.rmpeaks_btn.setEnabled(False)
            self.addpeak_cb.setEnabled(False)
            self.showexppeaks_cb.setEnabled(True)
            self.setth_cb.setEnabled(True)
            self.Y_rb.setEnabled(True)
            self.xmin_rb.setEnabled(True)
            self.xmax_rb.setEnabled(True)
            self.clear_btn.setEnabled(True)
            self.minus_sign_cb.setEnabled(True)
            self.importcalib_btn.setEnabled(True)
            self.withsb_cb.setEnabled(True)

            self.dataloaded = True
            self.threshxminBool = False
            self.threshxmaxBool = False
            self.threshyBool = False
            self.peaksfound = False

    def minus_sign_lr(self):
        if self.minus_sign_cb.isChecked():
            cts.minus_sign = True
        else:
            cts.minus_sign = False

        self.window().updateglobvar_fn()
        self.counts[:, 1] = (-1) * self.counts[:, 1]

        self.bgndremoved = False
        self.peaksfound = False

        self.reset_thresholds()

        self.rmpeaks_btn.setEnabled(False)
        self.tof2en_btn.setEnabled(False)
        if self.dataloaded:
            self.rmbgnd_btn.setEnabled(True)
        else:
            self.rmbgnd_btn.setEnabled(False)

        self.refreshplot_fn()

    def rmbgnd_lr(self):
        ''' "remove bgnd" button listener'''
        if self.counts.shape[0] < 1000: # SE10
            self.counts[:, 1] = self.counts[:, 1] - self.counts[10:20, 1].mean()
        else: # SE1
            self.counts[:, 1] = self.counts[:, 1] - self.counts[300:400, 1].mean()

        self.bgndremoved = True
        self.peaksfound = False

        self.reset_thresholds()

        self.rmpeaks_btn.setEnabled(False)
        self.tof2en_btn.setEnabled(False)
        self.rmbgnd_btn.setEnabled(False)
        self.refreshplot_fn()

    def reset_thresholds(self):
        self.threshyBool = False
        self.threshxminBool = False
        self.threshxmaxBool = False
        self.thxmin = 0
        self.thy = 0
        self.thxmax = self.counts[-1, 0]

    def update_cts_fn(self):
        ''' Updates the experimental parameters (and aguess, t0guess, cguess)'''

        try:
            cts.cur_nu = cts.C / (float(self.wvlength_le.text()) * 1e-9)
            cts.cur_Vp = float(self.retpot_le.text())
            cts.cur_L = float(self.toflength_le.text())
            cts.first_harm = int(self.firstharm_le.text())
            self.window().updateglobvar_fn()
            self.aguess = (0.5 * cts.ME * cts.cur_L ** 2 / cts.QE) / (cts.HEV * cts.cur_nu)
            # NB: cts.QE is used to convert the energy from eV to Joules. It's not the electron's charge
            self.t0guess = 5.8e-8
            self.cguess = (cts.cur_Ip + cts.cur_Vp) / (cts.HEV * cts.cur_nu)
            self.aguess_lb.setText("{:.3e}".format(self.aguess))
            self.t0guess_lb.setText("{:.3e}".format(self.t0guess))
            self.cguess_lb.setText("{:.3f}".format(self.cguess))
        except ValueError:
            print('Incorrect value. Must be a number')
        except ZeroDivisionError:
            print('ZeroDivisionError')

    ''' "find peaks" button listener '''
    def findpeaks_lr(self):

        nbpts_data = self.counts.shape[0]
        nbpts_pow = int(np.log(nbpts_data)/np.log(2)) # I take the highest power of two inferior to the number of points
        nbpts = 2**nbpts_pow
        print(nbpts)

        try:
            ip, dp, convolution = af.find_local_maxima(self.counts[:, 1], self.thy, self.thxmin, self.thxmax, nbpts,
                                                       int(self.sm1_le.text()), int(self.sm2_le.text()),
                                                       int(self.mindt_le.text()))
            #self.maximaIndices = (ip * cts.TOF_resolution).tolist()  # 1e-9 on SE1
            self.maximaIndices = self.counts[ip,0].tolist()
            self.maximaIntensity = dp.tolist()
            self.convolvedsignal = convolution.tolist()
            self.peaksfound = True
            self.rmpeaks_btn.setEnabled(True)
            self.addpeak_cb.setEnabled(True)
            self.tof2en_btn.setEnabled(True)
            self.withsb_cb.setEnabled(True)
            self.calibloaded = False
            self.refreshplot_fn()
        except ValueError:
            print('incorrect value for find local maxima paramaters')

    ''' "TOF to energy" button listener '''
    def tof2en_lr(self):
        if self.calibloaded:
            self.p_opt = [0, 0, 0]
            self.p_opt[0] = cts.afit
            self.p_opt[1] = cts.t0fit
            self.p_opt[2] = cts.cfit
            self.calibloaded = False
        else:

            qqi = float(self.firstharm_le.text())
            tofHHG = self.maximaIndices
            if self.withsb_bool:
                qq = np.arange(qqi, qqi + len(tofHHG), 1)
            else:
                qq = np.arange(qqi, qqi + 2 * len(tofHHG), 2)
            qq = qq[::-1]

            self.p_opt, self.pcov = opt.curve_fit(af.tof2EeV, tofHHG, qq,
                                                  np.array([self.aguess, self.t0guess, self.cguess]))
            toffit = np.linspace(tofHHG[0], tofHHG[-1], 100)
            self.fit_ax.plot(toffit, af.tof2EeV(toffit, self.p_opt[0], self.p_opt[1], self.p_opt[2]), label='Fit')
            self.fit_ax.plot(toffit, af.tof2EeV(toffit, self.aguess, self.t0guess, self.cguess), label='Guess')
            self.fit_ax.plot(tofHHG, qq, 'ys', label='TOF(peaks)')
        cts.afit = self.p_opt[0]
        cts.t0fit = self.p_opt[1]
        cts.cfit = self.p_opt[2]
        self.exportcalib_btn.setEnabled(True)

        self.calibBool = True

        self.update_fitpar_fn()
        self.update_envect_fn()
        self.window().updateglobvar_fn()

        tof = self.counts[:, 0]
        tof = tof[::-1]
        signal = self.counts[:, 1]
        self.Eevlin, self.signal = af.jacobian_transform(tof, signal)

        self.en_ax.set_xlabel("Energy (eV)")
        self.en_ax.set_ylabel("counts (arb. units)")
        self.en_ax.plot(self.Eevlin, self.signal)
        self.en_ax.set_xlim(0, 50)

        self.exportXUV_btn.setEnabled(True)

        self.fit_fc.draw()
        self.en_fc.draw()