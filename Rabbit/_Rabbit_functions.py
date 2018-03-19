# _Rabbit_functions.py

from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtCore import Qt
from glob import glob
import traceback, sys, os
import numpy as np

import glob_var as cts
import analysis_functions as af
from Rabbit import Rainbow_win as Rw


class Rabbit_functions_mixin:

    def importdata_lr(self):
        ''' "[Import] data" listener. Loads all the tof files, one by one, converts them to energy and regroup them into the
            cts.rabbit_mat array'''
        data_f = QFileDialog.getOpenFileName(self, 'Import data')
        data_filename = data_f[0]

        if (data_filename):
            fdir, fname = os.path.split(data_filename)
            fdir = fdir + '/'
            flist = glob(fdir + '*[0-9][0-9][0-9][0-9].txt')

            if (len(flist) != 0):
                cts.stepsnb = len(flist)
                cts.delay_vect = np.linspace(0, cts.scanstep_fs * cts.stepsnb, cts.stepsnb, endpoint=False)

                data = np.transpose(self.import_tof(data_filename))
                self.tof = data[0, :]
                self.toflength = data.shape[1]
                self.tof = self.tof[::-1]

                i = 0
                for fn in flist:
                    self.window().statusBar().showMessage("Processing " + str(i))

                    data = np.transpose(self.import_tof(fn))

                    counts = data[1, :] - data[1, 300:600].mean()
                    cts.energy_vect, counts2 = af.jacobian_transform(self.tof, counts)

                    self.elength = cts.energy_vect.shape[0]

                    if i == 0:
                        cts.rabbit_mat = np.zeros([cts.stepsnb, cts.energy_vect.shape[0]])
                        print(cts.rabbit_mat.shape)
                    cts.rabbit_mat[i, :] = counts2
                    self.data_tof.append(counts)
                    i += 1

                cts.rabbit_mat = np.array(cts.rabbit_mat)
                cts.rabbit_mat = cts.rabbit_mat[::-1, :]
                self.rab_widget.refreshplot_fn(cts.energy_vect, cts.delay_vect, cts.rabbit_mat)
                self.rab_widget.colorauto_cb.setEnabled(True)
                self.rab_widget.logscale_cb.setEnabled(True)
                self.window().statusBar().showMessage("Data loaded successfully")
                self.dataloaded = True
                self.importXUV_btn.setEnabled(True)
                self.selectbands_btn.setEnabled(True)
                self.smooth_btn.setEnabled(True)
                self.normalize_btn.setEnabled(True)
                self.exportrab_btn.setEnabled(True)

                self.window().updateglobvar_fn()
            else:
                self.window().statusBar().showMessage("Incorrect directory")

    def importrabparam_lr(self):
        ''' "[Import] RABBIT param" button listener. Loads RABBIT parameters from file'''
        rabp_fname = QFileDialog.getOpenFileName(self, 'Import RABBIT parameters')
        rabp_f = rabp_fname[0]
        if (rabp_f):
            try:
                cts.elow, cts.ehigh, cts.dE, cts.stepsnb, cts.scanstep_nm, \
                cts.afit, cts.t0fit, cts.cfit, cts.first_harm = np.loadtxt(rabp_f)

                cts.stepsnb = int(cts.stepsnb)
                cts.first_harm = int(cts.first_harm)
                self.scanparam_le.setText(str(cts.scanstep_nm))
                self.update_scanparam()
                self.elow_le.setText(str(cts.elow))
                self.ehigh_le.setText(str(cts.ehigh))
                self.dE_le.setText(str(cts.dE))
                cts.delay_vect = np.linspace(0, cts.scanstep_fs * cts.stepsnb, cts.stepsnb, endpoint=False)
                cts.energy_vect = np.arange(cts.elow, cts.ehigh, cts.dE)
                self.reset_btn()
                self.importrab_btn.setEnabled(True)

                self.window().updateglobvar_fn()
            except ValueError:
                self.window().statusBar().showMessage("Incorrect RABBIT parameters file")

    def importrab_lr(self):
        ''' "[Import] Rabbit" button listener. Loads RABBIT trace from file'''
        rab_fname = QFileDialog.getOpenFileName(self, 'Import RABBIT counts')
        rab_f = rab_fname[0]
        if (rab_f):
            try:
                cts.rabbit_mat = np.loadtxt(rab_f)

                self.rab_widget.refreshplot_fn(cts.energy_vect, cts.delay_vect, cts.rabbit_mat)
                self.dataloaded = True

                self.reset_btn()
                self.importrab_btn.setEnabled(True)
                self.exportrab_btn.setEnabled(True)
                self.selectbands_btn.setEnabled(True)
                self.normalize_btn.setEnabled(True)
                self.smooth_btn.setEnabled(True)
                self.importXUV_btn.setEnabled(True)
                self.rab_widget.colorauto_cb.setEnabled(True)
                self.rab_widget.logscale_cb.setEnabled(True)

                self.window().updateglobvar_fn()
            except ValueError:
                self.window().statusBar().showMessage("Incorrect RABBIT counts file")
            except TypeError:
                self.window().statusBar().showMessage("Incorrect RABBIT counts file")

    def importXUV_lr(self):
        ''' "[Import] XUV only" button listener. Loads XUV only from file (it is converted from tof to energy)'''
        self.update_envect_fn()
        xuv_fname = QFileDialog.getOpenFileName(self, "Import XUV only")
        xuv_f = xuv_fname[0]
        if (xuv_f):
            data = []
            try:
                data = np.loadtxt(xuv_f, unpack=True)
                self.tof = data[0, :]
                self.tof = self.tof[::-1]
                cts.xuvonly_vect = data[1, :] - data[1, 300:600].mean()
                cts.energy_vect, cts.xuvonly_vect = af.jacobian_transform(self.tof, cts.xuvonly_vect)
                cts.xuvonly_vect = np.array(cts.xuvonly_vect)
                self.xuvonlyloaded = True
                if self.bandselected:
                    self.subXUV_btn.setEnabled(True)

                self.window().updateglobvar_fn()

            except ValueError:
                self.window().statusBar().showMessage("Incorrect XUV only file")

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
        for i in range(cts.stepsnb):
            cts.rabbit_mat[i, :] = cts.rabbit_mat[i, :] / cts.rabbit_mat[i, :].sum()
        self.rab_widget.colorauto_cb.setCheckState(Qt.Checked)
        self.rab_widget.logscale_cb.setCheckState(Qt.Unchecked)

        self.normalize_btn.setText("normalized")
        self.normalize_btn.setEnabled(False)
        cts.rabnormalized = True
        self.window().updateglobvar_fn()

        self.rab_widget.refreshplot_fn(signal=cts.rabbit_mat)

    def normalrab_lr(self):
        ''' ANALYSIS - "Normal RABBIT" button listener'''
        try:
            cts.rabbitmode = "normal"
            self.FT_ax.cla()
            self.FT_ax.set_xlabel("freq")
            self.FT_ax.set_ylabel("ampl")

            hnu = cts.HEV * cts.cur_nu
            self.SBi = cts.first_harm + 1
            self.peak = np.zeros(cts.bandsnb)

            i = 0
            if cts.xuvsubstracted:
                signal = cts.rabbitxuvsub_mat
            else:
                signal = cts.rabbit_mat

            jx = []
            for xi, xf in cts.bands_vect:
                ji = np.argmin(abs(cts.energy_vect - xi))
                jf = np.argmin(abs(cts.energy_vect - xf)) + 1
                jx.append([ji, jf])

            for ji, jf in jx[:]:
                if jf - ji > 1: #if we integrate ie if we have more than one point
                    x_data = np.trapz(signal[:, ji:jf], cts.energy_vect[ji:jf], axis=1)
                else:
                    x_data = signal[:, ji]
                fpeak, peak, peak_phase = af.find_2w(x_data, cts.scanstep_fs)
                self.freqnorm, ampl, ang = af.FFT(x_data, cts.scanstep_fs)
                ampl = np.array(ampl)

                self.FT_ax.plot(self.freqnorm, ampl/ampl.max(), label='SB%d' %(self.SBi + 2*i), linewidth=1)
                self.peak[i] = np.absolute(peak)
                self.peak_phase.append(peak_phase)
                self.fpeak.append(fpeak)
                self.ampl.append(ampl)
                self.ang.append(ang)
                i += 1

            cts.FT_ampl = np.array(self.ampl)
            cts.FT_phase = np.array(self.ang)
            cts.fpeak = np.array(self.fpeak)
            cts.peak = np.array(self.peak)
            cts.peak_phase = - np.array(self.peak_phase)
            cts.peak_phase = np.unwrap(cts.peak_phase)

            fpeak_main_i = np.argmax(cts.peak[:]/cts.FT_ampl[:, 0]) #between 0 and cts.bandsnb - 1
            cts.fpeak_main = cts.fpeak[fpeak_main_i]
            self.fpeak_index = np.argmin(abs(self.freqnorm - cts.fpeak_main)) #between 0 and len(self.freqnorm)

            for i in range(cts.bandsnb):
                cts.peak[i] = cts.FT_ampl[i, self.fpeak_index]
                cts.peak_phase[i] = -1*(cts.FT_phase[i, self.fpeak_index] + cts.two_w_phioffset)
                cts.contrast[i, 2] = 2 * cts.peak[i] / cts.FT_ampl[i, 0]
            cts.peak_phase = np.unwrap(cts.peak_phase)

            leg = self.FT_ax.legend()
            for label in leg.get_texts():
                label.set_fontsize(8)
            for label in leg.get_lines():
                label.set_linewidth(1)

            self.FT_ax.relim()
            self.FT_ax.autoscale()
            self.FT_fc.draw()

            sborder = np.linspace(self.SBi, self.SBi + 2*(cts.bandsnb-1), cts.bandsnb)
            self.pa = np.polyfit(sborder, cts.peak_phase, 1)
            attochirp = self.pa[0]*sborder + self.pa[1]

            self.phase_ax.cla()
            self.phase_ax.set_xlabel("SB")
            self.phase_ax.set_ylabel("Phase Diff")
            self.phase_ax.plot(sborder, cts.peak_phase, 'm', marker = 's')
            self.phase_ax.plot(sborder, attochirp, 'k')
            self.phase_fc.draw()

            cts.chirp_as = self.pa[0]/(cts.cur_nu * 2*np.pi)/2 * 1e18
            cts.chirp_as_eV = float(cts.chirp_as/(hnu))

            self.rainbowrab_btn.setEnabled(True)
            self.FTcontrast_btn.setEnabled(True)
            self.plotSBvsdelay_btn.setEnabled(True)

            self.window().updateglobvar_fn()
        except:
            print(traceback.format_exception(*sys.exc_info()))

    def rainbowrab_lr(self):
        ''' ANALYSIS - "Rainbow RABBIT" button listener'''
        try:
            cts.rabbitmode = "rainbow"
            hnu = cts.HEV * cts.cur_nu
            self.SBi = cts.first_harm + 1
            i = 0

            jx = []
            for xi, xf in cts.bands_vect:
                ji = np.argmin(abs(cts.energy_vect - xi))
                jf = np.argmin(abs(cts.energy_vect - xf))
                jx.append([ji, jf])
            jx = np.array(jx)

            if cts.xuvsubstracted:
                signal = cts.rabbitxuvsub_mat
            else:
                signal = cts.rabbit_mat

            j=0
            for ii, jj in jx[:]:
                for x in range(ii,jj):
                    x_data = signal[:, x]
                    f2om, peak, peak_phase = af.find_2w(x_data, cts.scanstep_fs)
                    self.freqnorm, ampl, ang =  af.FFT(x_data, cts.scanstep_fs)

                    self.ang_rainbow.append(ang)
                    self.ampl_rainbow.append(ampl/ampl.max())
                    self.ampl_rainbow2.append(ampl/ampl.max())
                    self.peak_rainbow.append(np.absolute(peak))
                    self.energy_rainbow.append(cts.energy_vect[x])
                    self.energy_rainbow2.append(cts.energy_vect[x])
                    self.peak_phase_rainbow.append(peak_phase)
                    self.fpeak_rainbow.append(f2om)

                #filling between bands for plot
                if j<(cts.bandsnb - 1):
                    for k in range(jx[j, 1], jx[j+1, 0]):
                        self.energy_rainbow2.append(cts.energy_vect[k])
                        self.ampl_rainbow2.append([1e-10]*ampl.shape[0])
                j+=1
            cts.FT_ampl = np.array(self.ampl_rainbow)
            cts.FT_ampl = np.transpose(cts.FT_ampl)
            cts.FT_phase = np.array(self.ang_rainbow)
            cts.FT_phase = np.transpose(cts.FT_phase)

            self.ampl_rainbow2 = np.array(self.ampl_rainbow2)
            self.ampl_rainbow2 = np.transpose(self.ampl_rainbow2)
            cts.fpeak = np.array(self.fpeak_rainbow)
            cts.peak = np.array(self.peak_rainbow)
            cts.peak_phase = - np.array(self.peak_phase_rainbow)

            self.energy_rainbow = np.array(self.energy_rainbow)
            self.energy_rainbow2 = np.array(self.energy_rainbow2)

            fpeak_main_i = np.argmax(cts.peak)
            cts.fpeak_main = cts.fpeak[fpeak_main_i]
            self.fpeak_index = np.argmin(abs(self.freqnorm - cts.fpeak_main))

            for i in range(len(cts.peak)):
                cts.peak[i] = cts.FT_ampl[self.fpeak_index, i]
                cts.peak_phase[i] =  -1*cts.FT_phase[self.fpeak_index, i]
            cts.peak_phase = np.unwrap(cts.peak_phase)

            self.phase_ax.cla()
            self.phase_ax.xaxis.set_ticks(np.arange((cts.elow+0.5)//hnu, (cts.ehigh+0.5)//hnu), 2)
            self.phase_ax.set_xlabel("Harmonic order")
            self.phase_ax.set_ylabel("Phase (rad)")
            self.phase_ax.plot(self.energy_rainbow/hnu, cts.peak_phase, 'ms', markersize=2)
            self.phase_ax.plot(self.energy_rainbow/hnu, self.pa[0]*self.energy_rainbow/hnu + self.pa[1])
            self.phase_fc.draw()

            self.FT_ax.cla()
            self.FT_ax.set_xlabel("Harmonic order")
            self.FT_ax.set_ylabel("frequency")
            self.FT_ax.imshow(self.ampl_rainbow2, extent=(self.energy_rainbow2[0]/hnu,self.energy_rainbow2[-1]/hnu,
                                                          self.freqnorm[0], self.freqnorm[-1]), aspect='auto',
                                                        origin='lower', interpolation='nearest')
            self.FT_fc.draw()
            self.window().updateglobvar_fn()

            rw = Rw.RainbowWin(self)
        except Exception:
            print(traceback.format_exception(*sys.exc_info()))