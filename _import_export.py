# _import_export.py

import os, traceback, sys
from glob import glob
from PyQt5.QtWidgets import QFileDialog
from matplotlib.ticker import FormatStrFormatter

import glob_var as cts
import numpy as np
import analysis_functions as af


class Imp_Exp_Mixin:

    def import_tof(self, filename):
        with open(filename, 'r') as file:
            data = []
            try:
                data = [[float(digit) for digit in line.split()] for line in file]
                '''
                if cts.decimal_separ == 'dot':
                    data = [[float(digit) for digit in line.split()] for line in file]
                elif cts.decimal_separ == 'comma':
                    data = [[float(digit.replace(',', '.')) for digit in line.split()] for line in file]
                else:
                    print("Error in the code")'''

            except ValueError: #if the decimal separator is a comma
                data = [[float(digit.replace(',', '.')) for digit in line.split()] for line in file]

            except IndexError:
                print('Incorrect data file')
                self.window().statusBar().showMessage('Incorrect data file')

            except Exception:
                print(traceback.format_exception(*sys.exc_info()))

            finally:
                data_array = np.asarray(data)  # converting 1D list into 2D array

        return data_array


    # called by CalibWin or RabbitWin
    def importcalib_lr(self):
        #From a calibration file, loading afit, t0fit, cfit and the first harmonic
        # this function can be call either from the calibration tab or the rabbit tab. This means "self" represents
        # either an instance of CalibWin or RabbitWin.

        calib_tab = self.window().calib_tab
        rabbit_tab = self.window().rabbit_tab

        calib_fname = QFileDialog.getOpenFileName(self, 'Import calibration')
        calib_f = calib_fname[0]
        if (calib_f):
            with open(calib_f, 'r') as file:
                f = file.read().splitlines()
                try:
                    if (len(f) < 3):
                        raise IndexError
                    cts.afit, cts.t0fit, cts.cfit = [float(f[i]) for i in range(3)]

                    # calib_win part
                    calib_tab.calibloaded = True
                    calib_tab.calibBool = True
                    calib_tab.tof2en_btn.setEnabled(True)
                    calib_tab.update_fitpar_fn()

                    # rabbit_win part
                    rabbit_tab.reset_btn()
                    rabbit_tab.importdata_btn.setEnabled(True)

                    self.window().updateglobvar_fn()

                except IndexError:
                    print('Not enough data in the calib file. Needed: a_fit, t0_fit and c_fit')
                except ValueError:
                    print("Incorrect calibration data")

    # called by CalibWin
    def exportcalib_lr(self):
        # Exporting afit, t0fit, cfit and the first harmonic in a file after choosing its name and location

        filename = QFileDialog.getSaveFileName(self, 'Save XUV')
        fname = filename[0]
        if fname:
            fit_param = self.p_opt
            np.savetxt(fname, fit_param, fmt='%1.4e', delimiter='\t')

    # called by RabbitWin
    def exportrab_lr(self):
        ''' "[Export] RABBIT" button listener. Saves the rabbit in two files: *_counts for the data and *_param for
        the parameters'''
        rab_fname = QFileDialog.getSaveFileName(self)
        rab_f = rab_fname[0]
        try:
            if (rab_f):
                self.rabbit_param = np.array([cts.elow, cts.ehigh, cts.dE, cts.stepsnb, cts.scanstep_nm,
                                              cts.afit, cts.t0fit, cts.cfit, cts.first_harm])
                np.savetxt((rab_f + "_param.txt"), self.rabbit_param, fmt='%.5e', delimiter='\t')
                np.savetxt((rab_f + "_counts.txt"), cts.rabbit_mat, fmt='%.5e', delimiter='\t')
        except Exception:
            print(traceback.format_exception(*sys.exc_info()))

    # called by RainbowWin (_Rabbit_win_subwidgets.RainbowWin)
    def export_lr(self) -> None:
        filenamesave = QFileDialog.getSaveFileName(self)
        fname = filenamesave[0]
        try:
            if fname:
                en=[]
                amp=[]
                phase=[]
                header = "FT_padding = " + str(cts.FT_padding) + ", " + \
                         "FT_window = " + str(cts.FT_window) + ", " + \
                         "FT_zero_order = " + str(cts.FT_zero_order) + ", " + \
                         "FT_npad = " + str(cts.FT_npad) + "\n" + \
                         "two_w_average = " + str(cts.two_w_average) + ", " + \
                         "two_w_bfilter = " + str(cts.two_w_bfilter) + ", " + \
                         "two_w_integral = " + str(cts.two_w_integral) + ", " + \
                         "two_w_phioffset = " + str(cts.two_w_phioffset) + "\n\n" + \
                         "energy (eV)\t2w_amplitude\t2w_phase"
                k=0
                for i in range(self.par.energy_rainbow.shape[0]):
                    if (self.par.energy_rainbow[i] - self.par.energy_rainbow[i - 1] >\
                                    2 * cts.dE or i == self.par.energy_rainbow.shape[0] - 1):
                        l = np.array([en, amp, phase])
                        l = np.transpose(l)
                        np.savetxt((fname + "_SB" + str(cts.first_harm + 1 + 2 * k) + ".txt"), l, fmt='%.5e',
                                   delimiter='\t', header=header)
                        del en[:]
                        del amp[:]
                        del phase[:]
                        k += 1
                    en.append(self.par.energy_rainbow[i])
                    amp.append(cts.peak[i])
                    phase.append(cts.peak_phase[i])

        except Exception:
            print(traceback.format_exception(*sys.exc_info()))