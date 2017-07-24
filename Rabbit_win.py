import sys, os, traceback
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QPushButton, QApplication, QGridLayout, QHBoxLayout, QVBoxLayout, QGroupBox, QLabel
from PyQt5.QtWidgets import QLineEdit, QCheckBox, QRadioButton, QButtonGroup, QFileDialog, QDialog, QTableWidget
from PyQt5.QtWidgets import QTableWidgetItem, QTableView
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
from glob import glob
import numpy as np
from scipy import optimize as opt

# Homemade modules
import glob_var as cts
import analysis_functions as af
import other_widgets as ow

class RabbitWin(QWidget):
    """ Rabbit window

        For the names of the children widgets, I tried to put suffixes that indicate clearly their types:
        *_btn -> QPushButton,
        *_le -> QLineEdit,
        *_lb -> QLabel,
        *layout -> QHBoxLayout, QVBoxLayout or QGridLayout,
        *_box -> QGroupBox,
        *_cb -> QCheckBox,
        *_rb -> QRadioButton

        The functions that are connected to a widget's event have the suffix _lr (for 'listener'). For example,
        a button named test_btn will be connected to a function test_lr.
        Some functions may be connected to widgets but without the suffix _lr in their names. It means that they
        are not only called when interacting with the widget.
        """
    def __init__(self, parent=None) -> None:
        """Initialization of the window

                the main layout is called mainLayout. It is divided in two:
                    - graphLayout: the left part, contains all the figures
                    - commandLayout: the right part, contains all the buttons, fields, checkboxes...
                Both graphLayout and commandLayout are divided into sub-layouts.

                This function calls several functions to initialize each part of the window.
                The name of these functions has the shape 'init_*layout'."""
        super(RabbitWin, self).__init__(parent=parent)
        self.setWindowTitle("RABBIT")
        self.mainlayout = QHBoxLayout()
        self.graphlayout = QVBoxLayout()
        self.commandLayout = QVBoxLayout()
        self.commandLayout.setSpacing(10)

        self.init_var()
        self.init_importlayout()
        self.init_envectlayout()
        self.init_sigtreatmentlayout()
        self.init_rabbitlayout()
        self.init_exportlayout()
        self.init_plotbtnlayout()
        self.init_graphlayout()

        self.mainlayout.addLayout(self.graphlayout)
        self.mainlayout.addLayout(self.commandLayout)
        self.setLayout(self.mainlayout)
        self.show()

    ''' Initialization of instance attributes'''
    def init_var(self) -> None:
        self.dataloaded = False
        self.xuvonlyloaded = False
        self.bandselected = False

    ''' In commandLayout - Initialization of the "Import" section'''
    def init_importlayout(self) -> None:
        Importlayout = QGridLayout()
        Importlayout.setSpacing(10)

        Import_box = QGroupBox(self)
        Import_box.setTitle("Import")
        Import_box.setFixedSize(300, 100)

        self.importcalib_btn = QPushButton("calib", self)
        self.importdata_btn = QPushButton("data", self)
        self.importXUV_btn = QPushButton("XUV only", self)
        self.importrabparam_btn = QPushButton("RABBIT param", self)
        self.importrab_btn = QPushButton("RABBIT", self)

        self.importcalib_btn.clicked.connect(self.importcalib_lr)
        self.importdata_btn.clicked.connect(self.importdata_lr)
        self.importXUV_btn.clicked.connect(self.importXUV_lr)
        self.importrabparam_btn.clicked.connect(self.importrabparam_lr)
        self.importrab_btn.clicked.connect(self.importrab_lr)

        Importlayout.addWidget(self.importcalib_btn, 0, 0)
        Importlayout.addWidget(self.importdata_btn, 1, 0)
        Importlayout.addWidget(self.importXUV_btn, 0, 1)
        Importlayout.addWidget(self.importrabparam_btn, 0, 2)
        Importlayout.addWidget(self.importrab_btn, 1, 2)

        Import_box.setLayout(Importlayout)
        self.commandLayout.addWidget(Import_box)

        for widget in Import_box.children():
            if isinstance(widget, QPushButton):
                widget.setSizePolicy(0, 0)
                widget.setEnabled(False)
        self.importcalib_btn.setEnabled(True)
        self.importrabparam_btn.setEnabled(True)

    ''' In commandLayout - Initialization of the energy vector section, with elow, ehigh and dE'''
    def init_envectlayout(self) -> None:
        paramlayout = QHBoxLayout()
        envectlayout = QGridLayout()
        envectlayout.setSpacing(10)
        envect_box = QGroupBox()
        envect_box.setTitle("Energy vector parameters")
        envect_box.setSizePolicy(0, 0)

        self.elow_le = QLineEdit("{:.2f}".format(cts.elow), self)
        self.ehigh_le = QLineEdit("{:.2f}".format(cts.ehigh), self)
        self.dE_le = QLineEdit(str(cts.dE), self)

        self.elow_le.returnPressed.connect(self.update_envect_fn)
        self.ehigh_le.returnPressed.connect(self.update_envect_fn)
        self.dE_le.returnPressed.connect(self.update_envect_fn)

        envectlayout.addWidget(QLabel("E low (eV)"), 0, 0)
        envectlayout.addWidget(self.elow_le, 1, 0)
        envectlayout.addWidget(QLabel("E high (eV)"), 0, 1)
        envectlayout.addWidget(self.ehigh_le, 1, 1)
        envectlayout.addWidget(QLabel("dE (eV)"), 0, 2)
        envectlayout.addWidget(self.dE_le, 1, 2)

        envect_box.setLayout(envectlayout)

        for widget in envect_box.children():
            if isinstance(widget, QLabel) or isinstance(widget, QLineEdit):
                widget.setSizePolicy(0, 0)
                widget.setFixedSize(55, 20)

        paramlayout.addWidget(envect_box)

        scanparamlayout = QVBoxLayout()

        self.scanparam_le = QLineEdit(str(cts.scanstep_nm))
        self.scanparam_le.setSizePolicy(0, 0)
        self.scanparam_le.setFixedSize(55, 20)
        self.scanparam_le.returnPressed.connect(self.update_scanparam)
        label = QLabel("scan steps (nm)")
        label.setSizePolicy(0, 0)
        label.setFixedSize(80, 20)

        scanparamlayout.addWidget(label)
        scanparamlayout.addWidget(self.scanparam_le)

        paramlayout.addLayout(scanparamlayout)

        self.commandLayout.addLayout(paramlayout)

    ''' In commandLayout - Initialization of the "Signal Treatment" section'''
    def init_sigtreatmentlayout(self) -> None:
        cts.bandsnb = 5
        sigtreatment_box = QGroupBox("Signal Treatment", self)
        sigtreatmentlayout = QGridLayout()

        sigtreatment_box.setFixedSize(300, 100)

        self.smooth_btn = QPushButton("smooth", self)
        self.smooth_le = QLineEdit("2", self)
        self.normalize_btn = QPushButton("normalize", self)
        self.subXUV_btn = QPushButton("substract XUV", self)
        self.selectbands_btn = QPushButton("select bands")
        self.bandsnb_le = QLineEdit(str(cts.bandsnb), self)

        self.smooth_btn.clicked.connect(self.smoothrab_lr)
        self.normalize_btn.clicked.connect(self.normalizerab_lr)
        self.selectbands_btn.clicked.connect(self.selectbands_lr)
        self.bandsnb_le.returnPressed.connect(self.bandsnb_lr)
        self.subXUV_btn.clicked.connect(self.subXUV_lr)

        sigtreatmentlayout.addWidget(self.smooth_btn, 0, 0)
        sigtreatmentlayout.addWidget(self.smooth_le, 1, 0)
        sigtreatmentlayout.addWidget(self.normalize_btn, 0, 1)
        sigtreatmentlayout.addWidget(self.subXUV_btn, 1, 1)
        sigtreatmentlayout.addWidget(self.selectbands_btn, 0, 2)
        sigtreatmentlayout.addWidget(self.bandsnb_le, 1, 2)

        sigtreatment_box.setLayout(sigtreatmentlayout)

        for widget in sigtreatment_box.children():
            if isinstance(widget, QPushButton):
                widget.setSizePolicy(0, 0)
                widget.setEnabled(False)
            if isinstance(widget, QLineEdit):
                widget.setFixedSize(55, 20)

        self.commandLayout.addWidget(sigtreatment_box)

    ''' In commandLayout - Initialization of the "RABBIT" section'''
    def init_rabbitlayout(self) -> None:
        rabbitlayout = QGridLayout()
        rabbit_box = QGroupBox("RABBIT", self)
        rabbit_box.setFixedSize(300, 60)

        self.normalrab_btn = QPushButton("Normal", self)
        self.FTcontrast_btn = QPushButton("FT/Contrast", self)
        self.rainbowrab_btn = QPushButton("Rainbow", self)

        self.normalrab_btn.clicked.connect(self.normalrab_lr)
        self.FTcontrast_btn.clicked.connect(self.FTcontrast_lr)
        self.rainbowrab_btn.clicked.connect(self.rainbowrab_lr)

        rabbitlayout.addWidget(self.normalrab_btn, 0, 0)
        rabbitlayout.addWidget(self.FTcontrast_btn, 0, 1)
        rabbitlayout.addWidget(self.rainbowrab_btn, 0, 2)

        rabbit_box.setLayout(rabbitlayout)

        for widget in rabbit_box.children():
            if isinstance(widget, QPushButton):
                widget.setSizePolicy(0, 0)
                widget.setEnabled(False)

        self.commandLayout.addWidget(rabbit_box)

    ''' In commandLayout - Initialization of the "Export" section'''
    def init_exportlayout(self) -> None:
        exportlayout = QGridLayout()
        export_box = QGroupBox("Export", self)
        export_box.setFixedSize(300, 60)

        self.exportrab_btn = QPushButton("RABBIT", self)

        self.exportrab_btn.clicked.connect(self.exportrab_lr)

        exportlayout.addWidget(self.exportrab_btn)

        export_box.setLayout(exportlayout)

        for widget in export_box.children():
            if isinstance(widget, QPushButton):
                widget.setSizePolicy(0, 0)
                widget.setEnabled(False)
        self.commandLayout.addWidget(export_box)

    ''' In commandLayout - Initialization of the "Plot" section'''
    def init_plotbtnlayout(self) -> None:
        plotbtnlayout = QGridLayout()
        plotbtn_box = QGroupBox("Plot", self)
        plotbtn_box.setFixedSize(300, 60)

        self.plotSBvsdelay_btn = QPushButton("SB vs delay", self)

        self.plotSBvsdelay_btn.clicked.connect(self.plotSBvsdelay_lr)

        plotbtnlayout.addWidget(self.plotSBvsdelay_btn)
        plotbtn_box.setLayout(plotbtnlayout)

        for widget in plotbtn_box.children():
            if isinstance(widget, QPushButton):
                widget.setSizePolicy(0, 0)
                widget.setEnabled(False)
        self.commandLayout.addWidget(plotbtn_box)

    ''' In graphLayout - Initialization of the 3 figures'''
    def init_graphlayout(self) -> None:
        self.rab_widget = ow.plot3DWidget(self) # new class defined in other_widgets.py
        self.rab_widget.xlabel = "E (eV)"
        self.rab_widget.ylabel = "t (fs)"

        self.graphlayout.addWidget(self.rab_widget)

        self.phaseFTlayout = QHBoxLayout()
        self.phaselayout = QVBoxLayout()

        phase_fig = Figure(figsize=(2, 2), dpi=100)
        self.phase_fc = FigureCanvas(phase_fig)
        self.phase_fc.setSizePolicy(1, 0)
        self.phase_ax = self.phase_fc.figure.add_subplot(111)
        self.phase_ax.tick_params(labelsize = 8)
        nav = NavigationToolbar2QT(self.phase_fc, self)
        nav.setStyleSheet("QToolBar { border: 0px }")
        self.phaselayout.addWidget(self.phase_fc)
        self.phaselayout.addWidget(nav)
        self.phase_fc.draw()

        self.FTlayout = QVBoxLayout()
        FT_fig = Figure(figsize=(2, 2), dpi=100)
        self.FT_fc = FigureCanvas(FT_fig)
        self.FT_fc.setSizePolicy(1, 0)
        self.FT_ax = self.FT_fc.figure.add_subplot(111)
        self.FT_ax.tick_params(labelsize = 8)
        nav2 = NavigationToolbar2QT(self.FT_fc, self)
        nav2.setStyleSheet("QToolBar { border: 0px }")
        self.FTlayout.addWidget(self.FT_fc)
        self.FTlayout.addWidget(nav2)
        self.FT_fc.draw()

        self.phaseFTlayout.addLayout(self.phaselayout)
        self.phaseFTlayout.addLayout(self.FTlayout)

        self.graphlayout.addLayout(self.phaseFTlayout)

    ''' From a calibration file, loading afit, t0fit, cfit and the first harmonic'''
    def importcalib_lr(self) -> None:

        calib_fname = QFileDialog.getOpenFileName(self, 'Import calibration')
        calib_f = calib_fname[0]
        if (calib_f):
            with open(calib_f, 'r') as file:
                f = file.read().splitlines()
                try:
                    if (len(f) < 4):
                        raise IndexError
                    cts.afit, cts.t0fit, cts.cfit, cts.first_harm = [float(f[i]) for i in range(4)]
                    cts.first_harm = int(cts.first_harm)
                    self.window().updateglobvar_fn()
                    self.reset_btn()
                    self.importdata_btn.setEnabled(True)
                except IndexError:
                    print('Not enough data in the calib file. Needed: a_fit, t0_fit, c_fit and first harmonic')
                except ValueError:
                    print("Incorrect calibration data")

    ''' "[Import] data" listener. Loads all the tof files, one by one, converts them to energy and regroup them into the
    cts.rabbit_mat array'''
    def importdata_lr(self) -> None:
        data_f = QFileDialog.getOpenFileName(self, 'Import data')
        data_filename = data_f[0]

        if (data_filename):
            fdir, fname = os.path.split(data_filename)
            fdir = fdir + '/'
            flist = glob(fdir + '*delay_0*')

            if (len(flist) != 0):
                cts.stepsnb = len(flist)
                cts.delay_vect = np.linspace(0, cts.scanstep_fs * cts.stepsnb, cts.stepsnb, endpoint=False)

                try:
                    data = np.loadtxt(data_filename, unpack=True)
                except ValueError:  # if we don't choose a correct file but in the right directory
                    data = np.loadtxt(flist[0], unpack=True)

                self.tof = data[0, :]
                self.toflength = data.shape[1]
                self.tof = self.tof[::-1]
                self.data_tof = []
                cts.rabbit_mat = []

                i = 0
                for fn in flist:
                    self.window().statusBar().showMessage("Processing " + str(i))
                    data = np.loadtxt(fn, unpack=True)
                    counts = data[1, :] - data[1, 300:600].mean()
                    cts.energy_vect, counts2 = af.jacobian_transform(self.tof, counts)
                    self.elength = cts.energy_vect.shape[0]
                    cts.rabbit_mat.append(counts2)
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

    ''' "[Import] RABBIT param" button listener. Loads RABBIT parameters from file'''
    def importrabparam_lr(self) -> None:
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

    ''' "[Import] Rabbit" button listener. Loads RABBIT trace from file'''
    def importrab_lr(self) -> None:
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

    ''' "[Import] XUV only" button listener. Loads XUV only from file (it is converted from tof to energy)'''
    def importXUV_lr(self) -> None:
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

    ''' "smooth RABBIT" button listener'''
    def smoothrab_lr(self) -> None:
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

    ''' "normalize RABBIT" button listener'''
    def normalizerab_lr(self) -> None:
        for i in range(cts.stepsnb):
            cts.rabbit_mat[i, :] = cts.rabbit_mat[i, :] / cts.rabbit_mat[i, :].sum()
        self.rab_widget.colorauto_cb.setCheckState(Qt.Checked)
        self.rab_widget.logscale_cb.setCheckState(Qt.Unchecked)

        self.normalize_btn.setText("normalized")
        self.normalize_btn.setEnabled(False)
        cts.rabnormalized = True
        self.window().updateglobvar_fn()

        self.rab_widget.refreshplot_fn(signal=cts.rabbit_mat)

    ''' "select bands" button listener'''
    def selectbands_lr(self) -> None:
        try:
            cts.bandsnb = int(self.bandsnb_le.text())
            cts.bands_vect = np.zeros([cts.bandsnb, 2])
            self.subXUV_btn.setEnabled(False)
            self.normalrab_btn.setEnabled(False)
            self.rainbowrab_btn.setEnabled(False)
            self.FTcontrast_btn.setEnabled(False)
            self.plotSBvsdelay_btn.setEnabled(False)
            self.window().updateglobvar_fn()
            nw = selectBandsWin(self) # new class defined below
        except ValueError:
            self.window().statusBar().showMessage("Number of bands must be an integer")

    ''' called when pressing enter in the bandsnb_le object'''
    def bandsnb_lr(self) -> None:
        try:
            cts.bandsnb = int(self.bandsnb_le.text())
            self.window().updateglobvar_fn()
        except ValueError:
            self.window().statusBar().showMessage("Number of bands must be an integer")

    ''' "substract XUV" button listener. Opens a new window'''
    def subXUV_lr(self) -> None:
        cts.xuvsubstracted = False
        sw = subXUVWin(self) # new class defined below

    ''' ANALYSIS - "Normal RABBIT" button listener'''
    def normalrab_lr(self) -> None:
        try:
            cts.rabbitmode = "normal"
            self.FT_ax.cla()
            self.FT_ax.set_xlabel("frequency (units of f0)")
            self.FT_ax.set_ylabel("normalized amplitude (arb. u.)")

            hnu = cts.HEV * cts.cur_nu
            self.SBi = cts.first_harm + 1
            self.ampl = []
            self.ang = []
            self.fpeak = []
            self.peak = np.zeros(cts.bandsnb)
            self.peak_phase = []
            i = 0
            if cts.xuvsubstracted:
                signal = cts.rabbitxuvsub_mat
            else:
                signal = cts.rabbit_mat

            jx = []
            for xi, xf in cts.bands_vect:
                ji = np.argmin(abs(cts.energy_vect - xi))
                jf = np.argmin(abs(cts.energy_vect - xf))
                jx.append([ji, jf])

            for ji, jf in jx[:]:
                x_data = np.trapz(signal[:, ji:jf], cts.energy_vect[ji:jf], axis=1)
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
            self.phase_ax.set_xlabel("Harmonic order")
            self.phase_ax.set_ylabel("Phase (rad)")
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

    ''' ANALYSIS - "Rainbow RABBIT" button listener'''
    def rainbowrab_lr(self) -> None:
        try:
            cts.rabbitmode = "rainbow"
            hnu = cts.HEV * cts.cur_nu
            self.SBi = cts.first_harm + 1
            self.ampl_rainbow = []
            self.ampl_rainbow2 = []
            self.ang_rainbow = []
            self.fpeak_rainbow = []
            self.peak_rainbow = []
            self.peak_phase_rainbow = []
            self.energy_rainbow = []
            self.energy_rainbow2 = []
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

            rw = RainbowWin(self)
        except Exception:
            print(traceback.format_exception(*sys.exc_info()))

    ''' "FT/contrast" button listener. Opens a new window'''
    def FTcontrast_lr(self) -> None:
        try:
            w = FTContrastWin(self) # new class defined below
        except Exception:
            print(traceback.format_exception(*sys.exc_info()))

    ''' "[plot] SB vs delay" button listener. Opens a new window'''
    def plotSBvsdelay_lr(self) -> None:
        try:
            w = SBvsDelayWin(self) # new class defined below
        except Exception:
            print(traceback.format_exception(*sys.exc_info()))

    ''' Updates the values of the scan steps, in nm and fs'''
    def update_scanparam(self) -> None:
        try:
            cts.scanstep_nm = float(self.scanparam_le.text())
            cts.scanstep_fs = float(self.scanparam_le.text()) * 2 / (cts.C * 1e-6)
            self.window().updateglobvar_fn()
        except ValueError:
            print('Scan step must be a number')

    ''' Updates the energy vector parameters'''
    def update_envect_fn(self) -> None:

        cts.elow = float(self.elow_le.text())
        cts.ehigh = float(self.ehigh_le.text())
        cts.dE = float(self.dE_le.text())
        self.elow_le.setText("{:.2f}".format(cts.elow))
        self.ehigh_le.setText("{:.2f}".format(cts.ehigh))
        self.window().updateglobvar_fn()

    ''' "[Export] RABBIT" button listener. Saves the rabbit in two files: *_counts for the data and *_param for 
    the parameters'''
    def exportrab_lr(self) -> None:
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

    ''' "Reset" button listener. Resets the widgets, not the variables'''
    def reset_btn(self) -> None:
        self.importXUV_btn.setEnabled(False)
        self.importdata_btn.setEnabled(False)
        self.importrab_btn.setEnabled(False)
        self.smooth_btn.setEnabled(False)
        self.normalize_btn.setEnabled(False)
        self.subXUV_btn.setEnabled(False)
        self.selectbands_btn.setEnabled(False)
        self.normalrab_btn.setEnabled(False)
        self.FTcontrast_btn.setEnabled(False)
        self.rainbowrab_btn.setEnabled(False)
        self.exportrab_btn.setEnabled(False)
        self.plotSBvsdelay_btn.setEnabled(False)

        self.rab_widget.colorauto_cb.setEnabled(False)
        self.rab_widget.logscale_cb.setEnabled(False)

''' created when clicking on the "Select bands" button'''
class selectBandsWin(QDialog):
    def __init__(self, parent=RabbitWin):
        super(selectBandsWin, self).__init__(parent)
        self.layout = QHBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setGeometry(100, 100, 1000, 800)
        self.setWindowFlags(Qt.Window)

        win = selectBandsPlotWin(self) # new class defined below
        self.layout.addWidget(win)

        self.setLayout(self.layout)
        self.show()

''' used in the "selectBandsWin" object above'''
class selectBandsPlotWin(ow.plot3DWidget):
    def __init__(self, parent) -> None:
        super(selectBandsPlotWin, self).__init__(parent)

        self.init_vars()
        self.init_layout()

        self.logscale_cb.stateChanged.connect(self.refreshplot_sb)  # NB: logscale_cb still connected to logscale_lr in
        # the plot3DWidget __init__ function (and has to be).
        self.color_sld.valueChanged[int].disconnect() # for this widget, it doesn't work if I don't disconnect first from
        # the plot3DWidget function
        self.color_sld.valueChanged[int].connect(self.sbcolor_lr)
        self.fc.mpl_connect('button_press_event', self.onclick)
        self.refreshplot_fn()

        self.show()

    def init_vars(self) -> None:
        self.clickcount = 0
        self.x = cts.energy_vect
        self.y = cts.delay_vect
        self.signal = cts.rabbit_mat
        self.value = 25

    def init_layout(self) -> None:
        self.paramlayout = QHBoxLayout()

        self.color_bg = QButtonGroup(self)
        self.white_rb = QRadioButton("white", self)
        self.white_rb.value = "w"
        self.white_rb.toggled.connect(self.color_lr)

        self.black_rb = QRadioButton("black", self)
        self.black_rb.value = "k"
        self.black_rb.toggled.connect(self.color_lr)

        self.color_bg.addButton(self.black_rb)
        self.color_bg.addButton(self.white_rb)

        self.selectbands_cb = QCheckBox("Select bands", self)
        self.selectbands_cb.toggle()

        self.selectmode_bg = QButtonGroup(self)
        self.manualselect_rb = QRadioButton("manual", self)
        self.manualselect_rb.value = "manual"
        self.manualselect_rb.toggled.connect(self.selectmode_lr)

        self.FTselect_rb = QRadioButton("with FT", self)
        self.FTselect_rb.value = "FT"
        self.FTselect_rb.toggled.connect(self.selectmode_lr)

        self.selectmode_bg.addButton(self.manualselect_rb)
        self.selectmode_bg.addButton(self.FTselect_rb)

        self.FTthreshold_le = QLineEdit("0.2", self)
        self.FTthreshold_le.setSizePolicy(0, 0)
        self.FTthreshold_le.setFixedSize(60, 20)
        self.FTthreshold_le.returnPressed.connect(self.FTselect_fn)

        self.done_btn = QPushButton("Done", self)
        self.done_btn.setSizePolicy(0, 0)
        self.done_btn.clicked.connect(self.done_lr)

        self.cancel_btn = QPushButton("Cancel", self)
        self.cancel_btn.setSizePolicy(0, 0)
        self.cancel_btn.clicked.connect(self.cancel_lr)

        self.white_rb.toggle()
        self.manualselect_rb.toggle()

        self.paramlayout.addWidget(self.black_rb)
        self.paramlayout.addWidget(self.white_rb)
        self.paramlayout.addWidget(self.selectbands_cb)
        self.paramlayout.addWidget(self.manualselect_rb)
        self.paramlayout.addWidget(self.FTselect_rb)
        self.paramlayout.addWidget(self.FTthreshold_le)
        self.paramlayout.addStretch(1)
        self.paramlayout.addWidget(self.cancel_btn)
        self.paramlayout.addWidget(self.done_btn)


        self.optionslayout.addLayout(self.paramlayout)
        self.setLayout(self.mainlayout)
        self.mainlayout.setContentsMargins(10, 10, 10, 10)
    # Here I override a method in QDialog that closes the Dialog on pressing return or enter
    # It is willingly empty to avoid closing the window when changing the threshold of the FT selection
    def keyPressEvent(self, event) -> None:
        key = event.key()

    def sbcolor_lr(self, value) -> None:
        self.refreshplot_sb(value)
        self.value = value

    def selectmode_lr(self) -> None:
        rb = self.sender()
        self.selectmode = rb.value
        if self.selectmode == "manual":
            cts.bands_vect = np.zeros([cts.bandsnb, 2])
            self.FTthreshold_le.setEnabled(False)
            self.done_btn.setEnabled(False)
            self.refreshplot_sb(value=self.value)
        elif self.selectmode == "FT":
            self.FTthreshold_le.setEnabled(True)
            self.FTselect_fn()
            self.done_btn.setEnabled(True)

    def FTselect_fn(self) -> None:
        amp_threshold = float(self.FTthreshold_le.text())
        hnu = cts.HEV * cts.cur_nu
        self.SBi = cts.first_harm + 1

        self.ampl_rainbow = []
        self.ang_rainbow = []
        self.fpeak_rainbow = []
        self.peak_rainbow = []
        self.peak_phase_rainbow = []
        self.energy_rainbow = []

        for x in range(cts.energy_vect.shape[0]):
            if cts.xuvsubstracted:
                xdata = cts.rabbitxuvsub_mat[:, x]
            else:
                xdata = cts.rabbit_mat[:, x]

            f2om, peak, peak_phase = af.find_2w(xdata, cts.scanstep_fs)
            self.freqnorm, ampl, ang = af.FFT(xdata, cts.scanstep_fs)
            self.ampl_rainbow.append(ampl / ampl.max())
            self.peak_rainbow.append(np.absolute(peak))
            self.peak_phase_rainbow.append(peak_phase)

            self.fpeak_rainbow.append(f2om)
        self.ampl_rainbow = np.array(self.ampl_rainbow)
        self.ampl_rainbow = np.transpose(self.ampl_rainbow)
        self.fpeak_rainbow = np.array(self.fpeak_rainbow)
        self.peak_rainbow = np.array(self.peak_rainbow)
        self.peak_phase_rainbow = - np.array(self.peak_phase_rainbow)

        for i in range(cts.bandsnb):
            iemin = np.argmin(abs(cts.energy_vect - (self.SBi - 0.6 + 2 * i) * hnu))
            iemax = np.argmin(abs(cts.energy_vect - (self.SBi + 0.6 + 2 * i) * hnu))

            ampl = self.peak_rainbow[iemin:iemax]
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
            cts.bands_vect[i, 1] = cts.energy_vect[xi[-1] + iemin + 1]

        self.refreshplot_sb(value=self.value)

    def onclick(self, event) -> None:
        if self.selectbands_cb.isChecked() and self.selectmode == "manual":
            if self.clickcount < 2 * cts.bandsnb:
                if self.clickcount == 2 * cts.bandsnb - 1:
                    self.done_btn.setEnabled(True)
                cts.bands_vect[self.clickcount // 2, self.clickcount % 2] = event.xdata
                self.clickcount += 1
            # when all the lines are set, we can modify a line on clicking next to it
            else:
                k = np.argmin(abs(cts.bands_vect - event.xdata))
                cts.bands_vect[k // 2, k % 2] = event.xdata
            self.refreshplot_sb(value=self.value)

    def refreshplot_sb(self, value=25) -> None:
        try:
            self.refreshplot_fn(value=value)

            for i in range(2 * cts.bandsnb):
                xval = cts.bands_vect[i // 2, i % 2]
                if xval != 0:
                    line = np.full([cts.stepsnb, 2], xval)
                    line[:, 1] = cts.delay_vect[:]
                    self.ax.plot(line[:, 0], line[:, 1], color=self.color)

            self.fc.draw()
        except Exception:
            print(traceback.format_exception(*sys.exc_info()))

    def color_lr(self) -> None:
        rb = self.sender()
        self.color = rb.value
        self.refreshplot_sb(value=self.value)

    def done_lr(self) -> None:
        self.parent().parent().window().updateglobvar_fn()
        self.parent().parent().normalrab_btn.setEnabled(True)
        self.parent().parent().bandselected = True

        cts.contrast = np.zeros([cts.bandsnb, 3])
        for i in range(cts.bandsnb):
            cts.contrast[i, 0] = int(cts.first_harm + 2 * i + 1)
        self.parent().parent().window().updateglobvar_fn()

        if(self.parent().parent().xuvonlyloaded):
            self.parent().parent().subXUV_btn.setEnabled(True)
        self.parent().destroy()

    def cancel_lr(self) -> None:
        cts.bands_vect = np.zeros([cts.bandsnb, 2])
        self.parent().parent().window().updateglobvar_fn()
        self.parent().parent().plotSBvsdelay_btn.setEnabled(False)
        self.parent().destroy()

#I define these classes here so the class RabbitWin doesn't have to be called from another module.

''' created when clicking on the "Substract XUV" button"'''
class subXUVWin(QDialog):
    def __init__(self, parent=RabbitWin) -> None:
        super(subXUVWin, self).__init__(parent)
        self.layout = QHBoxLayout()
        #self.layout.setContentsMargins(0, 0, 0, 0)
        self.setWindowFlags(Qt.Window)

        cts.xuvsubstracted = False
        self.init_layout()
        self.rabgraph.refreshplot_fn(cts.energy_vect, cts.delay_vect, cts.rabbit_mat)

        self.setLayout(self.layout)
        self.show()

    def init_layout(self) -> None:
        self.rabgraph = ow.plot3DWidget(self)

        self.graph2layout = QVBoxLayout()

        fig2 = Figure(figsize=(4, 3), dpi=100)
        self.fc2 = FigureCanvas(fig2)
        self.ax2 = self.fc2.figure.add_subplot(111)
        self.fc2.draw()
        nav2 = NavigationToolbar2QT(self.fc2, self)
        nav2.setStyleSheet("QToolBar { border: 0px }")

        self.graph2layout.addWidget(self.fc2)
        self.graph2layout.addWidget(nav2)

        self.layout.addWidget(self.rabgraph)
        self.layout.addLayout(self.graph2layout)

        self.btn_layout = QHBoxLayout()

        self.substract_btn = QPushButton("Substract", self)
        self.cancel_btn = QPushButton("Cancel", self)
        self.done_btn = QPushButton("Done", self)

        self.substract_btn.setSizePolicy(0, 0)
        self.cancel_btn.setSizePolicy(0, 0)
        self.done_btn.setSizePolicy(0, 0)

        self.substract_btn.clicked.connect(self.substract_lr)
        self.cancel_btn.clicked.connect(self.cancel_lr)
        self.done_btn.clicked.connect(self.done_lr)

        self.btn_layout.addWidget(self.substract_btn)
        self.btn_layout.addWidget(self.cancel_btn)
        self.btn_layout.addWidget(self.done_btn)

        self.graph2layout.addLayout(self.btn_layout)

    def substract_lr(self) -> None:
        try:
            hnu = cts.HEV * cts.cur_nu
            cts.rabbitxuvsub_mat = np.array(cts.rabbit_mat) # the np.array() synthax must not be removed:
            # Here if we write a = b, changing a will change b whereas in a = np.array(b) b won't change

            if cts.rabnormalized:
                cts.xuvonly_vect = cts.xuvonly_vect/cts.xuvonly_vect.sum()
            if cts.rabsmoothed:
                sm = int(self.parent().smooth_le.text())
                cts.xuvonly_vect = af.smooth(cts.xuvonly_vect, sm)

            for i in range(cts.bandsnb):
                Harmonic_energy = (cts.first_harm + 2 * i)* hnu
                ii = np.argmin(abs(Harmonic_energy - hnu/2 - cts.energy_vect))
                jj = np.argmin(abs(Harmonic_energy + hnu/2 - cts.energy_vect))

                for k in range(cts.stepsnb):
                    maxrab = cts.rabbit_mat[k, np.where(abs(cts.energy_vect - Harmonic_energy) < hnu/2)].max()
                    maxxuv = cts.xuvonly_vect[np.where(abs(cts.energy_vect - Harmonic_energy) < hnu / 2)].max()
                    cts.rabbitxuvsub_mat[k, ii:jj] = cts.rabbit_mat[k, ii:jj] - cts.xuvonly_vect[ii:jj] * maxrab/maxxuv

            self.rabgraph.refreshplot_fn(signal=cts.rabbitxuvsub_mat)

            cts.xuvsubstracted = True

            self.ax2.set_xlabel("E (eV)")
            self.ax2.set_ylabel("t (fs)")
            integral_original = np.trapz(cts.rabbit_mat[:,:], axis=0)
            integral = np.trapz(cts.rabbitxuvsub_mat[:,:], axis=0)
            self.ax2.set_xlim([cts.energy_vect[0], cts.energy_vect[-1]])

            self.ax2.plot(cts.energy_vect, integral, label="substracted")
            self.ax2.plot(cts.energy_vect, integral_original, label="original")
            self.ax2.plot(cts.energy_vect, cts.xuvonly_vect*cts.stepsnb, label="XUV only")
            self.ax2.legend()
            self.fc2.draw()

        except Exception:
            print(traceback.format_exception(*sys.exc_info()))

    def cancel_lr(self) -> None:
        cts.rabbitxuvsub_mat = np.zeros([1,1])
        self.parent().window().updateglobvar_fn()
        self.destroy()

    def done_lr(self) -> None:
        self.parent().window().updateglobvar_fn()
        if cts.xuvsubstracted:
            self.parent().rab_widget.refreshplot_fn(cts.energy_vect, cts.delay_vect, cts.rabbitxuvsub_mat)
        else:
            self.parent().rab_widget.refreshplot_fn(cts.energy_vect, cts.delay_vect, cts.rabbit_mat)

        self.destroy()

''' created when clicking on the "FT/Contrast" button'''
class FTContrastWin(QDialog):
    def __init__(self, parent=RabbitWin) -> None:
        super(FTContrastWin, self).__init__(parent)

        self.setWindowFlags(Qt.Window)

        self.par = parent
        self.mainlayout = QHBoxLayout()

        self.init_var()
        self.init_graphLayout()
        self.init_two_wlayout()
        self.init_FTlayout()
        self.init_peaksTable()

        self.mainlayout.addLayout(self.commandLayout)
        self.setLayout(self.mainlayout)
        self.show()

    def init_var(self) -> None:
        self.SBi = int(cts.first_harm) + 1
        self.contrast = np.zeros(cts.bandsnb)
        self.contrast_int = np.zeros(cts.bandsnb)

        self.jx = []
        for xi, xf in cts.bands_vect:
            ji = np.argmin(abs(cts.energy_vect - xi))
            jf = np.argmin(abs(cts.energy_vect - xf))
            self.jx.append([ji, jf])

        self.fpeak = np.array(cts.fpeak)
        self.ampl = np.array(cts.FT_ampl)
        self.peak_phase = np.array(cts.peak_phase)
        self.peak = np.array(cts.peak)
        self.freqnorm = np.array(self.par.freqnorm)
        self.ang = np.array(cts.FT_phase)
        self.fpeak_main = cts.fpeak_main
        self.fpeak_index = self.par.fpeak_index

    def init_graphLayout(self) -> None:
        self.graphlayout = QVBoxLayout()

        fig = Figure(figsize=(4, 3), dpi=100)
        self.fc = FigureCanvas(fig)
        self.ax1 = self.fc.figure.add_subplot(211)
        self.ax2 = self.fc.figure.add_subplot(212)
        self.fc.draw()
        nav = NavigationToolbar2QT(self.fc, self)
        nav.setStyleSheet("QToolBar { border: 0px }")

        self.ax1.set_ylabel("Amplitude")
        self.ax1.tick_params(labelsize=8)

        self.ax2.set_xlabel("Frequency (units of \u03C9)")
        self.ax2.set_ylabel("Phase (rad)")
        self.ax2.tick_params(labelsize=8)

        self.graphlayout.addWidget(self.fc)
        self.graphlayout.addWidget(nav)

        self.cblayout = QHBoxLayout()

        self.cb = []
        for i in range(cts.bandsnb):
            cb = QCheckBox("SB " + str(self.SBi + 2 * i), self)
            cb.stateChanged.connect(self.cb_lr)
            self.cblayout.addWidget(cb)
            self.cb.append(cb)

        self.graphlayout.addLayout(self.cblayout)
        self.mainlayout.addLayout(self.graphlayout)

        for cb in self.cb:
            cb.setCheckState(Qt.Checked)

        self.refreshplot_fn()

    def init_two_wlayout(self) -> None:
        self.commandLayout = QVBoxLayout()

        self.two_wlayout = QGridLayout()
        self.two_w_box = QGroupBox()
        self.two_w_box.setTitle("Find 2w")
        self.two_w_box.setSizePolicy(0, 0)
        self.two_w_box.setFixedWidth(400)

        self.bfilter_cb = QCheckBox("ButterWorth filter", self)
        self.bfilter_cb.stateChanged.connect(self.bfilter_lr)

        self.average_cb = QCheckBox("average", self)
        self.average_cb.stateChanged.connect(self.average_lr)

        self.integral_cb = QCheckBox("integral", self)
        self.integral_cb.stateChanged.connect(self.integral_lr)

        self.phioffset_le = QLineEdit(str(cts.two_w_phioffset), self)
        self.phioffset_le.setSizePolicy(0, 0)
        self.phioffset_le.setFixedSize(55, 20)
        self.phioffset_le.returnPressed.connect(self.phioffset_lr)

        self.two_w_btn = QPushButton("Find 2w", self)
        self.setSizePolicy(0, 0)
        self.two_w_btn.clicked.connect(self.find_2w_lr)

        self.two_wlayout.addWidget(self.bfilter_cb, 1, 0)
        self.two_wlayout.addWidget(self.average_cb, 1, 1)
        self.two_wlayout.addWidget(self.integral_cb, 1, 2)
        self.two_wlayout.addWidget(QLabel("phi offset"), 0, 3)
        self.two_wlayout.addWidget(self.phioffset_le, 1, 3)
        self.two_wlayout.addWidget(self.two_w_btn, 1, 4)

        self.two_w_box.setLayout(self.two_wlayout)
        self.commandLayout.addWidget(self.two_w_box)

    def init_FTlayout(self) -> None:
        self.FTlayout = QGridLayout()
        self.FT_box = QGroupBox()
        self.FT_box.setTitle("FT")
        self.FT_box.setSizePolicy(0, 0)
        self.FT_box.setFixedWidth(400)

        self.FTpadding_cb = QCheckBox("0-padding", self)
        self.FTpadding_cb.stateChanged.connect(self.FTpaddingcb_lr)

        self.FTpadding_le = QLineEdit(str(cts.FT_npad), self)
        self.FTpadding_le.returnPressed.connect(self.FTpaddingle_lr)
        self.FTpadding_le.setSizePolicy(0, 0)
        self.FTpadding_le.setFixedSize(55, 20)

        self.FTwindow_cb = QCheckBox("window", self)
        self.FTwindow_cb.stateChanged.connect(self.FTwindow_lr)

        self.FTzeroorder_cb = QCheckBox("order 0", self)
        self.FTzeroorder_cb.stateChanged.connect(self.FTzeroorder_lr)
        self.FTzeroorder_cb.setCheckState(Qt.Checked)

        self.FTdt_le = QLineEdit("%.3f" %cts.scanstep_fs, self)
        self.FTdt_le.returnPressed.connect(self.FTdt_lr)
        self.FTdt_le.setSizePolicy(0, 0)
        self.FTdt_le.setFixedSize(55, 20)

        self.FT_btn = QPushButton("FT", self)
        self.FT_btn.setSizePolicy(0, 0)
        self.FT_btn.clicked.connect(self.FT_lr)

        self.FTlayout.addWidget(self.FTpadding_cb, 0, 0)
        self.FTlayout.addWidget(self.FTpadding_le, 1, 0)
        self.FTlayout.addWidget(self.FTwindow_cb, 0, 1)
        self.FTlayout.addWidget(self.FTzeroorder_cb, 0, 2)
        self.FTlayout.addWidget(QLabel("dt (fs)"), 0, 3)
        self.FTlayout.addWidget(self.FTdt_le, 1, 3)
        self.FTlayout.addWidget(self.FT_btn, 1, 4)

        self.FT_box.setLayout(self.FTlayout)
        self.commandLayout.addWidget(self.FT_box)

    def init_peaksTable(self) -> None:
        self.peaksTable = QTableWidget()
        self.peaksTable.setSelectionBehavior(QTableView.SelectRows)

        self.commandLayout.addWidget(self.peaksTable)

        self.refreshTable_fn()

    def FT_fn(self) -> None:
        try:
           self.ampl = []
           self.ang = []
           cts.FT_npad = int(self.FTpadding_le.text())
           self.FTdt_lr() # to update cts.scanstep_fs
           for ii, jj in self.jx[:]:
               x_data = np.trapz(cts.rabbit_mat[:, ii:jj], cts.energy_vect[ii:jj], axis=1)
               self.freqnorm, ampl, ang = af.FFT(x_data, cts.scanstep_fs)
               self.ampl.append(ampl)
               self.ang.append(ang)

           cts.FT_ampl = np.array(self.ampl)
           cts.FT_phase = np.array(self.ang)

        except ValueError:
            print("npad must be an integer")
        except Exception:
            print(traceback.format_exception(*sys.exc_info()))

    def FT_lr(self) -> None:
        self.FT_fn()
        self.refreshTable_fn()
        self.refreshplot_fn()

    def find_2w_fn(self) -> None:
        self.FT_fn()
        self.fpeak = []
        self.peak = np.zeros(cts.bandsnb)
        self.peak_phase = []

        cts.two_w_phioffset = float(self.phioffset_le.text())
        i=0
        for ii, jj in self.jx[:]:
            x_data=np.trapz(cts.rabbit_mat[:, ii:jj], cts.energy_vect[ii:jj], axis=1)
            fpeak, peak, peak_phase = af.find_2w(x_data, cts.scanstep_fs)
            self.peak[i] = np.absolute(peak)
            self.peak_phase.append(peak_phase)
            self.fpeak.append(fpeak)
            i+=1

        cts.fpeak = np.array(self.fpeak)
        cts.peak = np.array(self.peak)
        cts.peak_phase = np.array(self.peak_phase)

        fpeak_main_i = np.argmax(cts.peak[:]/cts.FT_ampl[:, 0])
        cts.fpeak_main = cts.fpeak[fpeak_main_i]
        self.fpeak_index = np.argmin(abs(self.freqnorm - cts.fpeak_main))

        for i in range(cts.bandsnb):
            cts.peak[i] = cts.FT_ampl[i, self.fpeak_index]
            cts.peak_phase[i] = -1*(cts.FT_phase[i, self.fpeak_index] + cts.two_w_phioffset)
        cts.peak_phase = np.unwrap(cts.peak_phase)

    def find_2w_lr(self) -> None:
        try:
            self.find_2w_fn()
            self.refreshTable_fn()
            self.refreshplot_fn()
        except Exception:
            print(traceback.format_exception(*sys.exc_info()))

    def refreshTable_fn(self) -> None:
        self.peaksTable.clear()
        self.peaksTable.setRowCount(cts.bandsnb + 1)
        self.peaksTable.setColumnCount(5)

        self.peaksTable.setItem(0, 0, QTableWidgetItem("SB"))
        self.peaksTable.setItem(0, 1, QTableWidgetItem("2*(2w/0)"))
        self.peaksTable.setItem(0, 2, QTableWidgetItem("integral"))
        self.peaksTable.setItem(0, 3, QTableWidgetItem("fpeak"))
        self.peaksTable.setItem(0, 4, QTableWidgetItem("phase at %.3f" % cts.fpeak_main))

        df = self.freqnorm[1] * cts.cur_nu

        for i in range(cts.bandsnb):
            self.contrast[i] = 2 * cts.peak[i] / cts.FT_ampl[i, 0]
            cts.contrast[i,2] = 2 * cts.peak[i] / cts.FT_ampl[i, 0]
            self.contrast_int[i] = 2 * np.trapz(cts.FT_ampl[i, self.fpeak_index - 1:self.fpeak_index + 1], dx=df) / \
                                   np.trapz(cts.FT_ampl[i, 0:2], dx=df)

            self.peaksTable.setItem(i+1, 0, QTableWidgetItem(str(self.SBi + 2 * i)))
            self.peaksTable.setItem(i+1, 1, QTableWidgetItem("{:.3f}".format(self.contrast[i])))
            self.peaksTable.setItem(i+1, 2, QTableWidgetItem("{:.3f}".format(self.contrast_int[i])))
            self.peaksTable.setItem(i+1, 3, QTableWidgetItem("{:.3f}".format(self.fpeak[i])))
            self.peaksTable.setItem(i+1, 4, QTableWidgetItem("{:.3f}".format(self.peak_phase[i])))

        w = 0
        for i in range(self.peaksTable.columnCount()):
            w = w + self.peaksTable.columnWidth(i)
        # print(w)
        h = 0
        for i in range(self.peaksTable.rowCount()):
            h = h + self.peaksTable.rowHeight(i)
        self.peaksTable.setFixedHeight(h + 25)
        self.peaksTable.setFixedWidth(w + 20)

    def refreshplot_fn(self) -> None:
        self.ax1.cla()
        self.ax2.cla()
        for i in range(cts.bandsnb):
            if self.cb[i].isChecked():
                self.ax1.plot(self.freqnorm, cts.FT_ampl[i, :] / cts.FT_ampl[i, 0], label="SB%d" % (self.SBi + 2 * i),
                              linewidth=1)
                self.ax2.plot(self.freqnorm, cts.FT_phase[i, :], label="SB%d" % (self.SBi + 2 * i), linewidth=1)

        leg = self.ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                            ncol=7, mode="expand", borderaxespad=0.)
        for label in leg.get_texts():
            label.set_fontsize(8)
        for label in leg.get_lines():
            label.set_linewidth(1)

        self.fc.draw()

    def phioffset_lr(self) -> None:
        cts.two_w_phioffset = float(self.phioffset_le.text())
        self.parent().window().updateglobvar_fn()

    def integral_lr(self) -> None:
        if self.integral_cb.isChecked():
            cts.two_w_integral = True
            if self.average_cb.isChecked():
                self.average_cb.setCheckState(Qt.Unchecked)
        else:
            cts.two_w_integral = False
        self.parent().window().updateglobvar_fn()

    def average_lr(self) -> None:
        if self.average_cb.isChecked():
            cts.two_w_average = True
            if self.integral_cb.isChecked():
                self.integral_cb.setCheckState(Qt.Unchecked)
        else:
            cts.two_w_average = False
        self.parent().window().updateglobvar_fn()

    def bfilter_lr(self) -> None:
        if self.bfilter_cb.isChecked():
            cts.two_w_bfilter = True
        else:
            cts.two_w_bfilter = False
        self.parent().window().updateglobvar_fn()

    def FTpaddingle_lr(self) -> None:
        try:
            cts.FT_npad = int(self.FTpadding_le.text())
            self.parent().window().updateglobvar_fn()
        except ValueError:
            print("npad must be an integer")

    def FTpaddingcb_lr(self) -> None:
        if self.FTpadding_cb.isChecked():
            cts.FT_padding = True
        else:
            cts.FT_padding = False
        self.parent().window().updateglobvar_fn()

    def FTwindow_lr(self) -> None:
        if self.FTwindow_cb.isChecked():
            cts.FT_window = True
        else:
            cts.FT_window = False
        self.parent().window().updateglobvar_fn()

    def FTzeroorder_lr(self) -> None:
        if self.FTzeroorder_cb.isChecked():
            cts.FT_zero_order = True
        else:
            cts.FT_zero_order = False
        self.parent().window().updateglobvar_fn()

    def FTdt_lr(self) -> None:
        cts.scanstep_fs = float(self.FTdt_le.text())
        cts.scanstep_nm = cts.scanstep_fs / 2 * (cts.C * 1e-6)
        self.parent().window().updateglobvar_fn()

    def cb_lr(self) -> None:
        self.refreshplot_fn()

    def keyPressEvent(self, event) -> None:
        key = event.key()

''' created at the end of the Rainbow Rabbit Treatment'''
class RainbowWin(QDialog):
    def __init__(self, parent=RabbitWin) -> None:
        super(RainbowWin, self).__init__(parent)
        self.setGeometry(300, 300, 600, 400)
        self.setWindowFlags(Qt.Window)

        self.par = parent
        self.mainlayout = QHBoxLayout()

        self.init_var()
        self.init_graphlayout()
        self.init_commandlayout()
        self.refreshplot_fn()

        self.setLayout(self.mainlayout)
        self.show()

    def init_var(self) -> None:
        self.plotachirp = True
        self.peak_phase = np.array(cts.peak_phase)

    def init_graphlayout(self) -> None:
        self.graphlayout = QVBoxLayout()

        fig = Figure(figsize=(4, 3), dpi=100)
        self.fc = FigureCanvas(fig)
        self.ax1 = self.fc.figure.add_subplot(211)
        self.ax2 = self.fc.figure.add_subplot(212)
        self.fc.draw()
        nav = NavigationToolbar2QT(self.fc, self)
        nav.setStyleSheet("QToolBar { border: 0px }")

        self.graphlayout.addWidget(self.fc)
        self.graphlayout.addWidget(nav)

        self.mainlayout.addLayout(self.graphlayout)

    def init_commandlayout(self) -> None:
        self.commandLayout = QHBoxLayout()

        self.subattochirp_btn = QPushButton("Sub. achirp", self)
        self.reset_btn = QPushButton("reset", self)
        self.export_btn = QPushButton("Export", self)

        self.subattochirp_btn.clicked.connect(self.subattochirp_lr)
        self.reset_btn.clicked.connect(self.reset_lr)
        self.export_btn.clicked.connect(self.export_lr)

        self.commandLayout.addWidget(self.subattochirp_btn)
        self.commandLayout.addWidget(self.reset_btn)
        self.commandLayout.addWidget(self.export_btn)
        self.mainlayout.addLayout(self.commandLayout)

    def subattochirp_lr(self) -> None:
        hnu = cts.HEV * cts.cur_nu
        self.peak_phase = cts.peak_phase - (self.par.pa[0] * self.par.energy_rainbow/hnu +\
                                           self.par.pa[1])
        self.plotachirp = False
        self.refreshplot_fn()
        self.subattochirp_btn.setEnabled(False)

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

    def reset_lr(self) -> None:
        self.peak_phase = np.array(cts.peak_phase)
        self.plotachirp = True
        self.subattochirp_btn.setEnabled(True)
        self.refreshplot_fn()

    def refreshplot_fn(self) -> None:
        hnu = cts.HEV * cts.cur_nu
        self.ax1.cla()
        self.ax2.cla()

        self.ax1.set_xlabel("Harmonic order")
        self.ax1.set_ylabel("Phase (rad)")
        self.ax1.tick_params(labelsize=8)

        self.ax2.set_ylabel("Amplitude (arb. u.)")
        self.ax2.tick_params(labelsize=8)

        self.ax1.plot(self.par.energy_rainbow/hnu, cts.peak, 'bs', markersize=2)
        self.ax2.plot(self.par.energy_rainbow/hnu, self.peak_phase, 'ms', markersize=2)
        if(self.plotachirp):
            self.ax2.plot(self.par.energy_rainbow/hnu, self.par.pa[0]*\
                          self.par.energy_rainbow/hnu + self.par.pa[1])
        self.fc.draw()

''' created when clicking on the "[plot] SB vs delay" button'''
class SBvsDelayWin(QDialog):
    def __init__(self, parent=RabbitWin) -> None:
        super(SBvsDelayWin, self).__init__(parent)
        self.setGeometry(300, 300, 1000, 500)
        self.setWindowFlags(Qt.Window)

        self.par = parent
        self.mainlayout = QHBoxLayout()

        self.init_var()
        self.init_graphlayout()
        self.init_commandlayout()
        self.refreshTable_fn()
        self.refreshplot_fn()

        self.setLayout(self.mainlayout)
        self.show()

    def init_var(self) -> None:

        self.jx = []
        for xi, xf in cts.bands_vect:
            ji = np.argmin(abs(cts.energy_vect - xi))
            jf = np.argmin(abs(cts.energy_vect - xf))
            self.jx.append([ji, jf])

        self.SB = []
        for ji, jf in self.jx[:]:
            x_data = np.trapz(cts.rabbit_mat[:, ji:jf], cts.energy_vect[ji:jf], axis=1)
            self.SB.append(x_data)
        self.SB = np.array(self.SB)

        for i in range(cts.bandsnb):
            popt, pcov = opt.curve_fit(self.cosfit_fn, cts.delay_vect, self.SB[i,:], np.array([1e-11, 1e-12, 0]))

            if popt[1] < 0:
                popt[1] = -1 * popt[1]
                popt[2] = popt[2] + np.pi

            cts.contrast[i, 1] = popt[1] / popt[0]

        self.parent().window().updateglobvar_fn()

    def init_graphlayout(self) -> None:
        self.graphlayout = QVBoxLayout()

        fig = Figure(figsize=(4, 3), dpi=100)
        self.fc = FigureCanvas(fig)
        self.ax = self.fc.figure.add_subplot(111)
        self.fc.draw()
        nav = NavigationToolbar2QT(self.fc, self)
        nav.setStyleSheet("QToolBar { border: 0px }")

        self.checkbox_layout = QHBoxLayout()
        self.checkbox_table = []
        for i in range(cts.bandsnb):
            cb = QCheckBox(str(cts.first_harm + 2*i +1), self)
            cb.setCheckState(Qt.Checked)
            cb.stateChanged.connect(self.refreshplot_fn)
            self.checkbox_layout.addWidget(cb)
            self.checkbox_table.append(cb)

        self.graphlayout.addWidget(self.fc)
        self.graphlayout.addWidget(nav)
        self.graphlayout.addLayout(self.checkbox_layout)

        self.mainlayout.addLayout(self.graphlayout)

    def init_commandlayout(self) -> None:
        self.commandLayout = QVBoxLayout()

        fig = Figure(figsize=(4, 3), dpi=100)
        self.contrast_fc = FigureCanvas(fig)
        self.contrast_ax = self.contrast_fc.figure.add_subplot(111)
        self.contrast_fc.draw()
        nav = NavigationToolbar2QT(self.contrast_fc, self)
        nav.setStyleSheet("QToolBar { border: 0px }")

        self.commandLayout.addWidget(self.contrast_fc)
        self.commandLayout.addWidget(nav)

        self.contrastTable = QTableWidget()
        self.contrastTable.setSelectionBehavior(QTableView.SelectRows)
        self.contrastTable.horizontalHeader().setDefaultSectionSize(90)
        self.contrastTable.verticalHeader().setDefaultSectionSize(30)

        self.commandLayout.addWidget(self.contrastTable)
        self.mainlayout.addLayout(self.commandLayout)

    def refreshTable_fn(self) -> None:
        self.contrastTable.clear()

        self.contrastTable.setRowCount(cts.bandsnb+1)
        self.contrastTable.setColumnCount(3)

        self.contrastTable.setItem(0, 0, QTableWidgetItem("SB"))
        self.contrastTable.setItem(0, 1, QTableWidgetItem("cos fit contrast"))
        self.contrastTable.setItem(0, 2, QTableWidgetItem("FT contrast"))

        for i in range(cts.bandsnb):
            self.contrastTable.setItem(i + 1, 0, QTableWidgetItem(str(cts.contrast[i, 0])))
            self.contrastTable.setItem(i + 1, 1, QTableWidgetItem("{:.3f}".format(cts.contrast[i, 1])))
            self.contrastTable.setItem(i + 1, 2, QTableWidgetItem("{:.3f}".format(cts.contrast[i, 2])))

        w = 0
        for i in range(self.contrastTable.columnCount()):
            w = w + self.contrastTable.columnWidth(i)
        h = 0
        for i in range(self.contrastTable.rowCount()):
            h = h + self.contrastTable.rowHeight(i)
        self.contrastTable.setFixedHeight(h + 25)
        self.contrastTable.setFixedWidth(w + 20)

    def refreshplot_fn(self) -> None:
       try:
            self.ax.cla()
            self.ax.set_xlabel("delay", fontsize=10)
            self.ax.set_ylabel("SB Intensity", fontsize=10)

            for i in range(cts.bandsnb):
                if self.checkbox_table[i].isChecked():
                    self.ax.plot(cts.delay_vect, self.SB[i,:], label="SB%d" % (cts.first_harm + 2 * i + 1))

            leg = self.ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                                  ncol=7, mode="expand", borderaxespad=0.)
            for label in leg.get_texts():
                label.set_fontsize(8)
            for label in leg.get_lines():
                label.set_linewidth(1)

            self.fc.draw()

            self.contrast_ax.cla()
            self.contrast_ax.set_xlabel("SB order", fontsize=10)
            self.contrast_ax.set_ylabel("Contrast", fontsize=10)

            self.contrast_ax.plot(cts.contrast[:, 0], cts.contrast[:, 1], 'rs', label="cos fit")
            self.contrast_ax.plot(cts.contrast[:, 0], cts.contrast[:, 2], 'b+', label="FT fit")

            contrast_leg = self.contrast_ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                                  ncol=7, mode="expand", borderaxespad=0.)
            for label in contrast_leg.get_texts():
                label.set_fontsize(8)
            for label in contrast_leg.get_lines():
                label.set_linewidth(1)
       except Exception:
           print(traceback.format_exception(*sys.exc_info()))

    def cosfit_fn(self, delay, a0, a1, phi):
        return a0 + a1 * np.cos(2* cts.fpeak_main * np.pi * cts.cur_nu * delay * 1e-15 + phi)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = RabbitWin()
    sys.exit(app.exec_())
