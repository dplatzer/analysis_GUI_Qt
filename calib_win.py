import sys
from PyQt5.QtCore import QRect, Qt
from PyQt5.QtWidgets import QWidget, QPushButton, QApplication, QGridLayout, QHBoxLayout, QVBoxLayout, QGroupBox, QLabel
from PyQt5.QtWidgets import QLineEdit, QCheckBox, QComboBox, QRadioButton, QButtonGroup, QFileDialog, QDialog, QTableWidget
from PyQt5.QtWidgets import QTableWidgetItem, QTableView
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
import numpy as np
from scipy import optimize as opt
from matplotlib.ticker import FormatStrFormatter

# Homemade modules
import glob_var as cts
import analysis_functions as af

class CalibWin(QWidget):
    """ Energy calibration window

    For the names of the children widgets, I tried to put suffixes that indicate clearly their types:
    *_btn -> QPushButton,
    *_le -> QLineEdit,
    *_lb -> QLabel,
    *layout -> QHBoxLayout, QVBoxLayout or QGridLayout,
    *_box -> QGroupBox,
    *_cb -> QCheckBox,
    *_rb -> QRadioButton

    The functions that are connected to a widget's event have the suffix _lr (for 'listener'). For example, a button named
    test_btn will be connected to a function test_lr.
    Some functions may be connected to widgets but without the suffix _lr in their names. It means that they are not only
    called when interacting with the widget.
    """
    def __init__(self, parent=None) -> None:
        """Initialization of the window

        the main layout is called mainLayout. It is divided in two:
            - graphLayout: the left part, contains all the figures
            - commandLayout: the right part, contains all the buttons, fields, checkboxes...
        Both graphLayout and commandLayout are divided into sub-layouts.

        This function calls several functions to initialize each part of the window. The name of these functions contains
        'init_*layout'.
        """
        super(CalibWin, self).__init__(parent=parent)

        self.setWindowTitle("Energy Calibration")
        self.mainlayout = QHBoxLayout()
        self.graphlayout = QVBoxLayout()
        self.commandlayout = QVBoxLayout()
        self.commandlayout.setSpacing(10)

        self.init_var()

        # initialization of all the widgets/layouts
        self.init_btnlayout()
        self.init_tof2evlayout()
        self.init_eparlayout()
        self.init_fitparlayout()
        self.init_envectlayout()
        self.init_tofgraphlayout()
        self.init_graphauxlayout()

        # making the buttons not resizable
        for widget in self.children():
            if isinstance(widget, QPushButton):
                widget.setSizePolicy(0, 0)

        self.mainlayout.addLayout(self.graphlayout)
        self.mainlayout.addLayout(self.commandlayout)
        self.setLayout(self.mainlayout)
        self.show()

    ''' Initialization of instance attributes'''
    def init_var(self) -> None:

        self.withsb_bool = False
        self.calibloaded = False
        self.dataloaded = False
        self.bgndremoved = False
        self.threshyBool = False
        self.threshxminBool = False
        self.threshxmaxBool = False
        self.showexppeaksBool = False
        self.calibBool = False

    ''' In commandLayout - Initialization of the top right part of the layout, with 6 buttons (see just below)'''
    def init_btnlayout(self) -> None:

        btnlayout = QGridLayout()
        btnlayout.setSpacing(10)

        self.load_btn = QPushButton("Load data", self)
        self.rmbgnd_btn = QPushButton("Remove bgnd", self)
        self.findpeaks_btn = QPushButton("Find peaks", self)
        self.rmpeaks_btn = QPushButton("Remove peaks", self)
        self.exportcalib_btn = QPushButton("Export calib", self)
        self.importcalib_btn = QPushButton("Import calib", self)
        self.exportXUV_btn = QPushButton("Export XUV", self)

        self.rmbgnd_btn.setEnabled(False)
        self.findpeaks_btn.setEnabled(False)
        self.rmpeaks_btn.setEnabled(False)
        self.exportcalib_btn.setEnabled(False)
        self.importcalib_btn.setEnabled(False)
        self.exportXUV_btn.setEnabled(False)

        self.load_btn.clicked.connect(self.loadfile_lr)
        self.rmbgnd_btn.clicked.connect(self.rmbgnd_lr)
        self.findpeaks_btn.clicked.connect(self.findpeaks_lr)
        self.importcalib_btn.clicked.connect(self.importcalib_lr)
        self.rmpeaks_btn.clicked.connect(self.removepeaks_lr)
        self.exportXUV_btn.clicked.connect(self.exportXUV_lr)
        self.exportcalib_btn.clicked.connect(self.exportcalib_lr)

        btnlayout.addWidget(self.load_btn, 0, 0)
        btnlayout.addWidget(self.rmbgnd_btn, 0, 1)
        btnlayout.addWidget(self.findpeaks_btn, 1, 0)
        btnlayout.addWidget(self.rmpeaks_btn, 1, 1)
        btnlayout.addWidget(self.exportcalib_btn, 1, 3)
        btnlayout.addWidget(self.importcalib_btn, 1, 2)
        btnlayout.addWidget(self.exportXUV_btn, 0, 3)
        self.commandlayout.addLayout(btnlayout)

    ''' In commandLayout - Initialization of the tof to eV section: parameters of the af.find_local_maxima function, 
    'TOF to energy' button and 'with sidebands checkbox' '''
    def init_tof2evlayout(self) -> None:

        tof2evlayout = QHBoxLayout()
        flmlayout = QGridLayout()
        flmlayout.setSpacing(10)
        self.flm_box = QGroupBox(self)
        self.flm_box.setTitle("Find local maxima parameters")

        self.sm1_le = QLineEdit("5", self)
        self.sm2_le = QLineEdit("100", self)
        self.mindt_le = QLineEdit("10", self)

        flmlayout.addWidget(QLabel("smooth1"), 0, 0)
        flmlayout.addWidget(self.sm1_le, 1, 0)
        flmlayout.addWidget(QLabel("smooth2"), 0, 1)
        flmlayout.addWidget(self.sm2_le, 1, 1)
        flmlayout.addWidget(QLabel("min dt"), 0, 2)
        flmlayout.addWidget(self.mindt_le, 1, 2)
        self.flm_box.setLayout(flmlayout)

        for widget in self.flm_box.children():
            if isinstance(widget, QLineEdit):
                widget.setSizePolicy(0, 0)
                widget.setFixedSize(50, 20)

        self.tof2en_btn = QPushButton("TOF to energy", self)
        self.tof2en_btn.clicked.connect(self.tof2en_lr)
        self.tof2en_btn.setEnabled(False)

        self.withsb_cb = QCheckBox("With sidebands", self)
        self.withsb_cb.stateChanged.connect(self.withsb_fn)

        tof2evlayout.addWidget(self.flm_box)
        tof2evlayout.addWidget(self.withsb_cb)
        tof2evlayout.addWidget(self.tof2en_btn)
        self.commandlayout.addLayout(tof2evlayout)

    ''' In commandLayout - Initialization of the experimental parameters section: Retarding potential, TOF length, 
    wavelength, gas and first harmonic expected to see.'''
    def init_eparlayout(self) -> None:

        gases = cts.GASLIST

        epar_box = QGroupBox(self)
        epar_box.setTitle("Experimental parameters")
        epar_box.setSizePolicy(0, 0)
        eparlayout = QGridLayout()
        eparlayout.setSpacing(10)

        self.retpot_le = QLineEdit("0", self)
        self.toflength_le = QLineEdit("2", self)
        self.wvlength_le = QLineEdit("800", self)
        self.gas_combo = QComboBox(self)
        self.gas_combo.addItems(gases)
        self.firstharm_le = QLineEdit("13", self)

        self.retpot_le.returnPressed.connect(self.update_cts_fn)
        self.toflength_le.returnPressed.connect(self.update_cts_fn)
        self.wvlength_le.returnPressed.connect(self.update_cts_fn)
        self.firstharm_le.returnPressed.connect(self.update_cts_fn)
        self.gas_combo.currentIndexChanged.connect(self.gas_combo_lr)

        eparlayout.addWidget(QLabel("Ret. pot. (V)"), 0, 0)
        eparlayout.addWidget(self.retpot_le, 1, 0)
        eparlayout.addWidget(QLabel("TOF length (m)"), 0, 1)
        eparlayout.addWidget(self.toflength_le, 1, 1)
        eparlayout.addWidget(QLabel("lambda (nm)"), 0, 2)
        eparlayout.addWidget(self.wvlength_le, 1, 2)
        eparlayout.addWidget(QLabel("gas"), 0, 3)
        eparlayout.addWidget(self.gas_combo, 1, 3)
        eparlayout.addWidget(QLabel("1st harm."), 0, 4)
        eparlayout.addWidget(self.firstharm_le, 1, 4)

        epar_box.setLayout(eparlayout)

        for widget in epar_box.children():
            if isinstance(widget, QLineEdit) or isinstance(widget, QComboBox):
                widget.setSizePolicy(0, 0)
                widget.setFixedSize(50, 20)

        self.commandlayout.addWidget(epar_box)

    ''' In commandLayout - Initialization of the fit parameter section with a, t0 and c. First line: guess values calculated 
    from the experimental parameters. Second line: fitted values.'''
    def init_fitparlayout(self) -> None:

        fitpar_box = QGroupBox()
        fitpar_box.setTitle("Calibration parameters")
        fitpar_box.setSizePolicy(0, 0)
        fitparlayout = QGridLayout()
        fitparlayout.setSpacing(10)

        self.aguess_lb = QLabel(self)
        self.afit_lb = QLabel(self)
        self.t0guess_lb = QLabel(self)
        self.t0fit_lb = QLabel(self)
        self.cguess_lb = QLabel(self)
        self.cfit_lb = QLabel(self)

        fitparlayout.addWidget(QLabel("a"), 0, 0)
        fitparlayout.addWidget(self.aguess_lb, 1, 0)
        fitparlayout.addWidget(self.afit_lb, 2, 0)
        fitparlayout.addWidget(QLabel("t0"), 0, 1)
        fitparlayout.addWidget(self.t0guess_lb, 1, 1)
        fitparlayout.addWidget(self.t0fit_lb, 2, 1)
        fitparlayout.addWidget(QLabel("c"), 0, 2)
        fitparlayout.addWidget(self.cguess_lb, 1, 2)
        fitparlayout.addWidget(self.cfit_lb, 2, 2)
        fitparlayout.addWidget(QLabel("guess"), 1, 3)
        fitparlayout.addWidget(QLabel("calc"), 2, 3)
        text = "q=a*1/(t-t0)²+c"
        fitparlayout.addWidget(QLabel(text), 1, 4)

        fitpar_box.setLayout(fitparlayout)

        for widget in fitpar_box.children():
            if isinstance(widget, QLabel) and widget.text() != text:
                widget.setSizePolicy(0, 0)
                widget.setFixedSize(55, 15)

        self.commandlayout.addWidget(fitpar_box)
        self.gas_combo.setCurrentIndex(2)  # Argon
        self.update_cts_fn()

    ''' In commandLayout - Initialization of the resulting energy vector section, with elow, ehigh and dE'''
    def init_envectlayout(self) -> None:

        cts.elow = (float(self.firstharm_le.text()) - 1) * cts.HEV * cts.cur_nu
        cts.ehigh = float(34 * cts.HEV * cts.cur_nu)
        cts.dE = 0.01

        envect_box = QGroupBox()
        envect_box.setTitle("Energy vector parameters")
        envect_box.setSizePolicy(0, 0)
        envectlayout = QGridLayout()
        envectlayout.setSpacing(10)

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

        self.commandlayout.addWidget(envect_box)
        self.update_envect_fn()

    ''' In graphLayout - Initialization of the top figure on the window, where the time of flight is plotted'''
    def init_tofgraphlayout(self) -> None:

        tof_fig = Figure(figsize=(4, 3), dpi=100)
        self.tof_fc = FigureCanvas(tof_fig)
        self.tof_fc.mpl_connect('button_press_event', self.onclick)
        self.tof_ax = self.tof_fc.figure.add_subplot(111)
        self.tof_fc.draw()
        tof_nav = NavigationToolbar2QT(self.tof_fc, self)
        tof_nav.setStyleSheet("QToolBar { border: 0px }")

        tgparalayout = QGridLayout()
        self.threshtype_bgroup = QButtonGroup()
        rblayout = QGridLayout()
        self.setth_cb = QCheckBox("Set threshold", self)
        self.setth_cb.setEnabled(False)
        self.setth_cb.stateChanged.connect(self.setth_lr)
        self.addpeak_cb = QCheckBox("Add peak", self)
        self.addpeak_cb.setEnabled(False)
        self.addpeak_cb.stateChanged.connect(self.addpeak_lr)
        self.showexppeaks_cb = QCheckBox("Show expected peaks", self)
        self.showexppeaks_cb.stateChanged.connect(self.showexppeaks_lr)
        self.showexppeaks_cb.setEnabled(False)

        self.clear_btn = QPushButton("Clear", self)
        self.clear_btn.clicked.connect(self.clear_lr)
        self.clear_btn.setEnabled(False)

        self.Y_rb = QRadioButton("Y", self)
        self.Y_rb.value = "Y"
        self.Y_rb.toggled.connect(self.threshtype_lr)
        self.Y_rb.toggle()
        self.Y_rb.setEnabled(False)

        self.xmin_rb = QRadioButton("X min", self)
        self.xmin_rb.value = "Xm"
        self.xmin_rb.toggled.connect(self.threshtype_lr)
        self.xmin_rb.setEnabled(False)

        self.xmax_rb = QRadioButton("X max", self)
        self.xmax_rb.value = "XM"
        self.xmax_rb.toggled.connect(self.threshtype_lr)
        self.xmax_rb.setEnabled(False)

        self.threshtype_bgroup.addButton(self.Y_rb)
        self.threshtype_bgroup.addButton(self.xmin_rb)
        self.threshtype_bgroup.addButton(self.xmax_rb)

        rblayout.addWidget(self.Y_rb, 0, 0)
        rblayout.addWidget(self.xmin_rb, 0, 1)
        rblayout.addWidget(self.xmax_rb, 0, 2)
        tgparalayout.addWidget(self.setth_cb, 0, 0)
        tgparalayout.addWidget(self.addpeak_cb, 0, 1)
        tgparalayout.addWidget(self.showexppeaks_cb, 0, 2)
        tgparalayout.addWidget(self.clear_btn, 0, 3)
        tgparalayout.addLayout(rblayout, 1, 0, 1, 3)

        self.graphlayout.addWidget(self.tof_fc)
        self.graphlayout.addWidget(tof_nav)
        self.graphlayout.addLayout(tgparalayout)

    ''' In graphLayout - Initialization the two bottom figures on the window'''
    def init_graphauxlayout(self) -> None:

        graphauxlayout = QHBoxLayout()
        ga1layout = QVBoxLayout()
        ga2layout = QVBoxLayout()

        fit_fig = Figure(figsize=(2, 2), dpi=100)
        self.fit_fc = FigureCanvas(fit_fig)
        self.fit_fc.setSizePolicy(1, 0)
        self.fit_ax = self.fit_fc.figure.add_subplot(111)
        self.fit_ax.tick_params(labelsize = 8)
        self.fit_fc.draw()
        ga1layout.addWidget(self.fit_fc)

        en_fig = Figure(figsize=(2, 2), dpi=100)
        self.en_fc = FigureCanvas(en_fig)
        self.en_fc.setSizePolicy(1, 0)
        self.en_ax = self.en_fc.figure.add_subplot(111)
        self.en_ax.tick_params(labelsize = 8)
        self.en_fc.draw()
        en_nav = NavigationToolbar2QT(self.en_fc, self)
        en_nav.setStyleSheet("QToolBar { border: 0px }")
        ga2layout.addWidget(self.en_fc)
        ga2layout.addWidget(en_nav)

        graphauxlayout.addLayout(ga1layout)
        graphauxlayout.addLayout(ga2layout)
        self.graphlayout.addLayout(graphauxlayout)

    ''' From a calibration file, loading afit, t0fit, cfit and the first harmonic'''
    def importcalib_lr(self) -> None:
        calib_fname = QFileDialog.getOpenFileName(self, 'Import calibration')
        calib_f = calib_fname[0]
        if(calib_f):
            with open(calib_f, 'r') as file:
                f = file.read().splitlines()
                try:
                    if(len(f) < 4):
                        raise IndexError
                    cts.afit, cts.t0fit, cts.cfit, cts.first_harm = [float(f[i]) for i in range(4)]
                    cts.first_harm = int(cts.first_harm)
                    self.calibloaded = True
                    self.calibBool = True
                    self.firstharm_le.setText(str(cts.first_harm))
                    self.update_fitpar_fn()
                    #print(repr(self.window()))
                    self.window().updateglobvar_fn()
                    self.tof2en_btn.setEnabled(True)

                except IndexError:
                    print('Not enough data in the calib file. Needed: a_fit, t0_fit, c_fit and first harmonic')
                except ValueError:
                    print("Incorrect calibration data")

    ''' Exporting afit, t0fit, cfit and the first harmonic in a file after choosing its name and location'''
    def exportcalib_lr(self) -> None:
        filename = QFileDialog.getSaveFileName(self, 'Save XUV')
        fname = filename[0]
        if fname:
            fit_param = np.append(self.p_opt, cts.first_harm)
            print('ok')
            np.savetxt(fname, fit_param, fmt='%1.4e', delimiter='\t')

    ''' Updating the fit parameters on the window'''
    def update_fitpar_fn(self) -> None:

        if (self.calibBool):
            self.afit_lb.setText("{:.3e}".format(cts.afit))
            self.t0fit_lb.setText("{:.3e}".format(cts.t0fit))
            self.cfit_lb.setText("{:.3f}".format(cts.cfit))

    ''' Updates the experimental parameters (and aguess, t0guess, cguess)'''
    def update_cts_fn(self) -> None:

        try:
            cts.cur_nu = cts.C / (float(self.wvlength_le.text()) * 1e-9)
            cts.cur_Vp = float(self.retpot_le.text())
            cts.cur_L = float(self.toflength_le.text())
            cts.first_harm = int(self.firstharm_le.text())
            self.window().updateglobvar_fn()
            self.aguess = (0.5 * cts.ME * cts.cur_L ** 2 / cts.QE) / (cts.HEV * cts.cur_nu)
            self.t0guess = 5.8e-8
            self.cguess = (cts.cur_Ip + cts.cur_Vp) / (cts.HEV * cts.cur_nu)
            self.aguess_lb.setText("{:.3e}".format(self.aguess))
            self.t0guess_lb.setText("{:.3e}".format(self.t0guess))
            self.cguess_lb.setText("{:.3f}".format(self.cguess))
        except ValueError:
            print('Incorrect value. Must be a number')
        except ZeroDivisionError:
            print('ZeroDivisionError')

    ''' Updating th energy vector parameters with the values written in the associated QLineEdit objects'''
    def update_envect_fn(self) -> None:

        cts.elow = float(self.elow_le.text())
        cts.ehigh = float(self.ehigh_le.text())
        cts.dE = float(self.dE_le.text())
        self.elow_le.setText("{:.2f}".format(cts.elow))
        self.ehigh_le.setText("{:.2f}".format(cts.ehigh))
        self.window().updateglobvar_fn()

    ''' Gas QCombobox listener'''
    def gas_combo_lr(self, i) -> None:
        cts.cur_Ip = cts.IPLIST[i]
        self.update_cts_fn()

    ''' Listener of the threshold radiobuttons: Y, Xmin or Xmax'''
    def threshtype_lr(self) -> None:

        rb = self.sender()
        self.threshtype = rb.value

    ''' "with sidebands" checkbox listener '''
    def withsb_fn(self, state) -> None:
        if state == Qt.Checked:
            self.withsb_bool = True
        else:
            self.withsb_bool = False

    ''' "load data" button listener'''
    def loadfile_lr(self) -> None:

        TOFfile = QFileDialog.getOpenFileName(self, 'Load data')
        fname = TOFfile[0]
        if(fname):
            with open(fname, 'r') as file:
                try:
                    data = [[float(digit) for digit in line.split()] for line in file]
                    data_array = np.asarray(data)  # converting 1D list into 2D array
                    self.counts = data_array

                    self.thxmin = 0
                    self.thy = 0
                    self.thxmax = self.counts.shape[0] * 1e-9

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
                    self.importcalib_btn.setEnabled(True)
                    self.withsb_cb.setEnabled(True)

                    self.dataloaded = True
                    self.threshxminBool = False
                    self.threshxmaxBool = False
                    self.threshyBool = False
                    self.peaksfound = False

                except ValueError:
                    print('Incorrect data file')
                except IndexError:
                    print('Incorrect data file')

    ''' "remove bgnd" button listener'''
    def rmbgnd_lr(self) -> None:

        self.counts[:,1] = self.counts[:,1] - self.counts[300:600,1].mean()

        self.bgndremoved = True
        self.peaksfound = False
        self.threshyBool = False
        self.threshxminBool = False
        self.threshxmaxBool = False
        self.thxmin = 0
        self.thy = 0
        self.thxmax = self.counts.shape[0] * 1e-9

        self.rmpeaks_btn.setEnabled(False)
        self.tof2en_btn.setEnabled(False)
        self.rmbgnd_btn.setEnabled(False)
        self.refreshplot_fn()

    ''' "show expected peaks" checkbox listener '''
    def showexppeaks_lr(self, state) -> None:
        if state == Qt.Checked:
            self.showexppeaksBool = True
        else:
            self.showexppeaksBool = False
        self.refreshplot_fn()

    ''' "set threshold" checkbox listener '''
    def setth_lr(self) -> None:
        if self.setth_cb.isChecked():
            self.addpeak_cb.setCheckState(Qt.Unchecked)

    ''' "add peak" checkbox listener '''
    def addpeak_lr(self) -> None:
        if self.addpeak_cb.isChecked():
            self.setth_cb.setCheckState(Qt.Unchecked)

    ''' "remove peaks" button listener '''
    def removepeaks_lr(self) -> None:
        rmp = rmPeaksDialog(self) # new class defined below

    ''' "find peaks" button listener '''
    def findpeaks_lr(self) -> None:

        try:
            ip, dp, convolution = af.find_local_maxima(self.counts[:, 1], self.thy, self.thxmin, self.thxmax,
                                                       int(self.sm1_le.text()), int(self.sm2_le.text()),
                                                       int(self.mindt_le.text()))
            self.maximaIndices = (ip*1e-9).tolist()
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

    ''' ANALYSIS - "TOF to energy" button listener '''
    def tof2en_lr(self) -> None:
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
                qq = np.arange(qqi, qqi+len(tofHHG),1)
            else:
                qq = np.arange(qqi, qqi + 2*len(tofHHG), 2)
            qq = qq[::-1]

            self.p_opt, self.pcov = opt.curve_fit(af.tof2EeV,tofHHG,qq,np.array([self.aguess,self.t0guess,self.cguess]))
            toffit = np.linspace(tofHHG[0], tofHHG[-1], 100)
            self.fit_ax.plot(toffit, af.tof2EeV(toffit,self.p_opt[0],self.p_opt[1],self.p_opt[2]),label='Fit')
            self.fit_ax.plot(toffit, af.tof2EeV(toffit,self.aguess,self.t0guess,self.cguess),label='Guess')
            self.fit_ax.plot(tofHHG, qq, 'ys',label='TOF(peaks)')
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

    ''' Updating the top left (TOF) graph'''
    def refreshplot_fn(self) -> None:
        xmin, xmax = self.tof_ax.get_xlim()
        ymin, ymax = self.tof_ax.get_ylim()
        self.tof_ax.cla()

        self.tof_ax.xaxis.set_major_formatter(FormatStrFormatter('%2.e'))
        self.tof_ax.yaxis.set_major_formatter(FormatStrFormatter('%1.e'))
        self.tof_ax.set_ylabel("counts (arb. units)")
        self.tof_ax.set_xlabel("TOF (s)")

        if self.dataloaded:
            self.tof_ax.plot(self.counts[:, 0], self.counts[:, 1])
        if self.threshyBool:
            self.tof_ax.plot(self.threshyline[:, 0], self.threshyline[:, 1], 'k')
        if self.threshxminBool:
            self.tof_ax.plot(self.threshxminline[:, 0], self.threshxminline[:, 1], 'k')
        if self.threshxmaxBool:
            self.tof_ax.plot(self.threshxmaxline[:, 0], self.threshxmaxline[:, 1], 'k')
        if self.peaksfound:
            self.tof_ax.plot(self.maximaIndices, self.maximaIntensity, 'ro')
            self.tof_ax.plot(self.counts[:, 0], self.convolvedsignal)

        if self.showexppeaksBool:
            qq2 = np.arange(cts.first_harm, 35, 1)

            y = np.linspace(ymax - (ymax - ymin) * 0.2, ymax, 100)
            for i in range((len(qq2))):
                if qq2[i] % 2 == 0:
                    c = 'r'
                else:
                    c = 'k'
                xval = float(np.math.sqrt(
                    0.5 * cts.ME * cts.cur_L ** 2 / cts.QE / (qq2[i] * cts.HEV * cts.cur_nu - cts.cur_Ip)) + 6e-8)
                x = np.full((100, 1), xval)
                self.tof_ax.plot(x, y, color=c, linewidth=1.0)

        if self.bgndremoved:
            self.bgndremoved = False #this means that when we remove bgnd, we don't keep the same scale
        else:
            self.tof_ax.set_ylim(ymin, ymax)

        self.tof_ax.set_xlim(xmin, xmax)
        self.tof_fc.draw()

    ''' called when double-clicking on the TOF graph'''
    def onclick(self, event) -> None:
        if self.dataloaded and self.addpeak_cb.isChecked(): #add a peak on clicking on the figure
            i=0
            ifound = False
            while(i<len(self.maximaIndices) and ifound == False):
                if(self.maximaIndices[i] < event.xdata):
                    i += 1
                else:
                    ifound = True
            self.maximaIndices.insert(i, event.xdata)
            self.maximaIntensity.insert(i, event.ydata)
            self.refreshplot_fn()

        if self.dataloaded and self.setth_cb.isChecked():
            if self.Y_rb.isChecked():
                self.thy = event.ydata
                self.threshyline = np.full((len(self.counts), 2), self.thy)
                self.threshyline[:, 0] = self.counts[:, 0]
                self.threshyBool = True
            elif self.xmin_rb.isChecked():
                self.thxmin = event.xdata
                self.threshxminline = np.full((len(self.counts), 2), self.thxmin)
                self.threshxminline[:, 1] = self.counts[:, 1]
                self.threshxminBool = True
            elif self.xmax_rb.isChecked():
                self.thxmax = event.xdata
                self.threshxmaxline = np.full((len(self.counts), 2), self.thxmax)
                self.threshxmaxline[:, 1] = self.counts[:, 1]
                self.threshxmaxBool = True
            self.refreshplot_fn()

    ''' "Export XUV" button listener. Saving the energy vector and the energy-converted TOF signal'''
    def exportXUV_lr(self) -> None:
        filename = QFileDialog.getSaveFileName(self, 'Save XUV')
        fname = filename[0]
        if fname:
            XUV_array = np.vstack((self.Eevlin, self.signal)).T
            np.savetxt(fname, XUV_array, delimiter='\t')

    ''' "Clear" button listener. Resets all the objects, but not the global variables'''
    def clear_lr(self) -> None:
        self.tof_ax.cla()
        self.fit_ax.cla()
        self.en_ax.cla()

        self.dataloaded = False
        self.threshyBool = False
        self.threshxminBool = False
        self.threshxmaxBool = False
        self.calibloaded = False
        self.calibBool = False

        self.Y_rb.setEnabled(False)
        self.xmin_rb.setEnabled(False)
        self.xmax_rb.setEnabled(False)

        for w in self.children():
            if isinstance(w, QPushButton) or isinstance(w, QCheckBox):
                w.setEnabled(False)
        self.load_btn.setEnabled(True)
        self.importcalib_btn.setEnabled(False)

        self.tof_fc.draw()
        self.fit_fc.draw()
        self.en_fc.draw()

        cts.clear_varlist()
        self.window().updateglobvar_fn()


''' created when clicking on the "remove peaks" button'''
class rmPeaksDialog(QDialog):
    def __init__(self, parent=CalibWin) -> None:

        super(rmPeaksDialog, self).__init__(parent)

        self.par = self.parent()
        self.currentindex = 0

        self.setGeometry(850, 500, 450, 300)
        self.setWindowTitle("Remove peaks")
        self.mainlayout = QHBoxLayout()
        self.rightlayout = QVBoxLayout()
        self.thirdlayout = QHBoxLayout()

        self.peaksTable = QTableWidget()
        self.peaksTable.setSelectionBehavior(QTableView.SelectRows)
        self.peaksTable.clicked.connect(self.changeindex_lr)
        self.peaksTable.verticalHeader().sectionClicked.connect(self.changeindex2_lr)

        self.refreshTable_fn()

        self.rmpeak_btn = QPushButton("Remove", self)
        self.ok_btn = QPushButton("Done", self)
        self.cancel_btn = QPushButton("Cancel", self)

        self.rmpeak_btn.clicked.connect(self.rmpeak_lr)
        self.ok_btn.clicked.connect(self.ok_lr)
        self.cancel_btn.clicked.connect(self.cancel_lr)

        self.thirdlayout.addWidget(self.ok_btn)
        self.thirdlayout.addWidget(self.cancel_btn)
        self.rightlayout.addWidget(self.rmpeak_btn)
        self.rightlayout.addLayout(self.thirdlayout)
        self.mainlayout.addWidget(self.peaksTable)
        self.mainlayout.addLayout(self.rightlayout)
        self.setLayout(self.mainlayout)

        self.show()

    ''' Updating the table'''
    def refreshTable_fn(self) -> None:
        self.peaksTable.clear()
        self.peaksTable.setRowCount(len(self.par.maximaIndices) + 1)
        self.peaksTable.setColumnCount(2)

        self.peaksTable.setItem(0, 0, QTableWidgetItem("xcoord"))
        self.peaksTable.setItem(0, 1, QTableWidgetItem("ycoord"))

        for i in range(len(self.par.maximaIndices)):
            self.peaksTable.setItem(i + 1, 0, QTableWidgetItem("{:.3e}".format(self.par.maximaIndices[i])))
            self.peaksTable.setItem(i + 1, 1, QTableWidgetItem("{:.3e}".format(self.par.maximaIntensity[i])))

    ''' Updating the row index on clicking on the vertical header'''
    def changeindex2_lr(self, logicalIndex) -> None:
        self.currentindex = logicalIndex - 1

    ''' Updating the row index on clicking anywhere in the row'''
    def changeindex_lr(self, clickedIndex) -> None:
        self.currentindex = clickedIndex.row()-1

    ''' "Remove peaks" button listener '''
    def rmpeak_lr(self) -> None:
        del self.par.maximaIndices[self.currentindex]
        del self.par.maximaIntensity[self.currentindex]
        self.par.refreshplot_fn()
        self.refreshTable_fn()
        if(self.currentindex>len(self.par.maximaIndices)-1):
            self.peaksTable.selectRow(self.currentindex)
            self.currentindex -= 1
        else:
            self.peaksTable.selectRow(self.currentindex+1)

    ''' "ok" button listener '''
    def ok_lr(self) -> None:
        self.destroy()

    ''' "cancel" button listener. Calling the parent function "findpeaks_lr" erases all the changes '''
    def cancel_lr(self) -> None:
        self.par.findpeaks_lr()
        self.par.refreshplot_fn()
        self.destroy()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = CalibWin()
    ex.setGeometry(300, 300, 1000, 600)
    sys.exit(app.exec_())