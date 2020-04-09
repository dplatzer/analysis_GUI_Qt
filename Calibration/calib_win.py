import sys, traceback
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QPushButton, QApplication, QGridLayout, QHBoxLayout, QVBoxLayout, QGroupBox, QLabel
from PyQt5.QtWidgets import QLineEdit, QCheckBox, QComboBox, QRadioButton, QButtonGroup, QFileDialog, QDialog, QTableWidget
from PyQt5.QtWidgets import QTableWidgetItem, QTableView, QFrame, QSplitter
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np

from matplotlib.ticker import FormatStrFormatter

# Homemade modules
import glob_var as cts
from Calibration import _calib_functions as _calf
import _import_export as _ie

class CalibWin(QWidget, _calf.calib_functions_mixin, _ie.Imp_Exp_Mixin):
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
    ##################################################################################
    ############################ Widget Initialization ###############################

    def __init__(self, parent=None):
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

        # initialization of all the widgets/layouts
        self.init_btnlayout()
        self.init_tof2evlayout()
        self.init_eparlayout()
        self.init_fitparlayout()
        self.init_envectlayout()
        self.init_exportlayout()
        self.init_tofgraphlayout()
        self.init_graphauxlayout()

        # making the buttons not resizable
        for widget in self.children():
            if isinstance(widget, QPushButton):
                widget.setSizePolicy(0, 0)

        self.mainlayout.addWidget(self.splitter2)
        self.mainlayout.addLayout(self.commandlayout)
        self.setLayout(self.mainlayout)

        self.init_var()

        self.show()

    def init_var(self):
        ''' Initialization of instance attributes'''

        self.withsb_bool = False
        self.calibloaded = False
        self.dataloaded = False
        self.bgndremoved = False
        self.threshyBool = False
        self.threshxminBool = False
        self.threshxmaxBool = False
        self.showexppeaksBool = False
        self.calibBool = False
        self.peaksfound = False

        self.thxmin = 0
        self.thy = 0
        self.thxmax = 0

        self.counts = []

        self.gas_combo.setCurrentIndex(cts.cur_gas_index)
        self.software_combo.setCurrentIndex(cts.software_version)

    def init_btnlayout(self):
        ''' In commandLayout - Initialization of the top right part of the layout, with 6 buttons (see just below)'''

        btnlayout = QGridLayout()
        btnlayout.setSpacing(10)
        btn_gb = QGroupBox()

        self.load_btn = QPushButton("Load data", self)
        self.rmbgnd_btn = QPushButton("Remove bgnd", self)
        self.findpeaks_btn = QPushButton("Find peaks", self)
        self.rmpeaks_btn = QPushButton("Remove peaks", self)

        self.scan_i_le = QLineEdit(str(cts.scan_i), self)
        self.spectrum_i_le = QLineEdit(str(cts.spectrum_i), self)
        self.skiplines_le = QLineEdit(str(cts.skiplines), self)
        self.filenametype_le = QLineEdit(str(cts.filenametype), self)

        # remove background average
        rmbg_avg_layout = QGridLayout()
        rmbg_avg_layout.setSpacing(10)
        self.rmbg_avg_box = QGroupBox(self)
        self.rmbg_avg_box.setTitle("Remove background")
        self.rmbg_avg_box.setFixedHeight(90)
        self.rmbg_avg_box.setSizePolicy(0, 0)

        self.rmbg_avg_min_le = QLineEdit(str(cts.rmbg_avg_min), self)
        self.rmbg_avg_max_le = QLineEdit(str(cts.rmbg_avg_max), self)

        rmbg_avg_layout.addWidget(QLabel("average between"), 0, 0)
        rmbg_avg_layout.addWidget(self.rmbg_avg_min_le, 0, 1)
        rmbg_avg_layout.addWidget(QLabel("and"), 0, 2)
        rmbg_avg_layout.addWidget(self.rmbg_avg_max_le, 0, 3)
        rmbg_avg_layout.addWidget(self.rmbgnd_btn, 1, 0)
        self.rmbg_avg_box.setLayout(rmbg_avg_layout)

        self.rmbg_avg_min_le.setSizePolicy(0, 0)
        self.rmbg_avg_min_le.setFixedSize(50, 20)
        self.rmbg_avg_max_le.setSizePolicy(0, 0)
        self.rmbg_avg_max_le.setFixedSize(50, 20)
        self.scan_i_le.setSizePolicy(0, 0)
        self.scan_i_le.setFixedSize(40, 20)
        self.spectrum_i_le.setSizePolicy(0, 0)
        self.spectrum_i_le.setFixedSize(40, 20)
        self.skiplines_le.setSizePolicy(0, 0)
        self.skiplines_le.setFixedSize(40, 20)

        self.rmbgnd_btn.setEnabled(False)
        self.findpeaks_btn.setEnabled(False)
        self.rmpeaks_btn.setEnabled(False)

        if cts.software_version != 2:
            self.spectrum_i_le.setEnabled(False)
            self.scan_i_le.setEnabled(False)

        self.load_btn.clicked.connect(self.loadfile_lr)
        self.rmbgnd_btn.clicked.connect(self.rmbgnd_lr)
        self.findpeaks_btn.clicked.connect(self.findpeaks_lr)
        self.rmpeaks_btn.clicked.connect(self.removepeaks_lr)

        self.rmbg_avg_min_le.returnPressed.connect(self.update_rmbg_avg_fn)
        self.rmbg_avg_max_le.returnPressed.connect(self.update_rmbg_avg_fn)
        self.spectrum_i_le.returnPressed.connect(self.update_spectrum_i)
        self.scan_i_le.returnPressed.connect(self.update_scan_i)
        self.skiplines_le.returnPressed.connect(self.update_skiplines)
        self.filenametype_le.returnPressed.connect(self.update_filenametype)

        self.software_combo = QComboBox()
        self.software_combo.addItems(cts.SOFTWARE_LIST)
        self.software_combo.currentIndexChanged.connect(self.software_combo_lr)

        btnlayout.addWidget(self.load_btn, 0, 0)
        btnlayout.addWidget(QLabel("Scan"), 1, 0)
        btnlayout.addWidget(self.scan_i_le, 1, 1)
        btnlayout.addWidget(QLabel("Spectrum"), 2, 0)
        btnlayout.addWidget(self.spectrum_i_le, 2, 1)
        btnlayout.addWidget(QLabel("Skip lines"), 1, 2)
        btnlayout.addWidget(self.skiplines_le, 2, 2)
        btnlayout.addWidget(QLabel("Software"), 1, 3)
        btnlayout.addWidget(self.software_combo, 2, 3)
        btnlayout.addWidget(QLabel("File name type"), 3, 2, 1, 2)
        btnlayout.addWidget(self.filenametype_le, 4, 2, 1, 2)
        btnlayout.addWidget(self.findpeaks_btn, 0, 2)
        btnlayout.addWidget(self.rmpeaks_btn, 0, 3)
        btnlayout.addWidget(self.rmbg_avg_box, 5, 0, 1, 4)

        btn_gb.setLayout(btnlayout)
        btn_gb.setSizePolicy(0, 0)

        self.commandlayout.addWidget(btn_gb)

    def init_tof2evlayout(self):
        ''' In commandLayout - Initialization of the tof to eV section: parameters of the af.find_local_maxima function,
            'TOF to energy' button and 'with sidebands checkbox' '''

        tof2evlayout = QHBoxLayout()
        tof2ev_box = QGroupBox()
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
        self.flm_box.setSizePolicy(0, 0)

        self.tof2en_btn = QPushButton("TOF to energy", self)
        self.tof2en_btn.clicked.connect(self.tof2en_lr)
        self.tof2en_btn.setEnabled(False)

        self.withsb_cb = QCheckBox("With sidebands", self)
        self.withsb_cb.setSizePolicy(0, 0)
        self.withsb_cb.stateChanged.connect(self.withsb_fn)

        tof2evlayout.addWidget(self.flm_box)
        tof2evlayout.addWidget(self.withsb_cb)
        tof2evlayout.addWidget(self.tof2en_btn)
        tof2ev_box.setLayout(tof2evlayout)
        tof2ev_box.setSizePolicy(0, 0)
        self.commandlayout.addWidget(tof2ev_box)

    def init_eparlayout(self):
        ''' In commandLayout - Initialization of the experimental parameters section: Retarding potential, TOF length,
            wavelength, gas and first harmonic expected to see.'''

        epar_box = QGroupBox(self)
        epar_box.setTitle("Experimental parameters")
        epar_box.setSizePolicy(0, 0)
        eparlayout = QGridLayout()
        eparlayout.setSpacing(10)

        self.retpot_le = QLineEdit(str(cts.cur_Vp), self)
        self.toflength_le = QLineEdit(str(cts.cur_L), self)
        self.wvlength_le = QLineEdit(str(cts.lambda_start), self)
        self.gas_combo = QComboBox(self)
        self.gas_combo.addItems(cts.GASLIST)
        self.firstharm_le = QLineEdit(str(cts.first_harm), self)

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
                widget.setFixedSize(50, 20)
                widget.setSizePolicy(0, 0)


        self.commandlayout.addWidget(epar_box)

    def init_fitparlayout(self):
        ''' In commandLayout - Initialization of the fit parameter section with a, t0 and c. First line: guess values calculated
            from the experimental parameters. Second line: fitted values.'''

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

    def init_envectlayout(self):
        ''' In commandLayout - Initialization of the resulting energy vector section, with elow, ehigh and dE'''

        self.envect_export_layout = QHBoxLayout()  # for both energy vector and export boxes
        self.envect_export_box = QGroupBox()

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

        self.envect_export_layout.addWidget(envect_box)
        self.update_envect_fn()

    def init_exportlayout(self):
        exportlayout = QHBoxLayout()

        self.exportcalib_btn = QPushButton("Export calib", self)
        self.exportXUV_btn = QPushButton("Export XUV", self)

        self.exportcalib_btn.setEnabled(False)
        self.exportXUV_btn.setEnabled(False)

        self.exportXUV_btn.clicked.connect(self.exportXUV_lr)
        self.exportcalib_btn.clicked.connect(self.exportcalib_lr)

        exportlayout.addWidget(self.exportcalib_btn)
        exportlayout.addWidget(self.exportXUV_btn)

        self.envect_export_layout.addLayout(exportlayout)
        self.envect_export_box.setLayout(self.envect_export_layout)
        self.envect_export_box.setSizePolicy(0, 0)
        self.commandlayout.addWidget(self.envect_export_box)

        empty = QWidget()
        empty.setSizePolicy(0, 1)
        self.commandlayout.addWidget(empty)

    def init_tofgraphlayout(self):
        ''' In graphLayout - Initialization of the top figure on the window, where the time of flight is plotted'''
        self.splitter1 = QSplitter(Qt.Horizontal)
        self.splitter2 = QSplitter(Qt.Vertical)

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

        self.minus_sign_cb = QCheckBox("*(-1)", self)
        self.minus_sign_cb.setEnabled(False)
        self.minus_sign_cb.stateChanged.connect(self.minus_sign_lr)

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
        tgparalayout.addWidget(self.minus_sign_cb, 1, 3)

        self.graphlayout.addWidget(self.tof_fc)
        self.graphlayout.addWidget(tof_nav)
        self.graphlayout.addLayout(tgparalayout)
        self.graph_frame = QFrame()
        self.graph_frame.setLayout(self.graphlayout)

    def init_graphauxlayout(self):
        ''' In graphLayout - Initialization the two bottom figures on the window'''

        graphauxlayout = QHBoxLayout()
        ga1layout = QVBoxLayout()
        self.ga1_frame = QFrame()
        self.ga1_frame.setLayout(ga1layout)
        ga2layout = QVBoxLayout()
        self.ga2_frame = QFrame()
        self.ga2_frame.setLayout(ga2layout)

        fit_fig = Figure(figsize=(2, 2), dpi=100)
        self.fit_fc = FigureCanvas(fit_fig)
        self.fit_fc.setSizePolicy(1, 0)
        self.fit_ax = self.fit_fc.figure.add_subplot(111)
        self.fit_ax.tick_params(labelsize = 8)
        ga1layout.addWidget(self.fit_fc)
        self.fit_fc.draw()

        en_fig = Figure(figsize=(2, 2), dpi=100)
        self.en_fc = FigureCanvas(en_fig)
        #self.en_fc.setSizePolicy(1, 0)
        self.en_ax = self.en_fc.figure.add_subplot(111)
        self.en_ax.tick_params(labelsize = 8)
        self.en_fc.draw()
        en_nav = NavigationToolbar2QT(self.en_fc, self)
        en_nav.setStyleSheet("QToolBar { border: 0px }")
        ga2layout.addWidget(self.en_fc)
        ga2layout.addWidget(en_nav)

        self.splitter1.addWidget(self.ga1_frame)
        self.splitter1.addWidget(self.ga2_frame)

        self.splitter2.addWidget(self.graph_frame)
        self.splitter2.addWidget(self.splitter1)

    #################################################################################
    ############################ Other methods ######################################

    def software_combo_lr(self, i):
        # Labview 2013
        if i == 0:
            cts.minus_sign = False
            cts.skiplines = 0  # number of lines in the header of data files. These lines are skipped when loading data
            cts.filenametype = '*delay_0*'
            self.scan_i_le.setEnabled(False)
            self.spectrum_i_le.setEnabled(False)

        # Labview 2016
        elif i == 1:
            cts.minus_sign = False
            cts.skiplines = 1200  # number of lines in the header of data files. These lines are skipped when loading data
            cts.filenametype = '*_PE_c00*.LCY'
            self.scan_i_le.setEnabled(False)
            self.spectrum_i_le.setEnabled(False)

        # PyMoDAQ
        elif i == 2:
            cts.minus_sign = False
            cts.skiplines = 0
            cts.filenametype = 'Not needed'
            self.scan_i_le.setEnabled(True)
            self.spectrum_i_le.setEnabled(True)

        else:
            print('Incorrect Labview version value')
            raise ValueError
        cts.software_version = i
        self.skiplines_le.setText(str(cts.skiplines))
        self.filenametype_le.setText(str(cts.filenametype))
        self.window().updateglobvar_fn()

    def update_filenametype(self):
        cts.filenametype = self.filenametype_le.text()

    def update_skiplines(self):
        try:
            cts.skiplines = int(self.skiplines_le.text())
        except ValueError:
            pass
        finally:
            self.skiplines_le.setText(str(cts.skiplines))
            self.window().updateglobvar_fn()

    def update_scan_i(self):
        try:
            cts.scan_i = int(self.scan_i_le.text())
            if cts.scan_i < 10:
                cts.str_scan_i = "".join(["00", str(cts.scan_i)])
            # I assume that we do less than 99 scans a day!
            else:
                cts.str_scan_i = "".join(["0", str(cts.scan_i)])
        except ValueError:
            cts.scan_i = 0
            cts.str_scan_i = '000'
        finally:
            self.scan_i_le.setText(str(cts.scan_i))
            self.window().updateglobvar_fn()

    def update_spectrum_i(self):
        try:
            cts.spectrum_i = int(self.spectrum_i_le.text())
        except ValueError:
            cts.spectrum_i = 0
        finally:
            self.spectrum_i_le.setText(str(cts.spectrum_i))

    def update_fitpar_fn(self):
        ''' Updating the fit parameters on the window'''

        if (self.calibBool):
            self.afit_lb.setText("{:.3e}".format(cts.afit))
            self.t0fit_lb.setText("{:.3e}".format(cts.t0fit))
            self.cfit_lb.setText("{:.3f}".format(cts.cfit))

    def update_envect_fn(self):
        ''' Updating th energy vector parameters with the values written in the associated QLineEdit objects'''

        cts.elow = float(self.elow_le.text())
        cts.ehigh = float(self.ehigh_le.text())
        cts.dE = float(self.dE_le.text())
        self.elow_le.setText("{:.2f}".format(cts.elow))
        self.ehigh_le.setText("{:.2f}".format(cts.ehigh))
        self.window().updateglobvar_fn()

    def gas_combo_lr(self, i):
        ''' Gas QCombobox listener'''
        cts.cur_Ip = cts.IPLIST[i]
        cts.first_harm = cts.FIRST_HARMLIST[i]
        self.firstharm_le.setText(str(cts.first_harm))
        cts.SBi = cts.first_harm + 1
        cts.elow = (cts.first_harm - 1) * cts.HEV * cts.cur_nu
        self.elow_le.setText("{:.2f}".format(cts.elow))
        self.update_cts_fn()

    def threshtype_lr(self):
        ''' Listener of the threshold radiobuttons: Y, Xmin or Xmax'''
        rb = self.sender()
        self.threshtype = rb.value

    def withsb_fn(self, state):
        ''' "with sidebands" checkbox listener '''
        if state == Qt.Checked:
            self.withsb_bool = True
        else:
            self.withsb_bool = False

    def showexppeaks_lr(self, state):
        ''' "show expected peaks" checkbox listener '''
        if state == Qt.Checked:
            self.showexppeaksBool = True
        else:
            self.showexppeaksBool = False
        self.refreshplot_fn()

    def setth_lr(self):
        ''' "set threshold" checkbox listener '''
        if self.setth_cb.isChecked():
            self.addpeak_cb.setCheckState(Qt.Unchecked)

    def addpeak_lr(self):
        ''' "add peak" checkbox listener '''
        if self.addpeak_cb.isChecked():
            self.setth_cb.setCheckState(Qt.Unchecked)

    def removepeaks_lr(self):
        ''' "remove peaks" button listener '''
        rmp = rmPeaksDialog(self) # new class defined below

    def refreshplot_fn(self):
        ''' Updating the top left (TOF) graph'''
        xmin, xmax = self.tof_ax.get_xlim()
        ymin, ymax = self.tof_ax.get_ylim()
        self.tof_ax.cla()

        # self.tof_ax.xaxis.set_major_formatter(FormatStrFormatter('%2.e'))
        # self.tof_ax.yaxis.set_major_formatter(FormatStrFormatter('%1.e'))
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
                try:
                    xval = float(np.math.sqrt(
                        0.5 * cts.ME * cts.cur_L ** 2 / cts.QE / (qq2[i] * cts.HEV * cts.cur_nu - cts.cur_Ip)) + 6e-8)
                    # NB: cts.QE is used to convert the energy from eV to Joules. It's not the electron's charge
                    x = np.full((100, 1), xval)
                except Exception:
                    print(traceback.format_exception(*sys.exc_info()))
                self.tof_ax.plot(x, y, color=c, linewidth=1.0)

        if self.bgndremoved:
            self.bgndremoved = False #this means that when we remove bgnd, we don't keep the same scale
        else:
            self.tof_ax.set_ylim(ymin, ymax)

        self.tof_ax.set_xlim(xmin, xmax)
        self.tof_fc.draw()

    def onclick(self, event):
        ''' called when double-clicking on the TOF graph'''
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

    def exportXUV_lr(self):
        ''' "Export XUV" button listener. Saving the energy vector and the energy-converted TOF signal'''
        filename = QFileDialog.getSaveFileName(self, 'Save XUV')
        fname = filename[0]
        if fname:
            XUV_array = np.vstack((self.Eevlin, self.signal)).T
            np.savetxt(fname, XUV_array, delimiter='\t')

    def clear_lr(self):
        ''' "Clear" button listener. Resets all the objects, but not the global variables'''
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
        self.minus_sign_cb.setEnabled(False)

        for w in self.children():
            if isinstance(w, QPushButton) or isinstance(w, QCheckBox):
                w.setEnabled(False)
        self.load_btn.setEnabled(True)

        self.minus_sign_cb.setEnabled(False)

        self.tof_fc.draw()
        self.fit_fc.draw()
        self.en_fc.draw()

        cts.clear_varlist()
        self.window().updateglobvar_fn()

class rmPeaksDialog(QDialog):
    ''' created when clicking on the "remove peaks" button'''
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
    def rmpeak_lr(self):
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