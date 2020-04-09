import sys, traceback
from PyQt5.QtWidgets import QWidget, QPushButton, QApplication, QGridLayout, QHBoxLayout, QVBoxLayout, QGroupBox, QLabel
from PyQt5.QtWidgets import QLineEdit, QSplitter, QFrame, QCheckBox, QButtonGroup, QRadioButton
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
import numpy as np

# Homemade modules
import glob_var as cts
import other_widgets as ow
import _import_export as ie
from Rabbit import FT_contrast_win as ftcw, _Rabbit_functions as Raf, SB_vs_delay_win as sbdw


class RabbitWin(QWidget, ie.Imp_Exp_Mixin, Raf.Rabbit_functions_mixin):
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

        In this file are defined the graphical objects and a few simple functions. The other functions
        (like the analysis ones) are defined in the "_Rabbit_functions.py" file
        """

    ##################################################################################
    ############################ Widget Initialization ###############################

    def __init__(self, parent=None):
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
        self.init_select_bands_layout()
        self.init_rabbitlayout()
        self.init_exportlayout()
        self.init_plotbtnlayout()
        self.init_graphlayout()

        self.mainlayout.addWidget(self.splitter2)
        self.mainlayout.addLayout(self.commandLayout)
        self.setLayout(self.mainlayout)
        self.show()

    def init_var(self):
        ''' Initialization of instance attributes'''
        self.dataloaded = False
        self.xuvonlyloaded = False
        self.bandselected = False

        self.data_tof = []
        self.tof = []
        self.toflength = 0

        self.elength = 0

        self.SBi = 0
        self.ampl = []
        self.ang = []
        self.fpeak = []
        self.peak = []
        self.peak_phase = []
        self.freqnorm = []
        self.pa = []

        self.ampl_rainbow = []
        self.ampl_rainbow2 = []
        self.ang_rainbow = []
        self.fpeak_rainbow = []
        self.peak_rainbow = []
        self.peak_phase_rainbow = []
        self.energy_rainbow = []
        self.energy_rainbow2 = []

        self.fpeak_index = 0

        self.clickcount = 0 # used in the onclick function

        self.rabbit_done = False

    def init_importlayout(self):
        ''' In commandLayout - Initialization of the "Import" section'''
        Importlayout = QGridLayout()
        Importlayout.setSpacing(10)

        Import_box = QGroupBox(self)
        Import_box.setTitle("Import")
        Import_box.setFixedSize(300, 120)

        self.importcalib_btn = QPushButton("calib", self)
        self.importdata_btn = QPushButton("data", self)
        self.importXUV_btn = QPushButton("XUV only", self)
        self.importRabbitVMI_btn = QPushButton("VMI RABBIT", self)
        self.vrep_le = QLineEdit(str(cts.vrep), self)

        self.need_calib_cb = QCheckBox("Calibrate VMI ?", self)

        self.importcalib_btn.clicked.connect(self.importcalib_lr)
        self.importdata_btn.clicked.connect(self.importdata_lr)
        self.importXUV_btn.clicked.connect(self.importXUV_lr)
        self.importRabbitVMI_btn.clicked.connect(self.importRabbitVMI_lr)
        self.vrep_le.returnPressed.connect(self.update_vrep)

        Importlayout.addWidget(self.importcalib_btn, 0, 0)
        Importlayout.addWidget(self.importdata_btn, 0, 1)
        Importlayout.addWidget(self.importXUV_btn, 0, 2)
        Importlayout.addWidget(self.importRabbitVMI_btn, 1, 0)
        Importlayout.addWidget(QLabel("Vrep:"), 1, 1)
        Importlayout.addWidget(self.vrep_le, 1, 2)
        Importlayout.addWidget(self.need_calib_cb, 2, 2)


        Import_box.setLayout(Importlayout)
        self.commandLayout.addWidget(Import_box)

        for widget in Import_box.children():
            if isinstance(widget, QPushButton):
                widget.setSizePolicy(0, 0)
                widget.setEnabled(False)

        self.vrep_le.setSizePolicy(0, 0)
        self.vrep_le.setFixedSize(40, 20)

        self.importcalib_btn.setEnabled(True)
        self.importRabbitVMI_btn.setEnabled(True)

    def init_envectlayout(self):
        ''' In commandLayout - Initialization of the energy vector section, with elow, ehigh and dE'''
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

    def init_sigtreatmentlayout(self):
        ''' In commandLayout - Initialization of the "Signal Treatment" section'''
        sigtreatment_box = QGroupBox("Signal Treatment", self)
        sigtreatmentlayout = QGridLayout()

        sigtreatment_box.setFixedSize(300, 100)

        self.smooth_btn = QPushButton("smooth", self)
        self.smooth_le = QLineEdit("2", self)
        self.normalize_btn = QPushButton("normalize", self)
        self.subXUV_btn = QPushButton("subtract XUV", self)

        self.smooth_btn.clicked.connect(self.smoothrab_lr)
        self.normalize_btn.clicked.connect(self.normalizerab_lr)
        self.subXUV_btn.clicked.connect(self.subXUV_lr)

        sigtreatmentlayout.addWidget(self.smooth_btn, 0, 0)
        sigtreatmentlayout.addWidget(self.smooth_le, 1, 0)
        sigtreatmentlayout.addWidget(self.normalize_btn, 0, 1)
        sigtreatmentlayout.addWidget(self.subXUV_btn, 0, 2)

        sigtreatment_box.setLayout(sigtreatmentlayout)

        for widget in sigtreatment_box.children():
            if isinstance(widget, QPushButton):
                widget.setSizePolicy(0, 0)
                widget.setEnabled(False)
            if isinstance(widget, QLineEdit):
                widget.setFixedSize(55, 20)

        self.commandLayout.addWidget(sigtreatment_box)

    def init_select_bands_layout(self):
        layout = QGridLayout()
        select_bands_box = QGroupBox("Select Bands", self)
        select_bands_box.setFixedSize(300, 150)

        self.color_bg = QButtonGroup(self)
        self.white_rb = QRadioButton("white", self)
        self.white_rb.value = "w"
        self.white_rb.toggle()
        self.white_rb.toggled.connect(self.color_lr)

        self.black_rb = QRadioButton("black", self)
        self.black_rb.value = "k"
        self.black_rb.toggled.connect(self.color_lr)

        self.color_bg.addButton(self.black_rb)
        self.color_bg.addButton(self.white_rb)

        self.selectbands_cb = QCheckBox("Select", self)
        self.selectbands_cb.stateChanged.connect(self.selectbands_lr) # defined in _Rabbit_functions

        self.showbands_cb = QCheckBox("Show", self)
        self.showbands_cb.stateChanged.connect(self.showbands_lr) # defined in _Rabbit_functions

        self.selectmode_bg = QButtonGroup(self)
        self.manualselect_rb = QRadioButton("manual", self)
        self.manualselect_rb.value = "manual"
        self.selectmode = "manual"
        self.manualselect_rb.toggle()
        self.manualselect_rb.toggled.connect(self.selectmode_lr)

        self.FTselect_rb = QRadioButton("with FT", self)
        self.FTselect_rb.value = "FT"
        self.FTselect_rb.toggled.connect(self.selectmode_lr)

        self.theoryselect_rb = QRadioButton("theory", self)
        self.theoryselect_rb.value = "theory"
        self.theoryselect_rb.toggled.connect(self.selectmode_lr)

        self.selectmode_bg.addButton(self.manualselect_rb)
        self.selectmode_bg.addButton(self.FTselect_rb)
        self.selectmode_bg.addButton(self.theoryselect_rb)

        self.FTthreshold_le = QLineEdit("0.7", self)
        self.FTthreshold_le.returnPressed.connect(self.FTselect_fn) # defined in _Rabbit_functions

        self.theorythreshold_le = QLineEdit("0.3", self)
        self.theorythreshold_le.returnPressed.connect(self.theoryselect_fn)

        self.integrate_cb = QCheckBox("integrate", self)
        self.integrate_cb.toggle()
        self.integrate_cb.stateChanged.connect(self.integrate_lr)

        self.bandsnb_le = QLineEdit(str(cts.bandsnb), self)
        self.bandsnb_le.returnPressed.connect(self.bandsnb_lr)

        self.SBi_le = QLineEdit(str(cts.first_harm + 1), self)
        self.SBi_le.returnPressed.connect(self.SBi_lr)

        self.done_btn = QPushButton("Done", self)
        self.done_btn.clicked.connect(self.done_lr) # defined in _Rabbit_functions

        self.clearbands_btn = QPushButton("Clear", self)
        self.clearbands_btn.clicked.connect(self.clearbands_lr) # defined in _Rabbit_functions

        layout.addWidget(self.black_rb, 3, 2)
        layout.addWidget(self.white_rb, 3, 3)
        layout.addWidget(self.selectbands_cb, 1, 1)
        layout.addWidget(self.showbands_cb, 1, 0)
        layout.addWidget(self.manualselect_rb, 0, 1)
        layout.addWidget(self.FTselect_rb, 0, 2)
        layout.addWidget(self.FTthreshold_le, 1, 2)
        layout.addWidget(self.theoryselect_rb, 0, 3)
        layout.addWidget(self.theorythreshold_le, 1, 3)
        layout.addWidget(self.integrate_cb, 0, 0)
        layout.addWidget(QLabel("number:"), 2, 0)
        layout.addWidget(self.bandsnb_le, 2, 1)
        layout.addWidget(QLabel("first SB:"), 2, 2)
        layout.addWidget(self.SBi_le, 2, 3)
        layout.addWidget(self.done_btn, 3, 0)
        layout.addWidget(self.clearbands_btn, 3, 1)

        select_bands_box.setLayout(layout)

        for widget in select_bands_box.children():
            if isinstance(widget, QPushButton):
                widget.setSizePolicy(0, 0)
                widget.setEnabled(False)
            if isinstance(widget, QLineEdit):
                widget.setFixedSize(55, 20)
                widget.setSizePolicy(0, 0)

        self.commandLayout.addWidget(select_bands_box)

    def init_rabbitlayout(self):
        ''' In commandLayout - Initialization of the "RABBIT" section'''
        rabbitlayout = QGridLayout()
        rabbit_box = QGroupBox("RABBIT", self)
        rabbit_box.setFixedSize(300, 100)

        self.normalrab_btn = QPushButton("Normal", self)
        self.FTcontrast_btn = QPushButton("FT/Contrast", self)
        self.rainbowrab_btn = QPushButton("Rainbow", self)
        self.clear_btn = QPushButton("Clear", self)
        self.lock_2w_cb = QCheckBox("Lock 2w", self)
        self.set_2w_le = QLineEdit("0", self)

        self.normalrab_btn.clicked.connect(self.normal_rabbit_lr)
        self.FTcontrast_btn.clicked.connect(self.FTcontrast_lr)
        self.rainbowrab_btn.clicked.connect(self.rainbow_rabbit_lr)
        self.clear_btn.clicked.connect(self.clear_lr)
        self.lock_2w_cb.stateChanged.connect(self.lock_2w_lr)  # defined in _Rabbit_functions
        self.set_2w_le.returnPressed.connect(self.set_2w_lr)  # defined in _Rabbit_functions

        rabbitlayout.addWidget(self.normalrab_btn, 0, 0)
        rabbitlayout.addWidget(self.FTcontrast_btn, 0, 1)
        rabbitlayout.addWidget(self.rainbowrab_btn, 0, 2)
        rabbitlayout.addWidget(self.clear_btn, 1, 0)
        rabbitlayout.addWidget(self.lock_2w_cb, 1, 1)
        rabbitlayout.addWidget(self.set_2w_le, 1, 2)

        self.set_2w_le.setFixedSize(55, 20)
        self.set_2w_le.setEnabled(False)
        self.lock_2w_cb.setEnabled(False)

        rabbit_box.setLayout(rabbitlayout)

        for widget in rabbit_box.children():
            if isinstance(widget, QPushButton):
                widget.setSizePolicy(0, 0)
                widget.setEnabled(False)

        self.commandLayout.addWidget(rabbit_box)

    def init_exportlayout(self):
        ''' In commandLayout - Initialization of the "Export" section'''
        exportlayout = QGridLayout()
        export_box = QGroupBox("Export", self)
        export_box.setFixedSize(300, 60)

        self.exportrab_btn = QPushButton("RABBIT", self)
        # self.export_2w_btn = QPushButton("2w", self)
        self.export_2w_btn = QPushButton("2w Rainbow", self)
        self.export_2w_std_btn = QPushButton("2w std RABBIT", self)

        self.exportrab_btn.clicked.connect(self.exportrab_lr)
        self.export_2w_btn.clicked.connect(self.export_2w_lr)
        self.export_2w_std_btn.clicked.connect(self.export_2w_std_lr)

        exportlayout.addWidget(self.exportrab_btn, 0, 0)
        exportlayout.addWidget(self.export_2w_btn, 0, 1)
        exportlayout.addWidget(self.export_2w_std_btn, 0, 2)

        export_box.setLayout(exportlayout)

        for widget in export_box.children():
            if isinstance(widget, QPushButton):
                widget.setSizePolicy(0, 0)
                widget.setEnabled(False)
        self.commandLayout.addWidget(export_box)

    def init_plotbtnlayout(self):
        ''' In commandLayout - Initialization of the "Plot" section'''
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

    def init_graphlayout(self):
        ''' In graphLayout - Initialization of the 3 figures'''
        self.splitter1 = QSplitter(Qt.Horizontal)
        self.splitter2 = QSplitter(Qt.Vertical)

        self.rab_widget = ow.plot3DWidget(self) # new class defined in other_widgets.py
        self.rab_widget.fc.mpl_connect('button_press_event', self.onclick)
        self.rab_widget.xlabel = "E [eV]"
        self.rab_widget.ylabel = "t [fs]"

        self.phaselayout = QVBoxLayout()
        self.phase_frame = QFrame()
        self.phase_frame.setLayout(self.phaselayout)

        phase_fig = Figure(figsize=(2, 2), dpi=100)
        self.phase_fc = FigureCanvas(phase_fig)
        self.phase_fc.setSizePolicy(1, 1)
        self.phase_fc.mpl_connect('button_press_event', self.onclick)
        self.phase_ax = self.phase_fc.figure.add_subplot(111)
        self.phase_ax.tick_params(labelsize = 10)
        nav = NavigationToolbar2QT(self.phase_fc, self)
        nav.setStyleSheet("QToolBar { border: 0px }")
        self.phaselayout.addWidget(self.phase_fc)
        self.phaselayout.addWidget(nav)
        self.phase_fc.draw()

        self.FTlayout = QVBoxLayout()
        self.FT_frame = QFrame()
        self.FT_frame.setLayout(self.FTlayout)
        FT_fig = Figure(figsize=(2, 2), dpi=100)
        self.FT_fc = FigureCanvas(FT_fig)
        self.FT_fc.setSizePolicy(1, 1)
        self.FT_fc.mpl_connect('button_press_event', self.onclick)
        self.FT_ax = self.FT_fc.figure.add_subplot(111)
        self.FT_ax.tick_params(labelsize = 10)
        nav2 = NavigationToolbar2QT(self.FT_fc, self)
        nav2.setStyleSheet("QToolBar { border: 0px }")
        self.FTlayout.addWidget(self.FT_fc)
        self.FTlayout.addWidget(nav2)
        self.int_rab_cb = QCheckBox("Delay integrated trace")
        self.int_rab_cb.setEnabled(False)
        self.int_rab_cb.stateChanged.connect(self.int_rab_lr)
        self.FTlayout.addWidget(self.int_rab_cb)
        self.FT_fc.draw()

        self.splitter1.addWidget(self.phase_frame)
        self.splitter1.addWidget(self.FT_frame)

        self.splitter2.addWidget(self.rab_widget)
        self.splitter2.addWidget(self.splitter1)

#################################################################################
############################ Other methods ######################################

    def update_vrep(self):
        try:
            cts.vrep = float(self.vrep_le.text())
        except ValueError:
            pass
        finally:
            self.vrep_le.setText(str(cts.vrep))

    def int_rab_lr(self):
        self.refreshplot(keep_lims=True)

    def bandsnb_lr(self):
        ''' called when pressing enter in the bandsnb_le object'''
        try:
            cts.bandsnb = int(self.bandsnb_le.text())
            self.window().updateglobvar_fn()
            cts.bands_vect = np.zeros([cts.bandsnb, 2])
        except ValueError:
            self.window().statusBar().showMessage("Number of bands must be an integer")

    def SBi_lr(self):
        try:
            cts.SBi = float(self.SBi_le.text())
            self.selectmode_fn()
        except ValueError:
            self.window().statusBar().showMessage("SBi must be a number")

    def subXUV_lr(self):
        ''' "substract XUV" button listener. Opens a new window'''
        cts.xuvsubstracted = False
        sw = subXUVWin(self) # new class defined below

    def FTcontrast_lr(self):
        ''' "FT/contrast" button listener. Opens a new window'''
        try:
            w = ftcw.FTContrastWin(self) # new class defined below
        except Exception:
            print(traceback.format_exception(*sys.exc_info()))

    def plotSBvsdelay_lr(self):
        ''' "[plot] SB vs delay" button listener. Opens a new window'''
        try:
            w = sbdw.SBvsDelayWin(self) # new class defined below
        except Exception:
            print(traceback.format_exception(*sys.exc_info()))

    def update_scanparam(self):
        ''' Updates the values of the scan steps, in nm and fs'''
        try:
            cts.scanstep_nm = float(self.scanparam_le.text())
            cts.scanstep_fs = float(self.scanparam_le.text()) * 2 / (cts.C * 1e-6)
            self.window().updateglobvar_fn()
        except ValueError:
            print('Scan step must be a number')

    def update_envect_fn(self):
        ''' Updates the energy vector parameters'''

        cts.elow = float(self.elow_le.text())
        cts.ehigh = float(self.ehigh_le.text())
        cts.dE = float(self.dE_le.text())
        self.elow_le.setText("{:.2f}".format(cts.elow))
        self.ehigh_le.setText("{:.2f}".format(cts.ehigh))
        self.window().updateglobvar_fn()

    def onclick(self, event):
        """ Updates cts.bands_vect with the energy value(s) obtained by clicking
         on one the three graphs (only if self.selectbands_cb checkbox is checked and
         the selection mode is 'manual'). The graphs are then replotted with the
         corresponding new energy line(s).

        :param event:
        :return: None
        """
        if self.selectbands_cb.isChecked() and self.selectmode == "manual":
            if self.clickcount == 0:
                self.clearbands_btn.setEnabled(True)
            if self.integrate_cb.isChecked():
                if self.clickcount < 2 * cts.bandsnb:
                    if self.clickcount == 2 * cts.bandsnb - 1:
                        self.done_btn.setEnabled(True)
                    cts.bands_vect[self.clickcount // 2, self.clickcount % 2] = event.xdata
                    self.clickcount += 1
                # when all the lines are set, we can modify a line on clicking next to it
                else:
                    k = np.argmin(abs(cts.bands_vect - event.xdata))
                    cts.bands_vect[k // 2, k % 2] = event.xdata
                self.refreshplot(keep_lims=True)
            else:
                if self.clickcount < cts.bandsnb:
                    if self.clickcount == cts.bandsnb - 1:
                        self.done_btn.setEnabled(True)
                    cts.bands_vect[self.clickcount, 0] = event.xdata
                    cts.bands_vect[self.clickcount, 1] = event.xdata
                    self.clickcount += 1
                # when all the lines are set, we can modify a line on clicking next to it
                else:
                    k = np.argmin(abs(cts.bands_vect[:,0] - event.xdata))
                    cts.bands_vect[k, 0] = event.xdata
                    cts.bands_vect[k, 1] = event.xdata
                self.refreshplot(keep_lims=True)

    def integrate_lr(self):
        self.selectmode_fn()

    def selectmode_lr(self):
        rb = self.sender()
        self.selectmode = rb.value
        self.selectmode_fn()

    def selectmode_fn(self):
        if self.selectmode == "manual":
            self.FTthreshold_le.setEnabled(False)
            self.theorythreshold_le.setEnabled(False)
            if self.clickcount == 2*cts.bandsnb:
                self.done_btn.setEnabled(True)
            else:
                self.done_btn.setEnabled(False)
            self.refreshplot(keep_lims=True)

        elif self.selectmode == "FT":
            if self.integrate_cb.isChecked():
                self.FTthreshold_le.setEnabled(True)
            else:
                self.FTthreshold_le.setEnabled(False)
            self.theorythreshold_le.setEnabled(False)
            self.FTselect_fn() # defined in _Rabbit_functions
            self.done_btn.setEnabled(True)

        elif self.selectmode == "theory":
            self.FTthreshold_le.setEnabled(False)
            if self.integrate_cb.isChecked():
                self.theorythreshold_le.setEnabled(True)
            else:
                self.theorythreshold_le.setEnabled(False)
            self.theoryselect_fn()
            self.done_btn.setEnabled(True)

    def color_lr(self) -> None:
        rb = self.sender()
        self.color = rb.value
        self.refreshplot(keep_lims=True)

    def reset_btn(self):
        ''' "Reset" button listener. Resets the widgets, not the variables'''
        self.importXUV_btn.setEnabled(False)
        self.importdata_btn.setEnabled(False)
        self.smooth_btn.setEnabled(False)
        self.normalize_btn.setEnabled(False)
        self.subXUV_btn.setEnabled(False)
        self.normalrab_btn.setEnabled(False)
        self.FTcontrast_btn.setEnabled(False)
        self.rainbowrab_btn.setEnabled(False)
        self.exportrab_btn.setEnabled(False)
        self.plotSBvsdelay_btn.setEnabled(False)

        self.rab_widget.colorauto_cb.setEnabled(False)
        self.rab_widget.logscale_cb.setEnabled(False)




if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = RabbitWin()
    sys.exit(app.exec_())
