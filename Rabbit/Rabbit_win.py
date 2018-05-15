import sys, traceback
from PyQt5.QtWidgets import QWidget, QPushButton, QApplication, QGridLayout, QHBoxLayout, QVBoxLayout, QGroupBox, QLabel
from PyQt5.QtWidgets import QLineEdit
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure

import numpy as np

# Homemade modules
import glob_var as cts
import other_widgets as ow
import _import_export as ie
from Rabbit import FT_contrast_win as ftcw, Rabbit_select_bands as Rsb, _Rabbit_functions as Raf, SB_vs_delay_win as sbdw


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
        """
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
        self.init_rabbitlayout()
        self.init_exportlayout()
        self.init_plotbtnlayout()
        self.init_graphlayout()

        self.mainlayout.addLayout(self.graphlayout)
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

    def init_importlayout(self):
        ''' In commandLayout - Initialization of the "Import" section'''
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

    def init_rabbitlayout(self):
        ''' In commandLayout - Initialization of the "RABBIT" section'''
        rabbitlayout = QGridLayout()
        rabbit_box = QGroupBox("RABBIT", self)
        rabbit_box.setFixedSize(300, 100)

        self.normalrab_btn = QPushButton("Normal", self)
        self.FTcontrast_btn = QPushButton("FT/Contrast", self)
        self.rainbowrab_btn = QPushButton("Rainbow", self)
        self.clear_btn = QPushButton("Clear", self)

        self.normalrab_btn.clicked.connect(self.normalrab_lr)
        self.FTcontrast_btn.clicked.connect(self.FTcontrast_lr)
        self.rainbowrab_btn.clicked.connect(self.rainbowrab_lr)
        self.clear_btn.clicked.connect(self.clear_lr)

        rabbitlayout.addWidget(self.normalrab_btn, 0, 0)
        rabbitlayout.addWidget(self.FTcontrast_btn, 0, 1)
        rabbitlayout.addWidget(self.rainbowrab_btn, 0, 2)
        rabbitlayout.addWidget(self.clear_btn, 1, 0)

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

        self.exportrab_btn.clicked.connect(self.exportrab_lr)

        exportlayout.addWidget(self.exportrab_btn)

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

    def selectbands_lr(self):
        ''' "select bands" button listener'''
        try:
            cts.bandsnb = int(self.bandsnb_le.text())
            cts.bands_vect = np.zeros([cts.bandsnb, 2])
            self.subXUV_btn.setEnabled(False)
            self.normalrab_btn.setEnabled(False)
            self.rainbowrab_btn.setEnabled(False)
            self.FTcontrast_btn.setEnabled(False)
            self.plotSBvsdelay_btn.setEnabled(False)
            self.window().updateglobvar_fn()
            nw = Rsb.selectBandsWin(self) # new class defined in Rabbit_select_bands.py
        except ValueError:
            self.window().statusBar().showMessage("Number of bands must be an integer")

    def bandsnb_lr(self):
        ''' called when pressing enter in the bandsnb_le object'''
        try:
            cts.bandsnb = int(self.bandsnb_le.text())
            self.window().updateglobvar_fn()
        except ValueError:
            self.window().statusBar().showMessage("Number of bands must be an integer")

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

    def reset_btn(self):
        ''' "Reset" button listener. Resets the widgets, not the variables'''
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

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = RabbitWin()
    sys.exit(app.exec_())
