# _Rabbit_select_bands.py

from PyQt5.QtWidgets import QHBoxLayout, QDialog, QButtonGroup, QRadioButton, QCheckBox, QLineEdit, QPushButton
from PyQt5.QtCore import Qt
import numpy as np
import traceback

import other_widgets as ow
import glob_var as cts
import analysis_functions as af


class selectBandsWin(QDialog):
    ''' created when clicking on the "Select bands" button'''
    def __init__(self, parent):
        super(selectBandsWin, self).__init__(parent)
        self.layout = QHBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setGeometry(100, 100, 1000, 800)
        self.setWindowFlags(Qt.Window)

        win = selectBandsPlotWin(self) # new class defined below
        self.layout.addWidget(win)

        self.setLayout(self.layout)
        self.show()

class selectBandsPlotWin(ow.plot3DWidget):
    ''' used in the "selectBandsWin" object above'''
    def __init__(self, parent):
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

        self.theoryselect_rb = QRadioButton("theory", self)
        self.theoryselect_rb.value = "theory"
        self.theoryselect_rb.toggled.connect(self.selectmode_lr)

        self.selectmode_bg.addButton(self.manualselect_rb)
        self.selectmode_bg.addButton(self.FTselect_rb)
        self.selectmode_bg.addButton(self.theoryselect_rb)

        self.FTthreshold_le = QLineEdit("0.2", self)
        self.FTthreshold_le.setSizePolicy(0, 0)
        self.FTthreshold_le.setFixedSize(60, 20)
        self.FTthreshold_le.returnPressed.connect(self.FTselect_fn)

        self.theorythreshold_le = QLineEdit("0.3", self)
        self.theorythreshold_le.setSizePolicy(0, 0)
        self.theorythreshold_le.setFixedSize(60, 20)
        self.theorythreshold_le.returnPressed.connect(self.theoryselect_fn)

        self.integrate_cb = QCheckBox("integrate", self)
        self.integrate_cb.toggle()
        self.integrate_cb.stateChanged.connect(self.integrate_lr)

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
        self.paramlayout.addWidget(self.theoryselect_rb)
        self.paramlayout.addWidget(self.theorythreshold_le)
        self.paramlayout.addWidget(self.integrate_cb)
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

    def integrate_lr(self) -> None:
        l=0
        self.selectmode_fn()

    def sbcolor_lr(self, value) -> None:
        self.refreshplot_sb(value)
        self.value = value

    def selectmode_lr(self) -> None:
        rb = self.sender()
        self.selectmode = rb.value
        self.selectmode_fn()

    def selectmode_fn(self) -> None:
        if self.selectmode == "manual":
            self.clickcount = 0
            cts.bands_vect = np.zeros([cts.bandsnb, 2])
            self.FTthreshold_le.setEnabled(False)
            self.theorythreshold_le.setEnabled(False)
            self.done_btn.setEnabled(False)
            self.refreshplot_sb(value=self.value)
        elif self.selectmode == "FT":
            if self.integrate_cb.isChecked():
                self.FTthreshold_le.setEnabled(True)
            else:
                self.FTthreshold_le.setEnabled(False)
            self.theorythreshold_le.setEnabled(False)
            self.FTselect_fn()
            self.done_btn.setEnabled(True)
        elif self.selectmode == "theory":
            self.FTthreshold_le.setEnabled(False)
            if self.integrate_cb.isChecked():
                self.theorythreshold_le.setEnabled(True)
            else:
                self.theorythreshold_le.setEnabled(False)
            self.theoryselect_fn()
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

        if self.integrate_cb.isChecked():
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
                cts.bands_vect[i, 1] = cts.energy_vect[xi[-1] + iemin]
        else:
            for i in range(cts.bandsnb):
                iemin = np.argmin(abs(cts.energy_vect - (self.SBi - 0.3 + 2 * i) * hnu))
                iemax = np.argmin(abs(cts.energy_vect - (self.SBi + 0.3 + 2 * i) * hnu))

                ampl = self.peak_rainbow[iemin:iemax]
                xi = np.argmax(ampl)
                cts.bands_vect[i, 0] = cts.energy_vect[xi + iemin]
                cts.bands_vect[i, 1] = cts.energy_vect[xi + iemin]

        self.refreshplot_sb(value=self.value)

    def theoryselect_fn(self) -> None:
        hnu = cts.HEV * cts.cur_nu
        self.SBi = cts.first_harm + 1


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

        self.refreshplot_sb(value=self.value)

    def onclick(self, event) -> None:
        if self.selectbands_cb.isChecked() and self.selectmode == "manual":
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
                self.refreshplot_sb(value=self.value)
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