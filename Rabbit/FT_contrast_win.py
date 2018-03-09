# FT_contrast_win.py

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QDialog, QHBoxLayout, QVBoxLayout, QCheckBox, QGridLayout, QGroupBox, QPushButton, QLineEdit
from PyQt5.QtWidgets import QLabel, QTableWidget, QTableView, QTableWidgetItem
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure

import numpy as np
import traceback, sys

import glob_var as cts
import analysis_functions as af


class FTContrastWin(QDialog):
    ''' created when clicking on the "FT/Contrast" button in the RabbitWin object'''
    def __init__(self, parent): # parent=RabbitWin
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
            jf = np.argmin(abs(cts.energy_vect - xf)) + 1
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
           for ji, jf in self.jx[:]:
               if jf - ji > 1:  # if we integrate ie if we have more than one point
                   x_data = np.trapz(cts.rabbit_mat[:, ji:jf], cts.energy_vect[ji:jf], axis=1)
               else:
                   x_data = cts.rabbit_mat[:, ji]

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
