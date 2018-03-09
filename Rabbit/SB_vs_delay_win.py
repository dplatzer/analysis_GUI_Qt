# SB_vs_delay_win

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QDialog, QHBoxLayout, QVBoxLayout, QCheckBox, QGridLayout
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT
from PyQt5.QtWidgets import QTableWidgetItem, QTableView, QTableWidget
from matplotlib.figure import Figure

import numpy as np
from scipy import optimize as opt

import glob_var as cts

class SBvsDelayWin(QDialog):
    ''' created when clicking on the "[plot] SB vs delay" button'''
    def __init__(self, parent): # parent=RabbitWin
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

    def init_var(self):

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
        self.a0 = []
        self.a1 = []
        self.phi = []

        for i in range(cts.bandsnb):
            popt, pcov = opt.curve_fit(self.cosfit_fn, cts.delay_vect, self.SB[i,:], np.array([1e-11, 1e-12, 0]))

            if popt[1] < 0:
                popt[1] = -1 * popt[1]
                popt[2] = popt[2] + np.pi

            cts.contrast[i, 1] = popt[1] / popt[0]
            self.a0.append(popt[0])
            self.a1.append(popt[1])
            self.phi.append(popt[2])

        print(self.a1)
        self.parent().window().updateglobvar_fn()

    def init_graphlayout(self) -> None:
        self.graphlayout = QVBoxLayout()

        fig = Figure(figsize=(4, 3), dpi=100)
        self.fc = FigureCanvas(fig)
        self.ax = self.fc.figure.add_subplot(111)
        self.fc.draw()
        nav = NavigationToolbar2QT(self.fc, self)
        nav.setStyleSheet("QToolBar { border: 0px }")

        self.checkbox_layout = QGridLayout()
        self.checkbox_table = []
        for i in range(cts.bandsnb):
            cb = QCheckBox(str(cts.first_harm + 2*i +1), self)
            cb2 = QCheckBox((str(cts.first_harm + 2*i +1) + "_cosfit"), self)

            cb.setCheckState(Qt.Checked)
            cb2.setCheckState(Qt.Unchecked)

            cb.stateChanged.connect(self.refreshplot_fn)
            cb2.stateChanged.connect(self.refreshplot_fn)

            self.checkbox_layout.addWidget(cb, 0, i)
            self.checkbox_layout.addWidget(cb2, 1, i)

            self.checkbox_table.append(cb)
            self.checkbox_table.append(cb2)

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

            for i in range(2*cts.bandsnb):
                if self.checkbox_table[i].isChecked() and i%2 == 0:
                    self.ax.plot(cts.delay_vect, self.SB[i//2,:], label="SB%d" % (cts.first_harm + 2 * i//2 + 1))
                if self.checkbox_table[i].isChecked() and i%2 != 0:
                    self.ax.plot(cts.delay_vect, self.cosfit_fn(cts.delay_vect, self.a0[i//2], self.a1[i//2], self.phi[i//2]),
                                 label="SB%d_cosfit" % (cts.first_harm + 2 * i//2 + 1))

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
