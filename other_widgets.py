
import sys, os, traceback
from PyQt5.QtCore import QRect, Qt
from PyQt5.QtWidgets import QWidget, QPushButton, QApplication, QGridLayout, QHBoxLayout, QVBoxLayout, QGroupBox, QLabel
from PyQt5.QtWidgets import QLineEdit, QCheckBox, QComboBox, QRadioButton, QButtonGroup, QFileDialog, QDialog, QTableWidget
from PyQt5.QtWidgets import QTableWidgetItem, QTableView, QLayout, QSlider
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
from glob import glob
import numpy as np

#homemade modules
import glob_var as cts

class plot3DWidget(QWidget):
    def __init__(self, parent):
        super(plot3DWidget, self).__init__(parent)
        #self.setGeometry(300, 300, 500, 300)
        self.setMinimumHeight(400)

        self.mainlayout = QVBoxLayout()
        self.mainlayout.setContentsMargins(0, 0, 0, 0)

        self.init_var()
        self.init_graphlayout()
        self.init_optionslayout()

        self.mainlayout.addLayout(self.graphlayout)
        self.mainlayout.addLayout(self.optionslayout)

        self.setLayout(self.mainlayout)
        self.show()

    def init_var(self):
        self.xlabel = ""
        self.ylabel = ""
        self.logscale = False
        self.x = np.zeros([1,1])
        self.y = np.zeros([1,1])
        self.signal = np.zeros([1,1])
        self.dataplotted = False

    def init_graphlayout(self):
        self.graphlayout = QVBoxLayout()

        fig = Figure(figsize=(4, 3), dpi=100)
        self.fc = FigureCanvas(fig)
        self.ax = self.fc.figure.add_subplot(111)
        self.fc.draw()
        nav = NavigationToolbar2QT(self.fc, self)
        nav.setStyleSheet("QToolBar { border: 0px }")

        self.graphlayout.addWidget(self.fc)
        self.graphlayout.addWidget(nav)

    def init_optionslayout(self):
        self.optionslayout = QVBoxLayout()

        self.colorlayout = QHBoxLayout()

        self.color_sld = QSlider(Qt.Horizontal, self)
        self.color_sld.setFocusPolicy(Qt.NoFocus)
        self.color_sld.setGeometry(0, 0, 100, 30)
        self.color_sld.valueChanged[int].connect(self.colorsld_lr)
        self.color_sld.setEnabled(False)

        self.colorauto_cb = QCheckBox("Color auto", self)
        self.colorauto_cb.stateChanged.connect(self.colorauto_lr)
        self.colorauto_cb.toggle()
        self.colorauto_cb.setEnabled(False)

        self.logscale_cb = QCheckBox("log scale", self)
        self.logscale_cb.stateChanged.connect(self.logscale_lr)
        self.logscale_cb.setEnabled(False)

        self.colorlayout.addWidget(self.color_sld)
        self.colorlayout.addWidget(self.colorauto_cb)
        self.colorlayout.addWidget(self.logscale_cb)

        self.optionslayout.addLayout(self.colorlayout)

    def colorauto_lr(self):
        if self.colorauto_cb.isChecked():
            self.color_sld.setValue(25)
            self.color_sld.setEnabled(False)
        else:
            self.color_sld.setEnabled(True)

    def logscale_lr(self):
        if self.logscale_cb.isChecked():
            self.logscale = True
            self.colorauto_cb.setCheckState(Qt.Checked)
            self.colorauto_cb.setEnabled(False)
            self.color_sld.setEnabled(False)
        else:
            self.logscale = False
            self.colorauto_cb.setEnabled(True)

        self.refreshplot_fn()

    def colorsld_lr(self, value):
        if(self.dataplotted):
            self.refreshplot_fn(value=value)

    def refreshplot_fn(self, x=None, y=None, signal=None, value=25):
        if x is None:
            x = self.x
        else:
            self.x = x
        if y is None:
            y = self.y
        else:
            self.y = y
        if signal is None:
            signal = self.signal
        else:
            self.signal = signal

        self.ax.cla()

        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)

        val = float(value)

        if self.logscale:
            sig = np.log(signal)
        else:
            sig = signal

        vmin = sig.min()
        #vmax = (sig.max()) * (1 + (25 - val) / 75) + sig.min()
        vmax = 1/75*((sig.min() - sig.max())*val + 100*sig.max() - 25*sig.min())

        if(self.colorauto_cb.isChecked()):
            im = self.ax.imshow(sig, extent=(x[0],x[-1], y[0], y[-1]), aspect='auto', origin='lower',
                                                    interpolation='nearest')
        else:
            im = self.ax.imshow(sig, extent=(x[0], x[-1], y[0], y[-1]), aspect='auto', origin='lower',
                           interpolation='nearest', vmin=vmin, vmax=vmax)
        self.fc.draw()
        if self.dataplotted == False:
            self.colorauto_cb.setEnabled(True)
            self.logscale_cb.setEnabled(True)
        self.dataplotted = True
