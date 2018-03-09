# sub_XUV_win.py

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QDialog, QHBoxLayout, QVBoxLayout, QPushButton
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
import numpy as np
import traceback, sys

import glob_var as cts
import analysis_functions as af
import other_widgets as ow

class subXUVWin(QDialog):
    ''' created when clicking on the "Substract XUV" button" in RabbitWin'''
    def __init__(self, parent): # parent=RabbitWin
        super(subXUVWin, self).__init__(parent)
        self.layout = QHBoxLayout()
        #self.layout.setContentsMargins(0, 0, 0, 0)
        self.setWindowFlags(Qt.Window)

        cts.xuvsubstracted = False
        self.init_layout()
        self.rabgraph.refreshplot_fn(cts.energy_vect, cts.delay_vect, cts.rabbit_mat)

        self.setLayout(self.layout)
        self.show()

    def init_layout(self):
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