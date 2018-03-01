
'''Develop'''
"""The program has to be runned here"""

import sys
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QApplication, QTabWidget, QHBoxLayout, QMainWindow, QSplitter, QTableWidget
from PyQt5.QtWidgets import QTableWidgetItem, QTableView, QDialog
import PyQt5.QtWidgets as qtw
import traceback
import numpy as np

# homemade modules
import calib_win, Rabbit_win
import glob_var as cts

'''The main object.'''
class mainWin(QMainWindow):
    def __init__(self) -> None:
        super(mainWin, self).__init__()
        self.setWindowTitle("Analysis Main")
        self.setGeometry(100, 100, 1100, 750)
        self.statusBar().showMessage("Statusbar - awaiting user control")
        self.show()

        self.container = QWidget()
        self.setCentralWidget(self.container)

        self.mainLayout = qtw.QHBoxLayout(self.container)
        self.mainLayout.setSpacing(0)
        self.mainLayout.setContentsMargins(5, 5, 5, 5)
        self.splitter = QSplitter(Qt.Horizontal)
        self.verticalwidgets()

        self.splitter.addWidget(self.var_table)
        self.splitter.addWidget(self.tabs_widget)

        self.mainLayout.addWidget(self.splitter)

    def verticalwidgets(self) -> None:
        # left side with the variable explorer
        self.var_table = QTableWidget()
        self.var_table.setSelectionBehavior(QTableView.SelectRows)
        self.var_table.doubleClicked.connect(self.var_table_lr)

        self.box_layout1 = qtw.QVBoxLayout()
        self.box_layout1.setContentsMargins(0, 0, 0, 0)
        self.box_layout1.setSpacing(0)

        Vheader = self.var_table.verticalHeader()
        Vheader.setDefaultSectionSize(20)
        Vheader.setVisible(False)

        self.var_table.setUpdatesEnabled(True)

        self.var_table.setColumnWidth(0, 50)
        self.var_table.setColumnCount(3)
        cts.update_varlist()
        self.var_table.setRowCount(len(cts.varlist))

        self.updateglobvar_fn()

        # right side with the tabs
        self.tabs_widget = QWidget()

        self.layout = qtw.QHBoxLayout(self.tabs_widget)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.main_panel = QTabWidget()
        self.tab1 = calib_win.CalibWin(self)
        self.tab2 = Rabbit_win.RabbitWin(self)

        self.main_panel.addTab(self.tab1, "Energy calibration")
        self.main_panel.addTab(self.tab2, "RABBIT")

        self.main_panel.currentChanged.connect(self.onTabChange)

        self.layout.addWidget(self.main_panel)

    ''' This function updates the values of the variable displayed on the left of the window. The list of the variables 
    is in the glob_var.py file'''
    def updateglobvar_fn(self) -> None:
        self.var_table.clear()
        self.var_table.setColumnCount(3)
        self.var_table.setRowCount(len(cts.varlist))

        self.var_table.setItem(0, 0, QTableWidgetItem("name"))
        self.var_table.setItem(0, 1, QTableWidgetItem("dtype"))
        self.var_table.setItem(0, 2, QTableWidgetItem("value"))

        cts.update_varlist() # updating the content of the variable list
        for i in range(len(cts.varlist)):
            var = cts.varlist[i][1]
            varname = cts.varlist[i][0]
            self.var_table.setItem(i, 0, QTableWidgetItem(cts.varlist[i][0]))
            type = self.custom_type_fn(var)
            self.var_table.setItem(i, 1, QTableWidgetItem(type))
            # these 'if's are here to make the variables easily readable
            if type == "int" or type == "float":
                if varname == "Ip" or varname == "Vp" or varname == "elow" or varname == "ehigh" or varname == "dE":
                    self.var_table.setItem(i, 2, QTableWidgetItem("{:.2f}".format(var)))
                elif varname == "L" or varname == "1st harm" or varname == "steps_nm" or varname == "steps_fs" \
                        or varname == "stepsnb" or varname == "bandsnb":
                    self.var_table.setItem(i, 2, QTableWidgetItem("{}".format(int(var))))
                else:
                    self.var_table.setItem(i, 2, QTableWidgetItem("{:.2e}" .format(var)))
            if type ==  "np.array":
                self.var_table.setItem(i, 2, QTableWidgetItem(str(var.shape)))
            if type == "list":
                self.var_table.setItem(i, 2, QTableWidgetItem(str(len(var))))
            if type == "bool":
                self.var_table.setItem(i, 2, QTableWidgetItem(str(var)))
            if type == "string":
                self.var_table.setItem(i, 2, QTableWidgetItem(var))

        self.var_table.resizeColumnsToContents()
        w = 0
        for i in range(self.var_table.columnCount()):
            w = w + self.var_table.columnWidth(i)
        self.var_table.setFixedWidth(w+5)

    def custom_type_fn(self, o) -> str:
        type = "other"
        if isinstance(o, int):
            type = "int"
        if isinstance(o, float):
            type = "float"
        if isinstance(o, np.ndarray):
            type = "np.array"
        if isinstance(o, list):
            type = "list"
        if isinstance(o, bool):
            type = "bool"
        if isinstance(o, str):
            type = "string"
        return type

    ''' function called when double-clicking on a variable'''
    def var_table_lr(self, doubleClickedIndex) -> None:
        rowindex = doubleClickedIndex.row()
        try:
            vardiag = varDialog(rowindex, self) # new class created below
        except Exception:
            print(traceback.format_exception(*sys.exc_info()))

    ''' Updating elow, ehigh and dE when going from one tab to the other'''
    def onTabChange(self) -> None:

        self.tab1.elow_le.setText("{:.2f}".format(cts.elow))
        self.tab1.ehigh_le.setText("{:.2f}".format(cts.ehigh))
        self.tab1.dE_le.setText("{:.2f}".format(cts.dE))
        self.tab2.elow_le.setText("{:.2f}".format(cts.elow))
        self.tab2.ehigh_le.setText("{:.2f}".format(cts.ehigh))
        self.tab2.dE_le.setText("{:.2f}".format(cts.dE))

''' The object that popups when double-clicking on a "environment" variable'''
class varDialog(QDialog):
    def __init__(self, index, parent=mainWin):
        super(varDialog, self).__init__(parent)

        self.varname = cts.varlist[index][0]
        self.var = cts.varlist[index][1]
        self.par = self.parent()

        self.setGeometry(850, 500, 450, 300)
        self.setWindowTitle(self.varname)
        self.mainlayout = QHBoxLayout()

        self.table = QTableWidget()
        self.table.setSelectionBehavior(QTableView.SelectRows)

        if(isinstance(self.var, int) or isinstance(self.var, float)):
            self.table.setColumnCount(1)
            self.table.setRowCount(1)
            self.table.setItem(0, 0, QTableWidgetItem(str(self.var)))

        if isinstance(self.var, str):
            self.table.setColumnCount(1)
            self.table.setRowCount(1)
            self.table.setItem(0, 0, QTableWidgetItem(self.var))

        if(isinstance(self.var, np.ndarray)):
            rownb = self.var.shape[0]
            try:
                colnb = self.var.shape[1]
            except IndexError:
                colnb = 1
            self.table.setColumnCount(colnb)
            self.table.setRowCount(rownb)

            for i in range(rownb):
                if colnb == 1:
                    self.table.setItem(i, 0, QTableWidgetItem("{:.3e}".format(self.var[i])))
                else:
                    for j in range(colnb):
                        self.table.setItem(i, j, QTableWidgetItem("{:.3e}".format(self.var[i,j])))


        self.mainlayout.addWidget(self.table)
        self.setLayout(self.mainlayout)
        self.show()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    m = mainWin()
    sys.exit(app.exec())
(app.exec())