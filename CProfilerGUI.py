#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os.path
import sys

from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QListWidgetItem, QTableWidgetItem
from PyQt5.QtCore import QRegExp, Qt
from PyQt5.QtGui import QDoubleValidator, QRegExpValidator
from PyQt5 import uic

uidir = os.path.join(os.path.dirname(__file__), "ui")
sys.path.append(uidir)

from simpledelegate import SimpleDelegate

class ExperimentWidget(QWidget):
    uifile = os.path.join(uidir, "experiment.ui")
    ui = uic.loadUiType(uifile)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)
        self.ui.barcodes_fw.setLabel("Forward barcodes:")
        self.ui.barcodes_rev.setLabel("Reverse barcodes:")
        self.ui.named_inserts.setLabel("Named variable sequences:")

        self.cellsdelegate = SimpleDelegate(QDoubleValidator(0, 1, 1000), self)
        self.ui.sortedCellsTbl.setItemDelegate(self.cellsdelegate)

        self.ui.insertSequence.setValidator(QRegExpValidator(QRegExp("^[atcg]+[n]+[atcg]+$", Qt.CaseInsensitive)))

        self.ui.barcodes_fw.rowAdded.connect(self.fwCodeAdded)
        self.ui.barcodes_fw.rowChanged.connect(self.fwCodeChanged)
        self.ui.barcodes_fw.rowRemoved.connect(self.fwCodeRemoved)

        self.ui.barcodes_rev.rowAdded.connect(self.revCodeAdded)
        self.ui.barcodes_rev.rowChanged.connect(self.revCodeChanged)
        self.ui.barcodes_rev.rowRemoved.connect(self.revCodeRemoved)

    def fwCodeAdded(self, row, text):
        nrow = self.ui.sortedCellsTbl.rowCount()
        if row != nrow:
            self.ui.sortedCellsTbl.insertRow(nrow)
            self.ui.sortedCellsTbl.setVerticalHeaderItem(nrow, QTableWidgetItem(text))

        if not self.ui.barcodes_rev.count() and not self.ui.sortedCellsTbl.columnCount():
            self.ui.sortedCellsTbl.insertColumn(0)
            self.ui.sortedCellsTbl.setHorizontalHeaderItem(0, QTableWidgetItem())

    def fwCodeChanged(self, row, text):
        self.ui.sortedCellsTbl.verticalHeaderItem(row).setText(text)

    def fwCodeRemoved(self, row):
        nrow = self.ui.sortedCellsTbl.rowCount()
        revcount = self.ui.barcodes_rev.count()
        if nrow > 1 or not revcount:
            self.ui.sortedCellsTbl.removeRow(row)
            if nrow == 1:
                self.ui.sortedCellsTbl.removeColumn(0)
        elif nrow == 1 and revcount:
            self.ui.sortedCellsTbl.verticalHeaderItem(0).setText('')

    def revCodeAdded(self, column, text):
        ncol = self.ui.sortedCellsTbl.columnCount()
        if column != ncol:
            self.ui.sortedCellsTbl.insertColumn(ncol)
            self.ui.sortedCellsTbl.setHorizontalHeaderItem(ncol, QTableWidgetItem(text))

        if not self.ui.barcodes_fw.count() and not self.ui.sortedCellsTbl.rowCount():
            self.ui.sortedCellsTbl.insertRow(0)
            self.ui.sortedCellsTbl.setVerticalHeaderItem(0, QTableWidgetItem())

    def revCodeChanged(self, column, text):
        self.ui.sortedCellsTbl.horizontalHeaderItem(column).setText(text)

    def revCodeRemoved(self, column):
        ncol = self.ui.sortedCellsTbl.columnCount()
        fwcount = self.ui.barcodes_fw.count()
        print(ncol, fwcount)
        if ncol > 1 or not fwcount:
            self.ui.sortedCellsTbl.removeColumn(column)
            if ncol == 1:
                self.ui.sortedCellsTbl.removeRow(0)
        elif ncol == 1 and fwcount:
            self.ui.sortedCellsTbl.horizontalHeaderItem(0).setText('')

class MainWidget(QWidget):
    uifile = os.path.join(uidir, "main.ui")
    ui = uic.loadUiType(uifile)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.experiments = 1
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)
        self.ui.splitter.setStretchFactor(0,0)
        self.ui.splitter.setStretchFactor(1,1)

        self.ui.addBtn.clicked.connect(self.addExperiment)
        self.ui.removeBtn.clicked.connect(self.removeExperiment)
        self.ui.listWidget.currentRowChanged.connect(self.ui.stackedWidget.setCurrentIndex)

    def addExperiment(self):
        item = QListWidgetItem("New Experiment %d" % (self.experiments), self.ui.listWidget)
        self.experiments += 1
        item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEditable | Qt.ItemIsEnabled)
        self.ui.stackedWidget.addWidget(ExperimentWidget(self.ui.stackedWidget))
        self.ui.listWidget.setCurrentItem(item)
        self.ui.listWidget.editItem(item)

    def removeExperiment(self):
        r = self.ui.listWidget.currentRow()
        self.ui.stackedWidget.removeWidget(self.ui.stackedWidget.currentWidget())
        item = self.ui.listWidget.takeItem(r)

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.widget = MainWidget(parent)
        self.setCentralWidget(self.widget)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main = MainWindow()
    main.resize(1024, 768)
    main.show()
    main.widget.ui.splitter.setSizes((10,1))
    app.exec()
