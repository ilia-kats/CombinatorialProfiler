#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os.path
import sys
import json

from pkg_resources import resource_stream

from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QListWidgetItem, QTableWidgetItem, QDialogButtonBox, QFileDialog, QMessageBox
from PyQt5.QtCore import QRegExp, Qt
from PyQt5.QtGui import QDoubleValidator, QRegExpValidator
from PyQt5 import uic

from .experimentswidget import ExperimentsWidget
from .experimentwidget import ExperimentWidget

class MainWidget(QWidget):
    ui = uic.loadUiType(resource_stream(__name__, "main.ui"))

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)

        self.ui.buttonBox.button(QDialogButtonBox.Save).clicked.connect(self.saveClicked)
        self.ui.buttonBox.button(QDialogButtonBox.Open).clicked.connect(self.openClicked)
        self.ui.buttonBox.button(QDialogButtonBox.Close).clicked.connect(self.closeClicked)

    def saveClicked(self):
        dlg = QFileDialog(self)
        dlg.setDefaultSuffix("json")
        dlg.setFileMode(QFileDialog.AnyFile)
        dlg.setAcceptMode(QFileDialog.AcceptSave)
        dlg.setNameFilters(("JSON files (*.json)", "All files (*)"))
        if dlg.exec():
            path = dlg.selectedFiles()[0]
            with open(path, 'w') as f:
                d = self.ui.settingsTab.serialize()
                d['experiments'] = self.ui.experimentsTab.serialize()
                json.dump(d, f, indent=4)

    def openClicked(self):
        dlg = QFileDialog(self)
        dlg.setDefaultSuffix("json")
        dlg.setFileMode(QFileDialog.ExistingFile)
        dlg.setAcceptMode(QFileDialog.AcceptOpen)
        dlg.setNameFilters(("JSON files (*.json)", "All files (*)"))
        if dlg.exec():
            path = dlg.selectedFiles()[0]
            try:
                with open(path, 'r') as f:
                    d = json.load(f)
                    self.ui.settingsTab.unserialize(d)
                    self.ui.experimentsTab.unserialize(d['experiments'])
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def closeClicked(self):
        QApplication.exit()


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.widget = MainWidget(parent)
        self.setCentralWidget(self.widget)

def main():
    app = QApplication(sys.argv)
    main = MainWindow()
    main.resize(1024, 768)
    main.show()
    return app.exec()
