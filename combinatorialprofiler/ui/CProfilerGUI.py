#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os.path
import sys
import json

from pkg_resources import resource_stream

from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QListWidgetItem, QTableWidgetItem, QDialogButtonBox, QFileDialog, QMessageBox, QStyle
from PyQt5.QtCore import QRegExp, Qt, pyqtRemoveInputHook
from PyQt5.QtGui import QDoubleValidator, QRegExpValidator
from PyQt5 import uic

from .experimentswidget import ExperimentsWidget
from .experimentwidget import ExperimentWidget

from combinatorialprofiler import version

class MainWidget(QWidget):
    ui = uic.loadUiType(resource_stream(__name__, "main.ui"))

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)

        style = QApplication.style()
        self.ui.saveBtn.setIcon(style.standardIcon(QStyle.SP_DialogSaveButton))
        self.ui.openBtn.setIcon(style.standardIcon(QStyle.SP_DialogOpenButton))
        self.ui.closeBtn.setIcon(style.standardIcon(QStyle.SP_DialogCloseButton))
        self.ui.aboutBtn.setIcon(style.standardIcon(QStyle.SP_MessageBoxInformation))

        self.ui.saveBtn.setEnabled(False)
        self.ui.saveBtn.clicked.connect(self.saveClicked)
        self.ui.openBtn.clicked.connect(self.openClicked)
        self.ui.closeBtn.clicked.connect(self.close)
        self.saved = True

        self.ui.aboutBtn.clicked.connect(self.about)

        self.ui.experimentsTab.valid.connect(self.ui.saveBtn.setEnabled)
        self.ui.experimentsTab.valid.connect(self.changed)
        self.ui.settingsTab.changed.connect(self.changed)

    def changed(self):
        self.saved = False

    def canQuit(self):
        if not self.saved:
            btn = QMessageBox.question(self, "Close now?", "The current configuration has not been saved. Quit anyway?", defaultButton = QMessageBox.No)
            if btn == QMessageBox.No:
                return False
            else:
                return True
        else:
            return True

    def close(self):
        QApplication.activeWindow().close()

    def about(self):
        QMessageBox.information(self, "About", """%s version %s.""" % (QApplication.applicationName(), version))

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
            self.saved = True

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
                self.saved = True
            except BaseException as e:
                QMessageBox.critical(self, "Error", str(e))


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.widget = MainWidget(parent)
        self.setCentralWidget(self.widget)

    def closeEvent(self, e):
        if not self.widget.canQuit():
            e.ignore()
        else:
            super().closeEvent(e)

def main():
    pyqtRemoveInputHook()
    app = QApplication(sys.argv)
    main = MainWindow()
    main.resize(1280, 1024)
    main.show()
    return app.exec()
