#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os.path
import sys
import json
from distutils.version import StrictVersion

from pkg_resources import resource_stream

from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QFileDialog, QMessageBox, QStyle, QAction
from PyQt5.QtCore import pyqtRemoveInputHook, pyqtSignal
from PyQt5.QtGui import QIcon, QKeySequence
from PyQt5 import uic

from .experimentswidget import ExperimentsWidget
from .experimentwidget import ExperimentWidget
from .aboutdialog import AboutDialog
from .resources import resources
from .util import WaitCursor, TempDefaultCursor

from combinatorialprofiler import jsonversion
from . import progname

class MainWidget(QWidget):
    ui = uic.loadUiType(resource_stream(__name__, "main.ui"))
    currentFileChanged = pyqtSignal('QString')
    modified = pyqtSignal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)

        style = QApplication.style()
        self.saveAction = QAction(style.standardIcon(QStyle.SP_DialogSaveButton), "Save", self)
        self.openAction = QAction(style.standardIcon(QStyle.SP_DialogOpenButton), "Open", self)
        self.closeAction = QAction(style.standardIcon(QStyle.SP_DialogCloseButton), "Close", self)
        self.aboutAction = QAction(style.standardIcon(QStyle.SP_MessageBoxInformation), "About", self)
        self.saveAction.setEnabled(False)
        self.saveAction.triggered.connect(self.saveClicked)
        self.openAction.triggered.connect(self.openClicked)
        self.closeAction.triggered.connect(self.close)
        self.isModified = False
        self.aboutAction.triggered.connect(self.about)

        self.saveAction.setShortcut(QKeySequence.Save)
        self.openAction.setShortcut(QKeySequence.Open)
        self.closeAction.setShortcut(QKeySequence.Close)

        self.ui.experimentsTab.valid.connect(self.saveAction.setEnabled)
        self.ui.experimentsTab.valid.connect(self.changed)
        self.ui.settingsTab.changed.connect(self.changed)

        self.ui.saveBtn.setDefaultAction(self.saveAction)
        self.ui.openBtn.setDefaultAction(self.openAction)
        self.ui.closeBtn.setDefaultAction(self.closeAction)
        self.ui.aboutBtn.setDefaultAction(self.aboutAction)

    def changed(self):
        self.setModified(True)

    def setModified(self, modified=True):
        self.isModified = modified
        self.modified.emit(self.isModified)

    def canQuit(self):
        if self.isModified:
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
        AboutDialog(self).exec()

    def saveClicked(self):
        dlg = QFileDialog(self)
        dlg.setDefaultSuffix("json")
        dlg.setFileMode(QFileDialog.AnyFile)
        dlg.setAcceptMode(QFileDialog.AcceptSave)
        dlg.setNameFilters(("JSON files (*.json)", "All files (*)"))
        if dlg.exec():
            with WaitCursor():
                path = dlg.selectedFiles()[0]
                with open(path, 'w') as f:
                    d = self.ui.settingsTab.serialize()
                    d['experiments'] = self.ui.experimentsTab.serialize()
                    d['version'] = str(jsonversion)
                    json.dump(d, f, indent=4)
                self.setModified(False)
                self.currentFileChanged.emit(path)

    def openClicked(self):
        dlg = QFileDialog(self)
        dlg.setDefaultSuffix("json")
        dlg.setFileMode(QFileDialog.ExistingFile)
        dlg.setAcceptMode(QFileDialog.AcceptOpen)
        dlg.setNameFilters(("JSON files (*.json)", "All files (*)"))
        if dlg.exec():
            path = dlg.selectedFiles()[0]
            try:
                with WaitCursor():
                    with open(path, 'r') as f:
                        d = json.load(f)
                        warn = None
                        if 'version' not in d:
                            warn = "The JSON configuration file does not contain a version number. It may be incompatible with this program."
                        else:
                            cversion = StrictVersion(d['version'])
                            if cversion.version[0] != jsonversion.version[0]:
                                raise BaseException("The JSON configuration file format is not compatible with this version of %s" % QApplication.applicationDisplayName())
                            elif cversion.version[1] != jsonversion.version[1]:
                                warn = "The JSON configuration file format version does not match this version of %s. Incompatibilites should be handled gracefully, but unexpected results may occur." % QApplication.applicationDisplayName()
                        if warn:
                            with TempDefaultCursor():
                                QMessageBox.warning(self, "Warning", warn)
                        self.ui.settingsTab.unserialize(d)
                        self.ui.experimentsTab.unserialize(d['experiments'])
                    self.setModified(False)
                    self.currentFileChanged.emit(path)
            except BaseException as e:
                QMessageBox.critical(self, "Error", str(e))


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        w = MainWidget(parent)
        self.setCentralWidget(w)
        w.modified.connect(self.setWindowModified)
        w.currentFileChanged.connect(self.fileChanged)
        self.fileChanged()

    def fileChanged(self, path=""):
        self.setWindowTitle(("%s[*] - " % os.path.basename(path)) + QApplication.applicationDisplayName())

    def closeEvent(self, e):
        if not self.centralWidget().canQuit():
            e.ignore()
        else:
            super().closeEvent(e)

def main():
    pyqtRemoveInputHook()
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon(":/icon.svg"))
    app.setApplicationName(progname)
    main = MainWindow()
    main.resize(1280, 1024)
    main.show()
    return app.exec()
