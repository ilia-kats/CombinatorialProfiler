# -*- coding: utf-8 -*-
from pkg_resources import resource_stream

from PyQt5.QtWidgets import QWidget, QListWidgetItem
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5 import uic

from .experimentwidget import ExperimentWidget

class ExperimentsWidget(QWidget):
    ui = uic.loadUiType(resource_stream(__name__, "experimentswidget.ui"))
    valid = pyqtSignal(bool)

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
        self.ui.listWidget.model().rowsRemoved.connect(self.experimentRemoved)
        self.ui.listWidget.itemChanged.connect(self.experimentNameChanged)

        self.evalid = {}

    def addExperiment(self, name=None, d=None):
        if not name:
            name = "New Experiment %d" % self.experiments
        item = QListWidgetItem(name, self.ui.listWidget)
        item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEditable | Qt.ItemIsEnabled)
        exp = ExperimentWidget(self.ui.stackedWidget)
        exp.valid.connect(self.experimentValid)
        self.ui.stackedWidget.addWidget(exp)
        self.ui.listWidget.setCurrentItem(item)
        if not d:
            self.experiments += 1
            self.ui.listWidget.editItem(item)
        else:
            exp.unserialize(d)
        exp.valid.emit(exp.isValid())

    def removeExperiment(self):
        r = self.ui.listWidget.currentRow()
        self.ui.listWidget.takeItem(r)

    def experimentRemoved(parent, first, last):
        for i in range(first, last + 1):
            self.ui.stackedWidget.removeWidget(self.ui.stackedWidget.widget(i))
        self.valid.emit(self.isValid())

    def experimentNameChanged(self):
        self.valid.emit(self.isValid())

    def serialize(self):
        return {self.ui.listWidget.item(i).text() : self.ui.stackedWidget.widget(i).serialize() for i in range(self.ui.listWidget.count())}

    def unserialize(self, d):
        self.ui.listWidget.clear()
        for n, e in d.items():
            self.addExperiment(n, e)

    def experimentValid(self, valid):
        self.evalid[self.sender()] = valid
        self.valid.emit(self.isValid())

    def isValid(self):
        seen = set()
        for i in range(self.ui.listWidget.count()):
            t = self.ui.listWidget.item(i).text()
            if not t or t in seen:
                return False
            seen.add(t)
        if not all(self.evalid.values()):
            return False
        return True
