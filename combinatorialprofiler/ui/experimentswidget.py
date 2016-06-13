from pkg_resources import resource_stream

from PyQt5.QtWidgets import QWidget, QListWidgetItem
from PyQt5.QtCore import Qt
from PyQt5 import uic

from .experimentwidget import ExperimentWidget

class ExperimentsWidget(QWidget):
    ui = uic.loadUiType(resource_stream(__name__, "experimentswidget.ui"))

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

        self.ui.listWidget.model().rowsInserted.connect(self.experimentAdded)
        self.ui.listWidget.model().rowsRemoved.connect(self.experimentRemoved)

    def addExperiment(self, name=None):
        if not name:
            name = "New Experiment %d" % self.experiments
        item = QListWidgetItem(name, self.ui.listWidget)
        item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEditable | Qt.ItemIsEnabled)
        self.ui.listWidget.setCurrentItem(item)
        self.ui.listWidget.editItem(item)

    def experimentAdded(self, parent, first, last):
        self.experiments += 1
        self.ui.stackedWidget.addWidget(ExperimentWidget(self.ui.stackedWidget))

    def experimentRemoved(self, parent, first, last):
        self.ui.stackedWidget.removeWidget(self.ui.stackedWidget.widget(first))

    def removeExperiment(self):
        r = self.ui.listWidget.currentRow()

        item = self.ui.listWidget.takeItem(r)

    def serialize(self):
        return {self.ui.listWidget.item(i).text() : self.ui.stackedWidget.widget(i).serialize() for i in range(self.ui.listWidget.count())}

