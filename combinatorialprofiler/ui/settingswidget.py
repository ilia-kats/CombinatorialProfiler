# -*- coding: utf-8 -*-
from pkg_resources import resource_stream

from PyQt5.QtWidgets import QWidget
from PyQt5.QtCore import pyqtSignal
from PyQt5 import uic

from enum import Enum

from .columnresizer import ColumnResizer

class SettingsWidget(QWidget):
    ui = uic.loadUiType(resource_stream(__name__, "settingswidget.ui"))

    changed = pyqtSignal()

    class BarcodeAlgo(Enum):
        hamming = 0
        seqlev = 1

        def __str__(self):
            return self.name

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)
        self.ui.mismatches.valueChanged.connect(self.changed)
        self.ui.hamming.changed.connect(self.changed)
        self.ui.seqlev.changed.connect(self.changed)

        self.ui.barcodeAlgoSettings.setCurrentIndex(self.ui.barcodeAlgo.currentIndex())

        self.resizer = ColumnResizer(self)
        self.resizer.addWidgetsFromLayout(self.ui.top_formLayout, 0)
        self.resizer.addWidgetsFromLayout(self.ui.bottom_formLayout, 0)
        self.resizer.addWidgetsFromLayout(self.ui.hamming.layout(), 0)
        self.resizer.addWidgetsFromLayout(self.ui.seqlev.layout(), 0)

    def serialize(self):
        d = {"insert_mismatches" : self.ui.mismatches.value(), "barcode_match_algo" : str(self.BarcodeAlgo(self.ui.barcodeAlgo.currentIndex())), "profile_plot_min_count" : self.ui.plot_minCount.value()}
        d.update(self.ui.barcodeAlgoSettings.currentWidget().serialize())
        return d

    def unserialize(self, d):
        self.ui.mismatches.setValue(d.get('insert_mismatches', 0))
        self.ui.barcodeAlgo.setCurrentIndex(self.BarcodeAlgo[d.get('barcode_match_algo', 'hamming')].value)
        self.ui.barcodeAlgoSettings.currentWidget().unserialize(d)
        self.ui.plot_minCount.setValue(d.get('profile_plot_min_count', 2))
