# -*- coding: utf-8 -*-
from pkg_resources import resource_stream

from PyQt5.QtWidgets import QWidget
from PyQt5.QtCore import pyqtSignal
from PyQt5 import uic

class BarcodeHammingWidget(QWidget):
    ui = uic.loadUiType(resource_stream(__name__, "barcodehammingwidget.ui"))

    changed = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)
        self.ui.barcodeLength.valueChanged.connect(self.changed)
        self.ui.barcodeMismatches.valueChanged.connect(self.changed)

    def serialize(self):
        return {"barcode_length" : self.ui.barcodeLength.value(), "barcode_mismatches" : round(self.ui.barcodeMismatches.value(), 3)}

    def unserialize(self, d):
        self.ui.barcodeLength.setValue(d.get('barcode_length', 0))
        self.ui.barcodeMismatches.setValue(d.get('barcode_mismatches', 0))
