# -*- coding: utf-8 -*-
from pkg_resources import resource_stream

from PyQt5.QtWidgets import QWidget
from PyQt5.QtCore import pyqtSignal
from PyQt5 import uic

class SettingsWidget(QWidget):
    ui = uic.loadUiType(resource_stream(__name__, "settingswidget.ui"))

    changed = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)
        self.ui.mismatches.valueChanged.connect(self.changed)
        self.ui.barcodeLength.valueChanged.connect(self.changed)

    def serialize(self):
        return {"insert_mismatches" : self.ui.mismatches.value(), "barcode_length" : self.ui.barcodeLength.value()}

    def unserialize(self, d):
        self.ui.mismatches.setValue(d.get('insert_mismatches', 0))
        self.ui.barcodeLength.setValue(d.get('barcode_length', 0))
