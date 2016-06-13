from pkg_resources import resource_stream

from PyQt5.QtWidgets import QWidget
from PyQt5 import uic

class SettingsWidget(QWidget):
    ui = uic.loadUiType(resource_stream(__name__, "settingswidget.ui"))

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)

    def serialize(self):
        return {"insert_mismatches" : self.ui.mismatches.value(), "barcode_length" : self.ui.barcodeLength.value()}
