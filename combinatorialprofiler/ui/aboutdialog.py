# -*- coding: utf-8 -*-
from pkg_resources import resource_stream

from PyQt5.QtWidgets import QApplication, QWidget, QDialog, QStyle
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtCore import Qt, QSize, pyqtSignal
from PyQt5 import uic

from .experimentwidget import ExperimentWidget
from .resources import resources

from combinatorialprofiler import version

class AboutDialog(QDialog):
    ui = uic.loadUiType(resource_stream(__name__, "aboutdialog.ui"))

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)
        self.setWindowTitle("About")
        style = QApplication.style()

        self.ui.okBtn.clicked.connect(self.accept)
        self.ui.okBtn.setIcon(style.standardIcon(QStyle.SP_DialogOkButton))
        self.ui.aboutText.setText("%s version %s." % (QApplication.applicationDisplayName(), version))

        logosize = self.ui.aboutText.fontMetrics().width(self.ui.aboutText.text())
        self.ui.knoplabLogo.setPixmap(QIcon(":/knoplablogo.svg").pixmap(logosize, logosize))
