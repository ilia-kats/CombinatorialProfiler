#-*- coding: utf-8 -*-

from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QCursor

class WaitCursor:
    def __init__(self, switchCursor=True):
        self._switch = switchCursor

    def __enter__(self):
        if self._switch:
            QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))

    def __exit__(self, *args):
        if self._switch:
            QApplication.restoreOverrideCursor()
