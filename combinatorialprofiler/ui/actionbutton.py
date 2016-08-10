#-*- coding: utf-8 -*-

from PyQt5.QtWidgets import QPushButton, QAction

class ActionButton(QPushButton):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._defaultAction = None

    def setDefaultAction(self, action):
        if self._defaultAction:
            self._defaultAction.changed.disconnect(self.updateStatus)
            self.clicked.disconnect(self._defaultAction.trigger)
            self.removeAction(self._defaultAction)

        self._defaultAction = action
        if self._defaultAction:
            self.addAction(self._defaultAction)
            self.updateStatus()
            self._defaultAction.changed.connect(self.updateStatus)
            self.clicked.connect(self._defaultAction.trigger)

    def defaultAction(self):
        return self._defaultAction

    def updateStatus(self):
        self.setText(self._defaultAction.text())
        self.setStatusTip(self._defaultAction.statusTip())
        self.setToolTip(self._defaultAction.toolTip())
        self.setIcon(self._defaultAction.icon())
        self.setEnabled(self._defaultAction.isEnabled())
        self.setCheckable(self._defaultAction.isCheckable())
        self.setChecked(self._defaultAction.isChecked())
