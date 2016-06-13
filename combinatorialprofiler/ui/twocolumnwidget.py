#-*- coding: utf-8 -*-
from pkg_resources import resource_stream

from PyQt5.QtWidgets import QWidget, QTableWidgetItem
from PyQt5.QtCore import QRegExp, Qt, pyqtSignal
from PyQt5.QtGui import QRegExpValidator
from PyQt5 import uic

from .simpledelegate import SimpleDelegate

class TwoColumnWidget(QWidget):
    ui = uic.loadUiType(resource_stream(__name__, "twocolumnwidget.ui"))

    rowAdded = pyqtSignal(int, 'QString')
    rowRemoved = pyqtSignal(int)
    rowChanged = pyqtSignal(int, 'QString')

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)
        self.ui.seqTbl.resizeColumnsToContents()

        self.seqdelegate = SimpleDelegate(QRegExpValidator(QRegExp( "^[atcg]+$", Qt.CaseInsensitive)), self)

        self.ui.addBtn.clicked.connect(self.addSequence)
        self.ui.removeBtn.clicked.connect(self.removeSequences)
        self.ui.seqTbl.cellChanged.connect(self.cellChanged)
        proto = QTableWidgetItem()
        proto.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEditable | Qt.ItemIsEnabled)
        self.ui.seqTbl.setItemPrototype(proto)
        self.ui.seqTbl.setItemDelegateForColumn(1, self.seqdelegate)

    def setLabel(self, label):
        self.ui.label.setText(label)

    def addSequence(self, name=None, sequence=None):
        nameitem = QTableWidgetItem(self.ui.seqTbl.itemPrototype())
        seqitem = QTableWidgetItem(self.ui.seqTbl.itemPrototype())

        if name:
            nameitem.setText(name)
        if sequence:
            seqitem.setText(sequence)

        nrow = self.ui.seqTbl.rowCount()
        self.ui.seqTbl.insertRow(nrow)
        self.rowAdded.emit(nrow + 1, nameitem.text())
        self.ui.seqTbl.setItem(nrow, 0, nameitem)
        self.ui.seqTbl.setItem(nrow, 1, seqitem)
        self.ui.seqTbl.setCurrentItem(nameitem)
        if not name:
            self.ui.seqTbl.editItem(nameitem)

    def removeSequences(self):
        selected = self.ui.seqTbl.selectionModel().selectedRows()
        while len(selected):
            row = selected[-1].row()
            self.ui.seqTbl.removeRow(row)
            self.rowRemoved.emit(row)
            selected = self.ui.seqTbl.selectionModel().selectedRows()

    def cellChanged(self, row, col):
        if col == 0:
            self.rowChanged.emit(row, self.ui.seqTbl.item(row, col).text())

    def count(self):
        return self.ui.seqTbl.rowCount()

    def serialize(self):
        return {self.ui.seqTbl.item(i, 0).text() : self.ui.seqTbl.item(i,1).text() for i in range(self.ui.seqTbl.rowCount())}

    def unserialize(self, d):
        self.ui.seqTbl.clear()
        for k in sorted(d.keys()):
            self.addSequence(k, d[k])
        self.ui.seqTbl.resizeColumnsToContents()
