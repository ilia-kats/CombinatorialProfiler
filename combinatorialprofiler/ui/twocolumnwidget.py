#-*- coding: utf-8 -*-
from io import StringIO

from llist import dllist

from pkg_resources import resource_stream

from PyQt5.QtWidgets import QApplication, QWidget, QTableWidgetItem, QStyle, QMessageBox, QFileDialog
from PyQt5.QtCore import QRegExp, Qt, pyqtSignal
from PyQt5.QtGui import QRegExpValidator, QIcon, QCursor
from PyQt5 import uic

from .simpledelegate import SimpleDelegate
from .util import WaitCursor
from ..csvio import readBarcodes

class TwoColumnWidget(QWidget):
    ui = uic.loadUiType(resource_stream(__name__, "twocolumnwidget.ui"))

    rowAdded = pyqtSignal(int, 'QString')
    rowRemoved = pyqtSignal(int)
    rowChanged = pyqtSignal(int, 'QString')
    valid = pyqtSignal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.invalid = [0, 0]
        self.byrow = (dllist(), dllist())
        self.unique = ({}, {})
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)
        style = QApplication.style()
        self.ui.fromFileBtn.setIcon(style.standardIcon(QStyle.SP_DirOpenIcon, None, self.ui.fromFileBtn))
        if QIcon.hasThemeIcon("edit-paste"):
            self.ui.fromClipboardBtn.setIcon(QIcon.fromTheme("edit-paste"))
        else:
            self.ui.fromClipboardBtn.setText("Paste")
        self.ui.seqTbl.resizeColumnsToContents()

        self.seqdelegate = SimpleDelegate(QRegExpValidator(QRegExp( "^[atcg]+$", Qt.CaseInsensitive)), self)

        self.ui.addBtn.clicked.connect(self.addSequence)
        self.ui.removeBtn.clicked.connect(self.removeSequences)
        self.ui.fromFileBtn.clicked.connect(self.fromFile)
        self.ui.fromClipboardBtn.clicked.connect(self.fromClipboard)
        self.ui.seqTbl.cellChanged.connect(self.cellChanged)
        proto = QTableWidgetItem()
        proto.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEditable | Qt.ItemIsEnabled)
        self.ui.seqTbl.setItemPrototype(proto)
        self.ui.seqTbl.setItemDelegateForColumn(1, self.seqdelegate)

        self.ui.seqTbl.model().rowsInserted.connect(self.rowsChanged)
        self.ui.seqTbl.model().rowsRemoved.connect(self.rowsChanged)

        self.rowsChanged()

    def rowsChanged(self):
        self.ui.removeBtn.setEnabled(self.ui.seqTbl.model().rowCount() > 0)

    def setLabel(self, label):
        self.ui.label.setText(label)

    def addSequence(self, name=None, sequence=None):
        nameitem = QTableWidgetItem(self.ui.seqTbl.itemPrototype())
        seqitem = QTableWidgetItem(self.ui.seqTbl.itemPrototype())
        nrow = self.ui.seqTbl.rowCount()
        if name:
            nameitem.setText(name)
        if sequence:
            seqitem.setText(sequence)
        n = nameitem.text()
        s = seqitem.text()
        self.byrow[0].append(n)
        self.byrow[1].append(s)
        self.unique[0][n] = self.unique[0].get(n, 0) + 1
        self.unique[1][s] = self.unique[1].get(s, 0) + 1

        if self.unique[0][n] > 1 or not name:
            self.invalid[0] += 1
        if self.unique[1][s] > 1 or not sequence:
            self.invalid[1] += 1
        self.ui.seqTbl.insertRow(nrow)
        self.rowAdded.emit(nrow + 1, nameitem.text())
        self.ui.seqTbl.setItem(nrow, 0, nameitem)
        self.ui.seqTbl.setItem(nrow, 1, seqitem)
        self.ui.seqTbl.setCurrentItem(nameitem)
        if not name:
            self.ui.seqTbl.editItem(nameitem)
        self.valid.emit(self.isValid())

    def _fromFile(self, f):
        try:
            with WaitCursor():
                for n, c in readBarcodes(f):
                    self.addSequence(n, c)
                self.ui.seqTbl.resizeColumnsToContents()
        except BaseException as e:
            QMessageBox.critical(self, "Error", str(e))

    def fromFile(self):
        dlg = QFileDialog(self)
        dlg.setFileMode(QFileDialog.ExistingFile)
        dlg.setAcceptMode(QFileDialog.AcceptOpen)
        dlg.setNameFilters(("CSV files (*.csv *.txt)", "All files (*)"))
        if dlg.exec():
            path = dlg.selectedFiles()[0]
            self._fromFile(path)

    def fromClipboard(self):
        f = StringIO(QApplication.clipboard().text())
        self._fromFile(f)

    def removeSequences(self):
        selected = self.ui.seqTbl.selectionModel().selection()
        with WaitCursor():
            for i in range(len(selected)):
                row = selected[i].bottom()
                currnode = [self.byrow[col].nodeat(row) for col in range(2)]
                prevnode = []
                for r in range(row, selected[i].top() - 1, -1):
                    for col in range(2):
                        prevnode.append(currnode[col].prev)
                        t = self.byrow[col].remove(currnode[col])
                        self.unique[col][t] -= 1
                        count = self.unique[col][t]
                        if count > 0 or not t:
                            self.invalid[col] -= 1
                        if not count:
                            del self.unique[col][t]
                    currnode = prevnode
                    prevnode = []
                    self.ui.seqTbl.removeRow(r)
                    self.rowRemoved.emit(r)
            self.valid.emit(self.isValid())

    def cellChanged(self, row, col):
        ot = self.byrow[col][row]
        nt = self.ui.seqTbl.item(row, col).text()
        if col == 0:
            self.rowChanged.emit(row, nt)
        self.unique[col][ot] -= 1
        c = self.unique[col].get(nt, 0) + 1
        self.unique[col][nt] = c
        oldcount = self.unique[col][ot]
        if c == 1 and (oldcount > 0 and ot != nt or oldcount == 0 and not ot and nt):
            self.invalid[col] -= 1
        if not oldcount:
            del self.unique[col][self.byrow[col][row]]
        self.byrow[col][row] = nt
        self.valid.emit(self.isValid())

    def count(self):
        return self.ui.seqTbl.rowCount()

    def serialize(self):
        return {self.ui.seqTbl.item(i, 0).text() : self.ui.seqTbl.item(i,1).text() for i in range(self.ui.seqTbl.rowCount())}

    def unserialize(self, d):
        self.ui.seqTbl.clearContents()
        for k in sorted(d.keys()):
            self.addSequence(k, d[k])
        self.ui.seqTbl.resizeColumnsToContents()

    def isValid(self):
        return not any(self.invalid)
