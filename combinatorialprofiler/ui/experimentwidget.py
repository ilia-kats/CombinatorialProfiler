#-*- coding: utf-8 -*-
from io import StringIO
from pkg_resources import resource_stream

from PyQt5.QtWidgets import QApplication, QWidget, QTableWidgetItem, QStyle, QMessageBox, QFileDialog
from PyQt5.QtCore import QRegExp, Qt, pyqtSignal
from PyQt5.QtGui import QDoubleValidator, QRegExpValidator, QIcon
from PyQt5 import uic

from .simpledelegate import SimpleDelegate
from ..csvio import readCellCounts

class ExperimentWidget(QWidget):
    ui = uic.loadUiType(resource_stream(__name__, "experiment.ui"))
    valid = pyqtSignal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)
        self.ui.barcodes_fw.setLabel("Forward barcodes:")
        self.ui.barcodes_rev.setLabel("Reverse barcodes:")
        self.ui.named_inserts.setLabel("Named variable sequences:")

        style = QApplication.style()
        self.ui.fromFileBtn.setIcon(style.standardIcon(QStyle.SP_DirOpenIcon, None, self.ui.fromFileBtn))
        if QIcon.hasThemeIcon("edit-paste"):
            self.ui.fromClipboardBtn.setIcon(QIcon.fromTheme("edit-paste"))
        else:
            self.ui.fromClipboardBtn.setText("Paste")

        self.cellsdelegate = SimpleDelegate(QDoubleValidator(0, 1, 1000), self)
        proto = QTableWidgetItem()
        proto.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEditable | Qt.ItemIsEnabled)
        proto.setData(Qt.EditRole, "0.0")
        self.ui.sortedCellsTbl.setItemPrototype(proto)
        self.ui.sortedCellsTbl.setItemDelegate(self.cellsdelegate)
        self.ui.sortedCellsTbl.cellChanged.connect(self._validChanged)

        self.ui.insertSequence.setValidator(QRegExpValidator(QRegExp("^[atcg]+[n]+[atcg]+$", Qt.CaseInsensitive)))

        self.ui.fromFileBtn.clicked.connect(self.fromFile)
        self.ui.fromClipboardBtn.clicked.connect(self.fromClipboard)

        self.ui.barcodes_fw.rowAdded.connect(self.fwCodeAdded)
        self.ui.barcodes_fw.rowChanged.connect(self.fwCodeChanged)
        self.ui.barcodes_fw.rowRemoved.connect(self.fwCodeRemoved)

        self.ui.barcodes_rev.rowAdded.connect(self.revCodeAdded)
        self.ui.barcodes_rev.rowChanged.connect(self.revCodeChanged)
        self.ui.barcodes_rev.rowRemoved.connect(self.revCodeRemoved)

        self.ui.barcodes_fw.valid.connect(self._validChanged)
        self.ui.barcodes_rev.valid.connect(self._validChanged)
        self.ui.named_inserts.valid.connect(self._validChanged)
        self.ui.insertSequence.textChanged.connect(self._validChanged)
        self.ui.dsiGrp.toggled.connect(self._validChanged)
        self.ui.dsiOnFw.toggled.connect(self._validChanged)
        self.ui.dsiOnRev.toggled.connect(self._validChanged)

    def fwCodeAdded(self, row, text):
        nrow = self.ui.sortedCellsTbl.rowCount()
        if row != nrow:
            self.ui.sortedCellsTbl.insertRow(nrow)
            self.ui.sortedCellsTbl.setVerticalHeaderItem(nrow, QTableWidgetItem(text))

        if not self.ui.barcodes_rev.count() and not self.ui.sortedCellsTbl.columnCount():
            self.ui.sortedCellsTbl.insertColumn(0)
            self.ui.sortedCellsTbl.setHorizontalHeaderItem(0, QTableWidgetItem())
        self.ui.sortedCellsTbl.resizeColumnsToContents()

    def fwCodeChanged(self, row, text):
        self.ui.sortedCellsTbl.verticalHeaderItem(row).setText(text)

    def fwCodeRemoved(self, row):
        nrow = self.ui.sortedCellsTbl.rowCount()
        revcount = self.ui.barcodes_rev.count()
        if nrow > 1 or not revcount:
            self.ui.sortedCellsTbl.removeRow(row)
            if nrow == 1:
                self.ui.sortedCellsTbl.removeColumn(0)
        elif nrow == 1 and revcount:
            self.ui.sortedCellsTbl.verticalHeaderItem(0).setText('')

    def revCodeAdded(self, column, text):
        ncol = self.ui.sortedCellsTbl.columnCount()
        if column != ncol:
            self.ui.sortedCellsTbl.insertColumn(ncol)
            self.ui.sortedCellsTbl.setHorizontalHeaderItem(ncol, QTableWidgetItem(text))

        if not self.ui.barcodes_fw.count() and not self.ui.sortedCellsTbl.rowCount():
            self.ui.sortedCellsTbl.insertRow(0)
            self.ui.sortedCellsTbl.setVerticalHeaderItem(0, QTableWidgetItem())
        self.ui.sortedCellsTbl.resizeColumnsToContents()

    def revCodeChanged(self, column, text):
        self.ui.sortedCellsTbl.horizontalHeaderItem(column).setText(text)

    def revCodeRemoved(self, column):
        ncol = self.ui.sortedCellsTbl.columnCount()
        fwcount = self.ui.barcodes_fw.count()
        if ncol > 1 or not fwcount:
            self.ui.sortedCellsTbl.removeColumn(column)
            if ncol == 1:
                self.ui.sortedCellsTbl.removeRow(0)
        elif ncol == 1 and fwcount:
            self.ui.sortedCellsTbl.horizontalHeaderItem(0).setText('')

    def _validChanged(self):
        self.valid.emit(self.isValid())

    def _fromFile(self, f):
        try:
            fwrows, revcols = self.getCellHeaderMapping()
            cells = readCellCounts(f, fwrows.keys(), revcols.keys())
            for ((fw, rev), val) in cells.groupby(level=('barcode_fw', 'barcode_rev')):
                if fw in fwrows and rev in revcols:
                    y = fwrows[fw]
                    x = revcols[rev]
                elif fw in fwrows and not self.ui.barcodes_rev.count():
                    y = fwrows[fw]
                    x = 0
                elif not self.ui.barcodes_fw.count() and rev in revcols:
                    y = 0
                    x = revcols[rev]
                else:
                    continue
                item = QTableWidgetItem(self.ui.sortedCellsTbl.itemPrototype())
                item.setData(Qt.DisplayRole, val.iloc[0].values.squeeze().item())
                self.ui.sortedCellsTbl.setItem(y,x, item)
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

    def serialize(self):
        d = {'insert': self.ui.insertSequence.text(), 'barcodes_fw': self.ui.barcodes_fw.serialize(), 'barcodes_rev': self.ui.barcodes_rev.serialize()}
        if self.ui.named_inserts.count():
            d['named_inserts'] = self.ui.named_inserts.serialize()
        if self.ui.dsiGrp.isChecked():
            if self.ui.dsiOnFw.isChecked():
                d['dsi'] = "forward"
            else:
                d['dsi'] = "reverse"
            cells = {}
            for r in range(self.ui.sortedCellsTbl.rowCount()):
                fwcode = self.ui.sortedCellsTbl.verticalHeaderItem(r).text()
                fwdict = {}
                for c in range(self.ui.sortedCellsTbl.columnCount()):
                    revcode = self.ui.sortedCellsTbl.horizontalHeaderItem(c).text()
                    item = self.ui.sortedCellsTbl.item(r, c)
                    if not item:
                        item = self.ui.sortedCellsTbl.itemPrototype()
                    fwdict[revcode] = float(item.data(Qt.EditRole))
                cells[fwcode] = fwdict
            d['sortedcells'] = cells
        return d

    def getCellHeaderMapping(self):
        fwrows = {}
        revcols = {}
        for r in range(self.ui.sortedCellsTbl.rowCount()):
            fwrows[self.ui.sortedCellsTbl.model().headerData(r, Qt.Vertical)] = r
        for c in range(self.ui.sortedCellsTbl.columnCount()):
            revcols[self.ui.sortedCellsTbl.model().headerData(c, Qt.Horizontal)] = c
        return (fwrows, revcols)

    def unserialize(self, d):
        self.ui.sortedCellsTbl.clear()

        self.ui.insertSequence.setText(d.get('insert', ''))
        self.ui.barcodes_fw.unserialize(d.get('barcodes_fw', {}))
        self.ui.barcodes_rev.unserialize(d.get('barcodes_rev', {}))
        self.ui.named_inserts.unserialize(d.get('named_inserts', {}))

        fwrows, revcols = self.getCellHeaderMapping()
        for fw, r in d.get('sortedcells', {}).items():
            for rev, f in r.items():
                item = QTableWidgetItem(self.ui.sortedCellsTbl.itemPrototype())
                item.setData(Qt.DisplayRole, str(f))
                self.ui.sortedCellsTbl.setItem(fwrows[fw], revcols[rev], item)
        self.ui.sortedCellsTbl.resizeColumnsToContents()

        if 'dsi' in d:
            self.ui.dsiGrp.setChecked(True)
            if d['dsi'] == 'forward':
                self.ui.dsiOnFw.setChecked(True)
            elif d['dsi'] == 'reverse':
                self.ui.dsiOnRev.setChecked(True)

    def isValid(self):
        if not self.ui.barcodes_fw.isValid():
            return False
        if not self.ui.barcodes_rev.isValid():
            return False
        if not self.ui.named_inserts.isValid():
            return False
        if not self.ui.insertSequence.hasAcceptableInput():
            return False
        if self.ui.dsiGrp.isChecked():
            if self.ui.dsiOnFw.isChecked():
                oloop = self.ui.sortedCellsTbl.columnCount
                iloop = self.ui.sortedCellsTbl.rowCount
                coords = lambda o,i: (i, o)
            else:
                oloop = self.ui.sortedCellsTbl.rowCount
                iloop = self.ui.sortedCellsTbl.columnCount
                coords = lambda o,i: (o, i)
            for o in range(oloop()):
                csum = 0.0
                for i in range(iloop()):
                    item = self.ui.sortedCellsTbl.item(*coords(o,i))
                    if item:
                        csum += float(item.text())
                if round(csum, 1) != 1:
                    return False
        return True
