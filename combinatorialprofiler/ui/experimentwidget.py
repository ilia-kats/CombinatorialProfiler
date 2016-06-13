from pkg_resources import resource_stream

from PyQt5.QtWidgets import QWidget, QTableWidgetItem
from PyQt5.QtCore import QRegExp, Qt
from PyQt5.QtGui import QDoubleValidator, QRegExpValidator
from PyQt5 import uic

from .simpledelegate import SimpleDelegate

class ExperimentWidget(QWidget):
    ui = uic.loadUiType(resource_stream(__name__, "experiment.ui"))

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = self.__class__.ui[0]()
        self.ui.setupUi(self)
        self.ui.barcodes_fw.setLabel("Forward barcodes:")
        self.ui.barcodes_rev.setLabel("Reverse barcodes:")
        self.ui.named_inserts.setLabel("Named variable sequences:")

        self.cellsdelegate = SimpleDelegate(QDoubleValidator(0, 1, 1000), self)
        proto = QTableWidgetItem()
        proto.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEditable | Qt.ItemIsEnabled)
        proto.setData(Qt.EditRole, "0.0")
        self.ui.sortedCellsTbl.setItemPrototype(proto)
        self.ui.sortedCellsTbl.setItemDelegate(self.cellsdelegate)

        self.ui.insertSequence.setValidator(QRegExpValidator(QRegExp("^[atcg]+[n]+[atcg]+$", Qt.CaseInsensitive)))

        self.ui.barcodes_fw.rowAdded.connect(self.fwCodeAdded)
        self.ui.barcodes_fw.rowChanged.connect(self.fwCodeChanged)
        self.ui.barcodes_fw.rowRemoved.connect(self.fwCodeRemoved)

        self.ui.barcodes_rev.rowAdded.connect(self.revCodeAdded)
        self.ui.barcodes_rev.rowChanged.connect(self.revCodeChanged)
        self.ui.barcodes_rev.rowRemoved.connect(self.revCodeRemoved)

    def fwCodeAdded(self, row, text):
        nrow = self.ui.sortedCellsTbl.rowCount()
        if row != nrow:
            self.ui.sortedCellsTbl.insertRow(nrow)
            self.ui.sortedCellsTbl.setVerticalHeaderItem(nrow, QTableWidgetItem(text))

        if not self.ui.barcodes_rev.count() and not self.ui.sortedCellsTbl.columnCount():
            self.ui.sortedCellsTbl.insertColumn(0)
            self.ui.sortedCellsTbl.setHorizontalHeaderItem(0, QTableWidgetItem())

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

    def serialize(self):
        d = {'insert': self.ui.insertSequence.text(), 'barcodes_fw': self.ui.barcodes_fw.serialize(), 'barcodes_rev': self.ui.barcodes_rev.serialize()}
        if self.ui.named_inserts.count():
            d['named_inserts'] = self.ui.named_inserts.serialize()
        if self.ui.ndsiGrp.isChecked():
            if self.ui.ndsiOnFw.isChecked():
                d['ndsi'] = "forward"
            else:
                d['ndsi'] = "reverse"
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

    def unserialize(self, d):
        self.ui.sortedCellsTbl.clear()

        self.ui.insertSequence.setText(d.get('insert', ''))
        self.ui.barcodes_fw.unserialize(d.get('barcodes_fw', {}))
        self.ui.barcodes_rev.unserialize(d.get('barcodes_rev', {}))
        self.ui.named_inserts.unserialize(d.get('named_inserts', {}))

        fwrows = {}
        revcols = {}
        for r in range(self.ui.sortedCellsTbl.rowCount()):
            fwrows[self.ui.sortedCellsTbl.model().headerData(r, Qt.Vertical)] = r
        for c in range(self.ui.sortedCellsTbl.columnCount()):
            revcols[self.ui.sortedCellsTbl.model().headerData(c, Qt.Horizontal)] = c
        for fw, r in d.get('sortedcells', {}).items():
            for rev, f in r.items():
                item = QTableWidgetItem(self.ui.sortedCellsTbl.itemPrototype())
                item.setData(Qt.DisplayRole, str(f))
                self.ui.sortedCellsTbl.setItem(fwrows[fw], revcols[rev], item)

        if 'ndsi' in d:
            self.ui.ndsiGrp.setChecked(True)
            if d['ndsi'] == 'forward':
                self.ui.ndsiOnFw.setChecked(True)
            elif d['ndsi'] == 'reverse':
                self.ui.ndsiOnRev.setChecked(True)
