# PyQt5 port of ColumnResizer: https://github.com/agateau/columnresizer
# Original C++ code licensed as LGPL 2.1

from PyQt5.QtCore import Qt, QEvent, QObject, QTimer, pyqtSignal
from PyQt5.QtWidgets import QFormLayout, QGridLayout, QWidgetItem

class FormLayoutWidgetItem(QWidgetItem):
    def __init__( self, widget, formLayout, itemRole ):
        QWidgetItem.__init__(self, widget)
        self.width = -1
        self.layout = formLayout
        self.itemRole = itemRole

    def sizeHint(self):
        size = QWidgetItem.sizeHint(self)
        if self.width != -1:
            size.setWidth(self.width)
        return size

    def minimumSize(self):
        size = QWidgetItem.minimumSize(self)
        if self.width != -1:
            size.setWidth(self.width)
        return size

    def maximumSize(self):
        size = QWidgetItem.maximumSize(self)
        if self.width != -1:
            size.setWidth(self.width)
        return size

    def setWidth(self, width):
        if width != self.width:
            self.width = width
            self.invalidate()

    def setGeometry(self, _rect):
        rect = _rect
        width = self.widget().sizeHint().width()
        if self.itemRole == QFormLayout.LabelRole and self.layout.labelAlignment() & Qt.AlignRight:
            rect.setLeft(rect.right() - width)
        QWidgetItem.setGeometry(self, rect)

    def formLayout(self):
        return self.layout

class ColumnResizer(QObject):
    class ColumnResizerPrivate:
        widgets = []
        wrWidgetItemList = []
        gridColumnInfoList = []

        def __init__( self, q ):
            self.q = q
            self.updateTimer = QTimer(q)
            self.updateTimer.setSingleShot(True)
            self.updateTimer.setInterval(0)
            self.updateTimer.timeout.connect(q.updateWidth)

        def scheduleWithUpdate(self):
            self.updateTimer.start()

    def __init__(self, parent = None):
        QObject.__init__(self, parent)
        self.d = self.ColumnResizerPrivate(self)

    def addWidget(self, widget):
        self.d.widgets.append( widget )
        widget.installEventFilter(self)
        self.d.scheduleWithUpdate()

    def updateWidth(self):
        width = 0
        for widget in self.d.widgets:
            width = max(widget.sizeHint().width(), width)
        for item in self.d.wrWidgetItemList:
            item.setWidth(width)
            item.formLayout().update()
        for (l, i) in self.d.gridColumnInfoList:
            l.setColumnMinimumWidth(i, width)

    def eventFilter(self, obj, event):
        if event.type() == QEvent.Resize:
            self.d.scheduleWithUpdate()
        return False

    def addWidgetsFromLayout(self, layout, column):
        assert(column >= 0)
        if isinstance(layout, QGridLayout):
            self.addWidgetsFromGridLayout(layout, column)
        elif isinstance(layout, QFormLayout):
            if column <= QFormLayout.SpanningRole:
                self.addWidgetsFromFormLayout(layout, column)

    def addWidgetsFromGridLayout(self, layout, column):
        for row in range( layout.rowCount() ):
            item = layout.itemAtPosition(row, column)
            if item is None:
                continue
            widget = item.widget()
            if widget is None:
                continue
            self.addWidget( widget )
        self.d.gridColumnInfoList.append((layout, column))

    def addWidgetsFromFormLayout(self, layout, role):
        for row in range(layout.rowCount()):
            item = layout.itemAt(row, role)
            if item is None:
                continue
            widget = item.widget()
            if widget is None:
                continue
            layout.removeItem(item)
            del item
            newItem = FormLayoutWidgetItem(widget, layout, role)
            layout.setItem(row, role, newItem)
            self.addWidget(widget)
            self.d.wrWidgetItemList.append(newItem)
