from PyQt5.QtWidgets import QStyledItemDelegate, QLineEdit
from PyQt5.QtCore import Qt

class SimpleDelegate(QStyledItemDelegate):
    def __init__(self, validator=None, parent=None):
        super().__init__(parent)
        self.validator = validator

    def createEditor(self, parent, option, index):
        line = QLineEdit(parent)
        line.setValidator(self.validator)
        return line

    def setEditorData(self, editor, index):
        editor.setText(index.data(Qt.EditRole))

    def setModelData(self, editor, model, index):
        model.setData(index, editor.text(), Qt.EditRole)
