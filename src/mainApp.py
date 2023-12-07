import sys
import os
from PyQt6.QtWidgets import (
    QApplication,
    QLabel,
    QWidget,
    QLineEdit,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QFileDialog,
    QComboBox,
    QSpinBox,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QPixmap

from .widgets.settingWidget import SettingWidget
from .widgets.propertyWidget import PropertyWidget


class MainApp(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        self.setWindowTitle("Molecular Visualizer")
        self.resize(600, 400)

        layout = QVBoxLayout()
        self.setLayout(layout)

        setting_menu = SettingWidget()

        layout.addWidget(setting_menu)

        property_menu = PropertyWidget(setting_menu.load_path_label)
        layout.addWidget(property_menu)


if __name__ == "__main__":
    # Create the application and run it
    app = QApplication(sys.argv)
    ex = MainApp()
    ex.show()
    sys.exit(app.exec())
