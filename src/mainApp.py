import sys
from PyQt6.QtWidgets import (
    QApplication,
    QMainWindow,
    QLabel,
    QWidget,
    QVBoxLayout,
    QPushButton,
    QFileDialog,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QPixmap


class App(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecular Visualizer")
        self.resize(500, 350)

        layout = QVBoxLayout()
        self.setLayout(layout)

        button = QPushButton("Create Pictures!", clicked=self.getfile)
        layout.addWidget(button)

        # l1 = QLabel()
        # l1.setPixmap(QPixmap("python.jpg"))

    def getfile(self):
        fname = QFileDialog.getOpenFileName(
            self, "Open file", "/", "Image files (*.jpg *.gif)"
        )
        print(fname)


if __name__ == "__main__":
    # Create the application and run it
    app = QApplication(sys.argv)
    ex = App()
    ex.show()
    sys.exit(app.exec())
