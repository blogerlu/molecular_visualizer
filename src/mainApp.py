import sys
from PyQt6.QtWidgets import (
    QApplication,
    QMainWindow,
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


class App(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecular Visualizer")
        self.resize(400, 200)

        layout = QVBoxLayout()
        self.setLayout(layout)

        # load_path_layout
        self.load_path_layout = QHBoxLayout()
        self.load_path_label = QLineEdit("Path to smiles file...")
        self.load_path_layout.addWidget(self.load_path_label)

        button_load_path = QPushButton("...", clicked=self.get_smi_file)
        self.load_path_layout.addWidget(button_load_path)

        layout.addLayout(self.load_path_layout)

        # save_path_layout
        self.save_path_layout = QHBoxLayout()
        self.save_path_label = QLineEdit("Path to save img...")
        self.save_path_layout.addWidget(self.save_path_label)

        button_save_path = QPushButton("...", clicked=self.get_save_path)
        self.save_path_layout.addWidget(button_save_path)

        layout.addLayout(self.save_path_layout)

        # name of file
        self.file_name_layout = QHBoxLayout()
        self.file_name_label = QLabel("Set file name: ")
        self.file_name_layout.addWidget(self.file_name_label)

        self.file_name_elabel = QLineEdit("default")
        self.file_name_layout.addWidget(self.file_name_elabel)

        layout.addLayout(self.file_name_layout)

        # annot format
        combo_annot = QComboBox()
        combo_annot.addItems(["Not Annotation", "Smiles Annotation", "Name Annotation"])
        combo_annot.currentTextChanged.connect(self.text_changed)

        layout.addWidget(combo_annot)

        # num rows set
        self.num_rows_layout = QHBoxLayout()
        self.file_name_label = QLabel("Set rows: ")
        self.num_rows_layout.addWidget(self.file_name_label)

        self.num_row_spine = QSpinBox()
        self.num_row_spine.setMinimum(1)
        self.num_row_spine.setMaximum(31)
        self.num_rows_layout.addWidget(self.num_row_spine)

        layout.addLayout(self.num_rows_layout)

        # create and save img
        button = QPushButton("Create Pictures!", clicked=self.getfile)
        layout.addWidget(button)

        # l1 = QLabel()
        # l1.setPixmap(QPixmap("python.jpg"))

    def get_save_path(self):
        fname = QFileDialog.getExistingDirectory(self, "Open file", "/")
        self.save_path_label.setText(fname[0])

    def get_smi_file(self):
        fname = QFileDialog.getOpenFileName(
            self, "Open file", "/", "Image files (*.jpg *.gif)"
        )
        self.load_path_label.setText(fname[0])

    def text_changed(self, s):
        print(s)

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
