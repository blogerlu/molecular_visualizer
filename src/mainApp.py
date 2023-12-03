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

from generate_img import processing_smi


class App(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecular Visualizer")
        self.resize(600, 200)

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
        self.combo_annot = QComboBox()
        self.combo_annot.addItems(
            ["Not Annotation", "Smiles Annotation", "Name Annotation"]
        )
        self.combo_annot.currentTextChanged.connect(self.combo_annot_text_changed)

        layout.addWidget(self.combo_annot)

        # num rows set
        self.num_rows_layout = QHBoxLayout()
        self.file_name_label = QLabel("Set columns: ")
        self.num_rows_layout.addWidget(self.file_name_label)

        self.num_row_spine = QSpinBox()
        self.num_row_spine.setMinimum(1)
        self.num_row_spine.setMaximum(31)
        self.num_rows_layout.addWidget(self.num_row_spine)

        layout.addLayout(self.num_rows_layout)

        # set img size
        self.img_size_layout = QHBoxLayout()
        self.img_size_width_label = QLabel("Img width: ")
        self.img_size_layout.addWidget(self.img_size_width_label)

        self.img_size_width_elabel = QLineEdit("200")
        self.img_size_layout.addWidget(self.img_size_width_elabel)

        self.img_size_height_label = QLabel("Img height: ")
        self.img_size_layout.addWidget(self.img_size_height_label)

        self.img_size_height_elabel = QLineEdit("200")
        self.img_size_layout.addWidget(self.img_size_height_elabel)

        layout.addLayout(self.img_size_layout)

        # create and save img
        button = QPushButton("Create Pictures!", clicked=self.save_img)
        layout.addWidget(button)

        # l1 = QLabel()
        # l1.setPixmap(QPixmap("python.jpg"))

    def get_save_path(self):
        fname = QFileDialog.getExistingDirectory(self, "Open file", "/")
        self.save_path_label.setText(fname[0])

    def get_smi_file(self):
        fname = QFileDialog.getOpenFileName(self, "Open file", "/", "Image files (*)")
        self.load_path_label.setText(fname[0])

    def combo_annot_text_changed(self, s):
        print(s)

    def save_img(self):
        path_to_smi = self.load_path_label.text()
        print("file", self.file_name_elabel.text())
        path_to_save = os.path.join(
            self.save_path_label.text(), self.file_name_elabel.text() + ".png"
        )
        mols_per_row = self.num_row_spine.value()
        sub_img_size = (
            int(self.img_size_width_elabel.text()),
            int(self.img_size_height_elabel.text()),
        )
        smiles_name = self.combo_annot.currentText()

        print(path_to_save)
        processing_smi(
            path_to_smi,
            path_to_save,
            smiles_name,
            mols_per_row,
            sub_img_size,
        )


if __name__ == "__main__":
    # Create the application and run it
    app = QApplication(sys.argv)
    ex = App()
    ex.show()
    sys.exit(app.exec())
