import sys
import os
from PyQt6.QtWidgets import QLabel, QWidget, QVBoxLayout, QPushButton, QHBoxLayout
from PyQt6.QtGui import QPixmap, QImage

import numpy as np
from PIL import Image

from ..generate_property import generate_property
from ..processing import mol_to_numpy


class PropertyWidget(QWidget):
    def __init__(self, load_path_label):
        super().__init__()

        self.current_index = 0

        self.load_path_label = load_path_label

        layout = QVBoxLayout()
        self.setLayout(layout)

        self.button_layout = QHBoxLayout()

        self.generate_button = QPushButton(
            "Generate Property!",
            clicked=self.create_property,
        )
        self.button_layout.addWidget(self.generate_button)

        self.left_button = QPushButton("<<", clicked=self.show_previous_image)
        self.left_button.setEnabled(False)
        self.button_layout.addWidget(self.left_button)

        self.right_button = QPushButton(">>", clicked=self.show_next_image)
        self.right_button.setEnabled(False)
        self.button_layout.addWidget(self.right_button)

        layout.addLayout(self.button_layout)

        print(os.listdir())
        file_path = r"src/imgs/default.png"
        image_array = np.array(Image.open(file_path))

        # load_path_layout
        # Convert NumPy array to QImage
        height, width, channel = image_array.shape
        bytes_per_line = 3 * width
        q_image = QImage(
            image_array.data,
            width,
            height,
            bytes_per_line,
            QImage.Format.Format_RGB888,
        )

        # Create QPixmap from QImage
        pixmap = QPixmap.fromImage(q_image)

        # Create QLabel and set the QPixmap
        self.mol_img = QLabel(self)
        self.mol_img.setPixmap(pixmap)

        layout.addWidget(self.mol_img)

        self.property_labels_layout = QVBoxLayout()
        self.name_label = QLabel("Name: ???")
        self.property_labels_layout.addWidget(self.name_label)

        self.weight_label = QLabel("Weight: ???")
        self.property_labels_layout.addWidget(self.weight_label)

        self.logP_label = QLabel("logP: ???")
        self.property_labels_layout.addWidget(self.logP_label)

        self.num_ha_acceptors_label = QLabel("NumHAcceptors: ???")
        self.property_labels_layout.addWidget(self.num_ha_acceptors_label)

        self.num_hd_donors_label = QLabel("NumHDonors: ???")
        self.property_labels_layout.addWidget(self.num_hd_donors_label)

        layout.addLayout(self.property_labels_layout)

    def create_property(self):
        path_to_smi = self.load_path_label.text()
        self.property = generate_property(path_to_smi)
        if len(self.property) == 0:
            return None
        self.smiles_list = list(self.property.keys())

        self.update(self.current_index)

        self.right_button.setEnabled(True)
        self.left_button.setEnabled(True)

    def update(self, index):
        self.update_image(index)
        self.update_property(index)

    def update_image(self, idx):
        image_array = mol_to_numpy(self.smiles_list[idx])

        height, width, channel = image_array.shape
        bytes_per_line = 3 * width
        q_image = QImage(
            image_array.data, width, height, bytes_per_line, QImage.Format.Format_RGB888
        )

        # Create QPixmap from QImage
        pixmap = QPixmap.fromImage(q_image)
        self.mol_img.setPixmap(pixmap)

    def update_property(self, idx):
        self.name_label.setText(
            f'Name: {self.property[self.smiles_list[idx]]["name"]}',
        )
        self.weight_label.setText(
            f'Weight: {self.property[self.smiles_list[idx]]["weight"]}',
        )
        self.logP_label.setText(
            f'logP: {self.property[self.smiles_list[idx]]["logP"]}',
        )
        self.num_ha_acceptors_label.setText(
            f'NumHAcceptors: {self.property[self.smiles_list[idx]]["num_ha_acceptors"]}',
        )
        self.num_hd_donors_label.setText(
            f'NumHDonors: {self.property[self.smiles_list[idx]]["num_hd_donors"]}',
        )

    def show_previous_image(self):
        self.current_index = max(self.current_index - 1, 0)
        self.update_image(self.current_index)
        self.update_property(self.current_index)

    def show_next_image(self):
        self.current_index = min(self.current_index + 1, len(self.smiles_list) - 1)
        self.update_image(self.current_index)
        self.update_property(self.current_index)
