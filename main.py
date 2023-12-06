import sys
import os
from PyQt6.QtWidgets import QApplication

from src.mainApp import MainApp

if __name__ == "__main__":
    # Create the application and run it
    app = QApplication(sys.argv)
    ex = MainApp()
    ex.show()
    sys.exit(app.exec())
