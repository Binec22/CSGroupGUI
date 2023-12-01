import sys
#sys.path.append(r'C:\Users\Binec\Desktop\ProjetCSGroup')
from interface import *

style_sheet = """
    QWidget {
        background-color: #333333;
        color: #ffffff;
    } 
    
    QPushButton {
        background-color: #555555;
        color: #ffffff;
        border: none;
        padding: 10px;
    }
    
    QPushButton:hover {
        background-color: #666666;
    }
"""


if __name__ == "__main__":
    
    app = QtWidgets.QApplication(sys.argv)
    app.setStyleSheet(open(os.path.abspath(os.path.join(os.path.dirname(__file__), 'StyleSheets/MacOS.qss')),'r').read())
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())