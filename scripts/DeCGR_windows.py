from PyQt5.QtWidgets import (
    QMainWindow, QLabel, QToolBar, QAction, QWidget, QVBoxLayout, QStackedWidget, QPushButton, QToolButton, QSplitter,
    QHBoxLayout, QTextEdit, QFrame, QApplication
)
from PyQt5.QtGui import QPixmap, QTextOption, QFont, QPainter, QColor, QDesktopServices
from PyQt5.QtCore import QUrl, Qt
import sys
from PyQt5 import QtCore, QtGui, QtWidgets
from decgr.DetermineFragment import FragmentAssembly
from decgr.simuProcess import simulation, get_cut_pos, plot_simu_matrix, read_assembly_result
from decgr.plot_reconstruct_map import global_expect, joint_SV, adjust_matrix, plot_assembly_matrix, read_one_assembly
from decgr.utils import *
import cooler, re
from decgr import correct
from decgr.plot_SV_region import extra_SV_region_inter, extra_SV_region_intra, plot_current_SV
from decgr.runHiSV import HiSV
import pkg_resources
os.environ['PYTHONIOENCODING'] = 'utf-8'

class CustomToolBar(QToolBar):
    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.fillRect(self.rect(), QColor(53, 58, 64))

class Interface(QMainWindow):
    def __init__(self, parent=None):
        super(Interface, self).__init__(parent)
        self.setupUi(self)
        self.stackedWidget = QStackedWidget()
        self.setCentralWidget(self.stackedWidget)
        self.path = os.path.split(os.path.abspath(__file__))[0]
        if not os.path.exists(self.path + '/result'):
            os.mkdir(self.path + '/result')
        if not os.path.exists(self.path + '/tmp'):
            os.mkdir(self.path + '/tmp')
        if not os.path.exists(self.path + '/cnv_norm_data'):
            os.mkdir(self.path + '/cnv_norm_data')

        self.setWindowTitle('')
        self.setWindowFlags(self.windowFlags() & ~QtCore.Qt.WindowMaximizeButtonHint)
        self.desktop = QApplication.desktop()
        self.screenRect = self.desktop.screenGeometry()
        self.screenwidth = self.screenRect.width()
        self.screenheight = self.screenRect.height()

        self.height = int(self.screenheight * 0.8)
        self.width = int(self.screenwidth * 0.8)
        self.setFixedSize(self.width, self.height)
        self.move(self.screenwidth * 0.1, self.screenheight * 0.1)

        self.setPalette(QtGui.QPalette(QtCore.Qt.white))
        self.setAutoFillBackground(True)

        self.setStyleSheet("""
            QPushButton {
                background-color: rgb(53, 58, 64);
                color: white;
                border-radius: 15px;
                border: 2px solid transparent;
            }
            QPushButton:hover {
                border: 2px solid white;
            }
            QPushButton:disabled {
                background-color: rgb(128, 128, 128);
                color: rgb(200, 200, 200);
                border: 2px solid rgb(128, 128, 128);
            }
        """)

        self.setStyleSheet('QAbstractButton {font-size: 20px;} .QLabel {font-size: 20px;} \
                     .QLineEdit {font-size: 14px;} #toLabel {font-size: 48px;} .QComboBox {font-size: 14px;} \
                     .QComboBox:disabled {color: rgb(0,0,0)} .QListWidget {font-size: 14px;}')
        # module 1
        # calling SV and visualization SV from contact map
        self.mainWidget = QWidget()
        self.splitter = QSplitter(QtCore.Qt.Horizontal)
        self.splitter.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.splitter.addWidget(self.inputPage)
        self.splitter.addWidget(self.page)
        self.splitter.setSizes([self.width // 3, self.width * 2 // 3])
        mainLayout = QHBoxLayout(self.mainWidget)
        mainLayout.addWidget(self.splitter)

        # module 2
        # assembling CSV from input HiC data
        self.mainWidget_1 = QWidget()
        self.splitter_1 = QSplitter(QtCore.Qt.Horizontal)
        self.splitter_1.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.splitter_1.addWidget(self.inputPage1)
        self.splitter_1.addWidget(self.page_1)
        self.splitter_1.setSizes([self.width // 3, self.width * 2 // 3])
        mainLayout_1 = QHBoxLayout(self.mainWidget_1)
        mainLayout_1.addWidget(self.splitter_1)

        # module 3
        # Validation CSV by simulating process
        self.mainWidget_2 = QWidget()
        self.splitter_2 = QSplitter(QtCore.Qt.Horizontal)
        self.splitter_2.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.splitter_2.addWidget(self.inputPage2)
        self.splitter_2.addWidget(self.page_2)
        self.splitter_2.setSizes([self.width // 3, self.width * 2 // 3])
        mainLayout_2 = QHBoxLayout(self.mainWidget_2)
        mainLayout_2.addWidget(self.splitter_2)

        # module 4
        # Visual CSV
        self.mainWidget_3 = QWidget()
        self.splitter_3 = QSplitter(QtCore.Qt.Horizontal)
        self.splitter_3.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.splitter_3.addWidget(self.inputPage3)
        self.splitter_3.addWidget(self.page_3)
        self.splitter_3.setSizes([self.width // 3, self.width * 2 // 3])
        mainLayout_3 = QHBoxLayout(self.mainWidget_3)
        mainLayout_3.addWidget(self.splitter_3)

        self.stackedWidget.addWidget(self.mainPage)
        self.stackedWidget.addWidget(self.mainWidget)
        self.stackedWidget.addWidget(self.mainWidget_1)
        self.stackedWidget.addWidget(self.mainWidget_2)
        self.stackedWidget.addWidget(self.mainWidget_3)

        toolbar = CustomToolBar("My main toolbar")
        toolbar.setFixedHeight(50)
        toolbar.setStyleSheet("background-color: rgb(53, 58, 64); color: white; font-size: 20px;")
        self.addToolBar(toolbar)
        self.current_page = 0
        self.tool_button_main = self.create_toolbar_button("DeCGR", 0, toolbar)
        self.tool_button_breakpoint = self.create_toolbar_button("Breakpoint Filtering", 1, toolbar)
        self.tool_button_fragment = self.create_toolbar_button("Fragment Assembly", 2, toolbar)
        self.tool_button_validation = self.create_toolbar_button("Validation CGRs", 3, toolbar)
        self.tool_button_reconstruct = self.create_toolbar_button("Reconstruct Hi-C Map", 4, toolbar)
        self.update_toolbar_buttons()
        # page(detect and visual SV) connect button
        self.runHiSVButton.clicked.connect(self.HiSV_parameter_setting)
        self.hicFolderButton.clicked.connect(lambda: self.inputFilePath('hicfile'))
        self.pageSubmitButton.clicked.connect(self.visualSV_step)
        self.pageResetButton.clicked.connect(self.clear_plotSV)
        self.hicsvLoadButton.clicked.connect(self.set_sv_resolution)
        # page1(assemble CSV) connect button
        self.tumorFolderButton.clicked.connect(lambda: self.inputFilePath('tumor'))
        self.svFileButton.clicked.connect(self.load_chrom_from_file)
        re = QtCore.QRegExp('^[1-9][0-9]*$')
        validator = QtGui.QRegExpValidator(re)
        self.pre_cnvButton.clicked.connect(self.pre_cnv_norm_parameter_setting)

        self.page1SubmitButton.clicked.connect(self.step1)
        self.LoadButton.clicked.connect(self.set_tumor_resolution)
        # page2(validation CSV) connect button
        self.conFolderButton.clicked.connect(lambda: self.inputFilePath('control'))
        self.SimutumorFolderButton.clicked.connect(lambda: self.inputFilePath('simutumor'))
        self.simuLoadButton.clicked.connect(self.set_simusample_resolution)
        self.assemblyFileButton.clicked.connect(self.load_chrom_from_assemfile)
        self.page2SubmitButton.clicked.connect(self.step2)
        # page3
        self.rehicFolderButton.clicked.connect(lambda: self.inputFilePath('reconstruct'))
        self.reLoadButton.clicked.connect(self.set_reconstruct_resolution)
        self.cnvButton.clicked.connect(self.cnv_norm_parameter_setting)
        self.reAssemblyFileButton.clicked.connect(self.load_chrom_from_CGRfile)
        self.page3SubmitButton.clicked.connect(self.step3)
        self.exportResult.clicked.connect(self.saveResult)

    def create_toolbar_button(self, label, page_index, toolbar):
        button_action = QAction(label, self)
        button_action.setStatusTip(f"Go to Page {page_index}")
        tool_button = QToolButton()
        tool_button.setDefaultAction(button_action)
        tool_button.clicked.connect(lambda: self.switch_page(page_index))
        toolbar.addWidget(tool_button)
        return tool_button

    def switch_page(self, page_index):
        self.stackedWidget.setCurrentIndex(page_index)
        self.current_page = page_index
        self.update_toolbar_buttons()

    def update_toolbar_buttons(self):
        buttons = [self.tool_button_main, self.tool_button_breakpoint, 
                   self.tool_button_fragment, self.tool_button_validation, 
                   self.tool_button_reconstruct]
        for i, button in enumerate(buttons):
            if i == self.current_page:
                button.setStyleSheet("color: rgb(245, 195, 66); font-weight: bold;")
            else:
                button.setStyleSheet("color: white; font-weight: normal;")

    def show_error_message(self, message):
        # error
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText(message)
        msg.setWindowTitle("Error")
        msg.exec_()

    def reset_all_inputs(self):
        for widget in self.findChildren(QtWidgets.QLineEdit):
            widget.clear()
        for widget in self.findChildren(QtWidgets.QComboBox):
            widget.setCurrentIndex(0)
        for widget in self.findChildren(QtWidgets.QCheckBox):
            widget.setChecked(False)

    def inputFilePath(self, tag):
        # show the input path
        if tag == 'hicfile':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', "")[0]
            if path == '':
                self.show_error_message("This is a error file path.")
            self.hicFolder.setText(path)
        if tag == 'tumor':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', "")[0]
            if path == '':
                self.show_error_message("This is a error file path.")
            self.tumorFolder.setText(path)
        elif tag == 'simutumor':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', "")[0]
            if path == '':
                self.show_error_message("This is a error file path.")
            self.simutumorFolder.setText(path)
        elif tag == 'sv':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', "")[0]
            if path == '':
                self.show_error_message("This is a error file path.")
            self.svFile.setText(path)
        elif tag == 'control':
            # path = QtWidgets.QFileDialog.getExistingDirectory(self, 'Input Directory Path', self.path)
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', "")[0]
            if path == '':
                self.show_error_message("This is a error file path.")
            self.conFolder.setText(path)
        elif tag == 'assembly':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', "")[0]
            if path == '':
                self.show_error_message("This is a error file path.")
            self.assemblyFile.setText(path)
        elif tag == 're':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', "")[0]
            if path == '':
                self.show_error_message("This is a error file path.")
            self.reAssemblyFile.setText(path)
        elif tag == 'reconstruct':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', "")[0]
            if path == '':
                self.show_error_message("This is a error file path.")
            self.rehicFolder.setText(path)


    def set_simusample_resolution(self):
        tumor_path = self.simutumorFolder.text()
        control_path = self.conFolder.text()
        _, tumor_file_extension = os.path.splitext(tumor_path)
        _, control_file_extension = os.path.splitext(control_path)
        if tumor_file_extension.lower() in ['.cool', '.mcool'] and control_file_extension.lower() in ['.cool', '.mcool'] :
            if 'mcool' in tumor_path:
                resolutions = cooler.fileops.list_coolers(tumor_path)
                tumor_resolutions = [res.split('/')[-1] for res in resolutions if '/resolutions/' in res]
            else:
                clr = cooler.Cooler(tumor_path)
                cur_resolutions = clr.info['bin-size']
                tumor_resolutions = [str(cur_resolutions)]
            if 'mcool' in control_path:
                resolutions = cooler.fileops.list_coolers(control_path)
                control_resolutions = [res.split('/')[-1] for res in resolutions if '/resolutions/' in res]
            else:
                clr = cooler.Cooler(control_path)
                cur_resolutions = clr.info['bin-size']
                control_resolutions = [str(cur_resolutions)]
            common_resolutions = list(set(control_resolutions) & set(tumor_resolutions))
            common_resolutions = sorted(common_resolutions, key=int)
            if len(common_resolutions) == 0:
                self.show_error_message("The two cool files do not have common resolutions.")
            self.simuresCombox.clear()
            self.simuresCombox.addItems(common_resolutions)
            self.assemblyFileButton.setEnabled(True)
        else:
            self.show_error_message("This is a wrong file, it should be a cool format file.")

    def set_tumor_resolution(self):
        tumor_path = self.tumorFolder.text()
        _, tumor_file_extension = os.path.splitext(tumor_path)
        if tumor_file_extension.lower() in ['.cool', '.mcool']:
            if 'mcool' in tumor_path:
                resolutions = cooler.fileops.list_coolers(tumor_path)
                tumor_resolutions = [res.split('/')[-1] for res in resolutions if '/resolutions/' in res]
            else:
                clr = cooler.Cooler(tumor_path)
                cur_resolutions = clr.info['bin-size']
                tumor_resolutions = [str(cur_resolutions)]
            
            tumor_resolutions = sorted(tumor_resolutions, key=int)
            if len(tumor_resolutions) == 0:
                self.show_error_message("The two cool files do not have common resolutions.")
            self.resCombox.clear()
            self.resCombox.addItems(tumor_resolutions)
            self.svFileButton.setEnabled(True)
        else:
            self.show_error_message("This is a wrong file, it should be a cool format file.")

    def set_sv_resolution(self):
        tumor_path = self.hicFolder.text()
        _, tumor_file_extension = os.path.splitext(tumor_path)
        if tumor_file_extension.lower() in ['.cool', '.mcool']:
            if 'mcool' in tumor_path:
                resolutions = cooler.fileops.list_coolers(tumor_path)
                tumor_resolutions = [res.split('/')[-1] for res in resolutions if '/resolutions/' in res]
            else:
                clr = cooler.Cooler(tumor_path)
                cur_resolutions = clr.info['bin-size']
                tumor_resolutions = [str(cur_resolutions)]
            tumor_resolutions = sorted(tumor_resolutions, key=int)
            if len(tumor_resolutions) == 0:
                self.show_error_message("The two cool files do not have common resolutions.")
            self.svresCombox.clear()
            self.svresCombox.addItems(tumor_resolutions)
            self.runHiSVButton.setEnabled(True)
            self.pageSubmitButton.setEnabled(True)
        else:
            self.show_error_message("This is a wrong file, it should be a cool format file.")

    def set_reconstruct_resolution(self):
        tumor_path = self.rehicFolder.text()
        _, tumor_file_extension = os.path.splitext(tumor_path)
        if tumor_file_extension.lower() in ['.cool', '.mcool']:
            if 'mcool' in tumor_path:
                resolutions = cooler.fileops.list_coolers(tumor_path)
                tumor_resolutions = [res.split('/')[-1] for res in resolutions if '/resolutions/' in res]
            else:
                clr = cooler.Cooler(tumor_path)
                cur_resolutions = clr.info['bin-size']
                tumor_resolutions = [str(cur_resolutions)]
            tumor_resolutions = sorted(tumor_resolutions, key=int)
            if len(tumor_resolutions) == 0:
                self.show_error_message("The two cool files do not have common resolutions.")
            self.reresCombox.clear()
            self.reresCombox.addItems(tumor_resolutions)
            self.reAssemblyFileButton.setEnabled(True)
        else:
            self.show_error_message("This is a wrong file, it should be a cool format file.")
        
    def HiSV_parameter_setting(self):
        self.dialog = QtWidgets.QDialog(self)
        self.dialog.setWindowTitle("Parameter Input")
        self.param1Label = QLabel("Chrom List:")
        self.param1Input = QtWidgets.QLineEdit(self.dialog)
        self.param1Input.setPlaceholderText("e.g. chr1,chr2,chr3")
        self.param2Label = QLabel("window:")
        self.param2Input = QtWidgets.QLineEdit(self.dialog)
        self.param2Input.setText("10")
        self.param3Label = QLabel("weight:")
        self.param3Input = QtWidgets.QLineEdit(self.dialog)
        self.param3Input.setText("0.2")
        self.param4Label = QLabel("cutoff:")
        self.param4Input = QtWidgets.QLineEdit(self.dialog)
        self.param4Input.setText("0.6")
        self.okButton = QPushButton("Run", self.dialog)
        self.okButton.clicked.connect(self.torunHiSV)
        layout = QVBoxLayout()
        param1Layout = QHBoxLayout()
        param1Layout.addWidget(self.param1Label)
        param1Layout.addWidget(self.param1Input)
        layout.addLayout(param1Layout)
        param2Layout = QHBoxLayout()
        param2Layout.addWidget(self.param2Label)
        param2Layout.addWidget(self.param2Input)
        layout.addLayout(param2Layout)
        param3Layout = QHBoxLayout()
        param3Layout.addWidget(self.param3Label)
        param3Layout.addWidget(self.param3Input)
        layout.addLayout(param3Layout)
        param4Layout = QHBoxLayout()
        param4Layout.addWidget(self.param4Label)
        param4Layout.addWidget(self.param4Input)
        layout.addLayout(param4Layout)
        layout.addWidget(self.okButton)
        self.dialog.setLayout(layout) 
        self.dialog.exec_()

    def cnv_norm_parameter_setting(self):
        ## is 'sweight'
        hic_path = self.rehicFolder.text()
        res = int(self.reresCombox.currentText())
        if 'mcool' in hic_path:
            hic_file = hic_path + '::/resolutions/' + str(res)
        else:
            hic_file = hic_path
        self.cnv_norm_hic_file = hic_file
        
        self.dialog1 = QtWidgets.QDialog(self)
        self.dialog1.setWindowTitle("Parameter Input")
        # Ref genome
        self.norm_param1Label = QLabel("Ref Genome:")
        self.norm_param1Input = QtWidgets.QComboBox(self.dialog1)
        self.norm_param1Input.addItems(['hg19', 'hg38', 'mm9', 'mm10'])
        # Enzyme
        self.norm_param2Label = QLabel("Enzyme:")
        self.norm_param2Input = QtWidgets.QComboBox(self.dialog1)
        self.norm_param2Input.addItems(['HindIII', 'MboI', 'DpnII', 'BglII', 'Arima', 'uniform'])
        # nproc
        self.norm_param3Label = QLabel("nproc:")
        self.norm_param3Input = QtWidgets.QLineEdit(self.dialog1)
        self.norm_param3Input.setText("10")
        self.norm_okButton = QPushButton("Run", self.dialog1)
        self.norm_okButton.clicked.connect(self.cnv_correct)
        # setting
        layout = QVBoxLayout()
        norm_param1Layout = QHBoxLayout()
        norm_param1Layout.addWidget(self.norm_param1Label)
        norm_param1Layout.addWidget(self.norm_param1Input)
        layout.addLayout(norm_param1Layout)
        norm_param2Layout = QHBoxLayout()
        norm_param2Layout.addWidget(self.norm_param2Label)
        norm_param2Layout.addWidget(self.norm_param2Input)
        layout.addLayout(norm_param2Layout)
        norm_param3Layout = QHBoxLayout()
        norm_param3Layout.addWidget(self.norm_param3Label)
        norm_param3Layout.addWidget(self.norm_param3Input)
        layout.addLayout(norm_param3Layout)        
        layout.addWidget(self.norm_okButton)
        self.dialog1.setLayout(layout) 
        self.dialog1.exec_()

    def pre_cnv_norm_parameter_setting(self):
        hic_path = self.tumorFolder.text()
        res = int(self.resCombox.currentText())
        if 'mcool' in hic_path:
            hic_file = hic_path + '::/resolutions/' + str(res)
        else:
            hic_file = hic_path
        self.pre_cnv_norm_hic_file =  hic_file
        self.dialog2 = QtWidgets.QDialog(self)
        self.dialog2.setWindowTitle("Parameter Input")
        # Ref genome
        self.pre_param1Label = QLabel("Ref Genome:")
        self.pre_param1Input = QtWidgets.QComboBox(self.dialog2)
        self.pre_param1Input.addItems(['hg19', 'hg38', 'mm9', 'mm10'])
        # Enzyme
        self.pre_param2Label = QLabel("Enzyme:")
        self.pre_param2Input = QtWidgets.QComboBox(self.dialog2)
        self.pre_param2Input.addItems(['HindIII', 'MboI', 'DpnII', 'BglII', 'Arima', 'uniform'])
        # nproc
        self.pre_param3Label = QLabel("nproc:")
        self.pre_param3Input = QtWidgets.QLineEdit(self.dialog2)
        self.pre_param3Input.setText("10")
        self.pre_okButton = QPushButton("Run", self.dialog2)
        self.pre_okButton.clicked.connect(self.pre_cnv_correct)
        # setting
        layout = QVBoxLayout()
        pre_param1Layout = QHBoxLayout()
        pre_param1Layout.addWidget(self.pre_param1Label)
        pre_param1Layout.addWidget(self.pre_param1Input)
        layout.addLayout(pre_param1Layout)
        pre_param2Layout = QHBoxLayout()
        pre_param2Layout.addWidget(self.pre_param2Label)
        pre_param2Layout.addWidget(self.pre_param2Input)
        layout.addLayout(pre_param2Layout)
        pre_param3Layout = QHBoxLayout()
        pre_param3Layout.addWidget(self.pre_param3Label)
        pre_param3Layout.addWidget(self.pre_param3Input)
        layout.addLayout(pre_param3Layout)        
        layout.addWidget(self.pre_okButton)
        self.dialog2.setLayout(layout) 
        self.dialog2.exec_()

    def torunHiSV(self):
        self.okButton.setEnabled(False)
        chromList = self.param1Input.text().split(',')
        window = int(self.param2Input.text())
        weight = float(self.param3Input.text())
        cutoff = float(self.param4Input.text())
        hic_path = self.hicFolder.text()
        res = int(self.svresCombox.currentText())
        if 'mcool' in hic_path:
            hic_file = hic_path + '::/resolutions/' + str(res)
        else:
            hic_file = hic_path
        clr = cooler.Cooler(hic_file)
        # get chrom_list
        self.chrom_info = clr.chroms()[:]
        chrom_names = self.chrom_info['name'].tolist()
        is_subset = all(item in chrom_names for item in chromList)
        if is_subset == False:
            self.show_error_message("Invalid chromosome.")
            return
        if weight <= 0 or cutoff <= 0 or cutoff >= 1:
            self.show_error_message("Invalid parameter for HiSV.")
            return
        result = HiSV(clr, chromList, window, weight, cutoff)
        self.okButton.setEnabled(True)
        self.dialog.accept()
        output_text = ""
        output_text += ("Position1" + '\t' + "Position2" + '\n')
        for k in range(len(result['pos1_start'])):
            cur_pos1 = str(result['chrom1'][k]) + ":" + str(result['pos1_start'][k]) + "-" + str(result['pos1_end'][k])
            cur_pos2 = str(result['chrom2'][k]) + ":" + str(result['pos2_start'][k]) + "-" + str(result['pos2_end'][k])
            output_text += (cur_pos1 + '\t' + cur_pos2 + '\n')
        self.HiSV_result.setText(output_text)


    def region2position(self, region):
        if ":" not in region or "-" not in region:
            self.show_error_message("Invalid genome region!")
            chrom = start = end = None   
        chrom, range_part = region.split(":")
        start, end = range_part.split("-")
        chrom_names = self.chrom_info['name'].tolist()
        if chrom not in chrom_names:
            self.show_error_message("Invalid chromosome!")
            chrom = start = end = None      
        else:
            chrom_length = self.chrom_info.loc[self.chrom_info['name'] == chrom, 'length'].values[0]
            start = int(start)
            end = int(end)
            if start < 0 or start > chrom_length or end < 0 or end > chrom_length or end < start:
                self.show_error_message("Invalid genome region!")
                chrom = start = end = None  
        return chrom, start, end
    
    def visualSV_step(self):
        # chr4:160350000-160500000; chr4:163600000-164100000
        self.pageSubmitButton.setEnabled(False)
        hic_path = self.hicFolder.text()
        res = int(self.svresCombox.currentText())
        if 'mcool' in hic_path:
            hic_file = hic_path + '::/resolutions/' + str(res)
        else:
            hic_file = hic_path
        clr = cooler.Cooler(hic_file)
        resolution = clr.info['bin-size']
        # get chrom_list
        self.chrom_info = clr.chroms()[:]
        pos1 = self.position1Input.text()
        if pos1 == "":
            self.show_error_message("Error genome region!")
            self.pageSubmitButton.setEnabled(True)
            return
        chrom1, start1, end1 = self.region2position(pos1)
        if chrom1 is None:
            self.show_error_message("Error genome region!")
            self.pageSubmitButton.setEnabled(True)
            return
        pos2 = self.position2Input.text()
        if pos2 == "":
            self.show_error_message("Error genome region!")
            self.pageSubmitButton.setEnabled(True)
            return
        chrom2, start2, end2 = self.region2position(pos2)
        if chrom2 is None:
            self.show_error_message("Error genome region!")
            self.pageSubmitButton.setEnabled(True)
            return
        add_region = 5000000
        chrom_list = []
        start_list = []
        end_list = []
        SV_list = []
        if chrom1 == chrom2:
            # intra SV
            if start1 < start2:
                region_start = start1-add_region
            else:
                region_start = start2-add_region
            if region_start < 0:
                region_start = 0
            if end1 < end2:
                region_end = end2+add_region
            else:
                region_end = end1+add_region
            chrom_length = self.chrom_info.loc[self.chrom_info['name'] == chrom1, 'length'].values[0]
            if region_end > chrom_length:
                region_end = chrom_length
            # SV block
            SV_list.append((start1-region_start)//resolution)
            SV_list.append((end1-region_start)//resolution)
            SV_list.append((start2-region_start)//resolution)
            SV_list.append((end2-region_start)//resolution)
            chrom_list.append(chrom1)
            start_list.append(region_start)
            end_list.append(region_end)
            region = chrom1 + ':' + str(region_start) + '-'  + str(region_end)
            matrix, interval = extra_SV_region_intra(clr, region)
            canvas = plot_current_SV(matrix, interval, chrom_list, start_list, end_list, SV_list, resolution, self.Interaction_show)
        else:
            # inter SV
            chrom_length = self.chrom_info.loc[self.chrom_info['name'] == chrom1, 'length'].values[0]
            region_start = max(0, start1-add_region)
            region_end = min(end1+add_region, chrom_length)
            region1 = chrom1 + ':' + str(region_start) + '-'  + str(region_end)
            chrom_list.append(chrom1)
            start_list.append(region_start)
            end_list.append(region_end)
            SV_list.append((start1-region_start)//resolution)
            SV_list.append((end1-region_start)//resolution)
            chrom_length = self.chrom_info.loc[self.chrom_info['name'] == chrom2, 'length'].values[0]
            region_start = max(0, start2-add_region)
            region_end = min(end2+add_region, chrom_length)
            region2 = chrom2 + ':' + str(region_start) + '-'  + str(region_end)
            chrom_list.append(chrom2)
            start_list.append(region_start)
            end_list.append(region_end)
            matrix, interval = extra_SV_region_inter(clr, region1, region2)
            # SV block
            SV_list.append((start2-region_start)//resolution+interval[0])
            SV_list.append((end2-region_start)//resolution+interval[0])
            canvas = plot_current_SV(matrix, interval, chrom_list, start_list, end_list, SV_list, resolution, self.Interaction_show)
        self.hicmap_SV = canvas
        self.plotSV()
        # set button enable
        self.hicFolderButton.setEnabled(False)
        self.hicsvLoadButton.setEnabled(False)
        self.runHiSVButton.setEnabled(False)
        self.pageSubmitButton.setEnabled(False)

    def plotSV(self):
        parent_width = self.page.width()
        parent_height = self.page.height()
        self.verticalLayout.addWidget(self.hicmap_SV)
        self.hicmap_SV.setFixedSize(parent_width, parent_height)
        self.pageResetButton.setEnabled(True)

    def clear_plotSV(self):
        self.position1Input.clear()
        self.position2Input.clear()
        self.Interaction_show.clear()
        while self.verticalLayout.count():
            item = self.verticalLayout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()
        self.verticalLayout.update()
        self.pageSubmitButton.setEnabled(True)

    def load_chrom_from_file(self):
        file_name, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select CGR file", "", "All Files (*)" )
        if file_name:
            chr_list = []
            self.svFile.setText(file_name)
            with open(file_name, 'r') as f:
                for line in f.readlines():
                    linelist = line.strip().split('\t')
                    if len(linelist) == 4:
                        chr_list.append(linelist[0])
                        chr_list.append(linelist[2])
                    else:    
                        self.show_error_message("This is a wrong rearrangement breakpoint file.")
                        return
            chr_list = list(set(chr_list))
            chr_list = sorted(chr_list)
            self.chrom_list = chr_list
            chrom_text = ""
            for cur_chrom in chr_list:
                chrom_text += (cur_chrom)
                chrom_text += '\n'
            self.chromlist_show.setText(chrom_text)
            self.pre_cnvButton.setEnabled(True)
            self.tumorFolderButton.setEnabled(False)
            self.LoadButton.setEnabled(False)
            self.svFileButton.setEnabled(False)
        else:
            self.show_error_message(f"Error reading the rearrangement breakpoint file.")
            return
        
    def load_chrom_from_assemfile(self):
        # load chromosomes from input assemble file
        output_text = ""
        file_name, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select assemble file", "", "All Files (*)" )
        if file_name:
            chr_list = []
            self.assemblyFile.setText(file_name)
            with open(file_name, 'r') as f:
                f.readline()
                for line in f.readlines():
                    output_text += line
                    linelist = line.strip().split('\t')
                    if len(linelist) == 6:
                        chr_list.append(linelist[1])
                    else:    
                        self.show_error_message("This is a wrong CGR file.")
                        return
            chr_list = list(set(chr_list))
            self.simuchrom_list = chr_list
            self.simuchromlist_show.setText(output_text)
            self.page2SubmitButton.setEnabled(True)
            self.simuLoadButton.setEnabled(False)
            self.conFolderButton.setEnabled(False)
            self.assemblyFileButton.setEnabled(False)
            self.SimutumorFolderButton.setEnabled(False)
        else:
            self.show_error_message(f"Error reading the CGR file.")
            return

    def load_chrom_from_CGRfile(self):
        # load chromosomes from input CGR file
        file_name, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select CGR file", "", "All Files (*)" )
        if file_name:
            chr_list = []
            self.reAssemblyFile.setText(file_name)
            output_text = "" 
            with open(file_name, 'r') as f:
                head = f.readline()
                output_text += (head) 
                for line in f.readlines():
                    linelist = line.strip().split('\t')
                    if len(linelist) == 5:
                        output_text += (line)
                        chr_list.append(linelist[0])
                    else:    
                        self.show_error_message("This is a wrong CGR file.")
                        return
            chr_list = list(set(chr_list))
            self.rechrom_list = chr_list
            self.fragment_show.setText(output_text) 
            self.rehicFolderButton.setEnabled(False)
            self.reLoadButton.setEnabled(False)
            self.reAssemblyFileButton.setEnabled(False)
            self.cnvButton.setEnabled(True)
        else:
            self.show_error_message(f"Error reading the CGR file.")
            return
        
    def chr_sort_key(self, chr_name):
        match = re.match(r'chr(\d+|X|Y)', chr_name)
        if match:
            part = match.group(1)
            # Map X and Y to higher numbers to sort after numeric chromosomes
            return int(part) if part.isdigit() else (100 if part == 'X' else 101)
        return float('inf') 

    def step1(self):
        hic_file = self.pre_cnv_norm_hic_file
        SV_file = self.svFile.text()
        res = int(self.resCombox.currentText())
        clr = cooler.Cooler(hic_file)
        chrom_list = self.chrom_list
        chrom_list = sorted(chrom_list, key=self.chr_sort_key)
        dict_count = {}
        dict_norm = {}
        for i in range(len(chrom_list)):
            for j in range(i, len(chrom_list)):
                name = str(chrom_list[i]) + '_' + str(chrom_list[j])
                chr_pair = [chrom_list[i], chrom_list[j]]
                cur_matrix = load_norm_matrix(clr, chr_pair)
                cur_matrix *= 1000
                dict_count[chrom_list[i]] = np.shape(cur_matrix)[0]
                dict_norm[name] = cur_matrix
        canvas_matrix, canvas_assembly, CGR_info = FragmentAssembly(SV_file,
                                            res, dict_norm, chrom_list, dict_count, self.path)
        result_file = self.path + "/result/assembly_result.txt"

        output_text = ""
        output_text += ("num" + '\t' + 'chrom' + '\t' + "start" + '\t' + "end" + '\t' + "node" + '\t' + "orient" + '\n')
        with open(result_file, 'w') as out:
            out.write("num" + '\t' + 'chrom' + '\t' + "start" + '\t' + "end" + '\t' + "node" + '\t' + "orient" + '\n')
            for CGR_num, fragment_chrom, fragment_start, fragment_end, fragment_orient, fragment_node in zip(CGR_info['num'], CGR_info['chrom'], CGR_info['start'], CGR_info['end'], CGR_info['orient'], CGR_info['node']):
                out.write(
                        str(int(CGR_num) + 1) + '\t' + fragment_chrom+ '\t' + str(int(fragment_start) * res) + '\t' +
                        str(int(fragment_end) * res) + '\t' + fragment_node + '\t' + fragment_orient + '\n')
                output_text += (
                        str(int(CGR_num) + 1) + '\t' + fragment_chrom+ '\t' + str(int(fragment_start) * res)  + '\t' +
                        str(int(fragment_end) * res) + '\t' + fragment_node + '\t' + fragment_orient + '\n')
        self.CGR_result.setText(output_text)
        self.hicmap1 = canvas_matrix
        # self.link = canvas_link
        self.assembly = canvas_assembly
        self.plot1()
        self.pre_cnvButton.setEnabled(False)
        self.LoadButton.setEnabled(False)
        self.tumorFolderButton.setEnabled(False)
        self.svFileButton.setEnabled(False)
        self.page1SubmitButton.setEnabled(False)

    def plot1(self):
        parent_width = self.page_1.width()
        parent_height = self.page_1.height()
        assembly_width = int(parent_width * 0.2)
        assembly_height = int(parent_height * 0.2)
        self.assembly.setGeometry(0, 0, assembly_width, assembly_height)
        self.verticalLayout_1.addWidget(self.assembly)
        hicmap1_height = int(parent_height * 0.8)
        self.verticalLayout_1.addWidget(self.hicmap1)
        self.hicmap1.setFixedSize(parent_width, hicmap1_height)


    def step2(self):
        result_file = self.assemblyFile.text()
        res = int(self.simuresCombox.currentText())
        tumor_path = self.simutumorFolder.text()
        control_path = self.conFolder.text()
        # load cool file
        if 'mcool' in tumor_path:
            tumor_hic_file = tumor_path + '::/resolutions/' + str(res)
        else:
            tumor_hic_file = tumor_path
        tumor_clr = cooler.Cooler(tumor_hic_file)
        if 'mcool' in control_path:
            control_hic_file = control_path + '::/resolutions/' + str(res)
        else:
            control_hic_file = control_path
        control_clr = cooler.Cooler(control_hic_file)
        simuchrom_list = sorted(self.simuchrom_list, key=lambda x: int(x[3:]))
        dict_matrix = {}
        dict_count = {}
        control_dict = {}
        control_len = []
        for i in range(len(simuchrom_list)):
            for j in range(i, len(simuchrom_list)):
                name = str(simuchrom_list[i]) + '_' + str(simuchrom_list[j])
                chr_pair = [simuchrom_list[i], simuchrom_list[j]]
                cur_control_matrix = load_matrix(control_clr, chr_pair)
                control_dict[name] = cur_control_matrix
                cur_tumor_matrix = load_matrix(tumor_clr, chr_pair)
                dict_matrix[name] = cur_tumor_matrix
                dict_count[simuchrom_list[i]] = np.shape(cur_tumor_matrix)[0]
                if i == j:
                    control_len.append(cur_control_matrix.shape[0])      
        dis2cnt = runFitHiC(self.path, control_dict, simuchrom_list, control_len, res)
        assembly_dict = read_assembly_result(result_file, res)
        simu_dict = simulation(assembly_dict, dis2cnt, control_dict, simuchrom_list,res)
        for key in simu_dict:
            SV_matrix = simu_dict[key]
            filename = self.path + "/result/simu_" + key + '_matrix.txt'
            np.savetxt(filename, SV_matrix, fmt='%.2f', delimiter='\t')
        dict_start, dict_end = get_cut_pos(assembly_dict, res, dict_count)
        canvas = plot_simu_matrix(dict_matrix, simu_dict, dict_start, dict_end, simuchrom_list, assembly_dict,
                                  res, self.simuInteraction_show)
        self.simuMap = canvas
        self.plot2()
        self.SimutumorFolderButton.setEnabled(False)
        self.simuLoadButton.setEnabled(False)
        self.conFolderButton.setEnabled(False) 
        self.assemblyFileButton.setEnabled(False)
        self.page2SubmitButton.setEnabled(False)   
        self.page2visualResetButton.setEnabled(True)

    def plot2(self):
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.addWidget(self.simuMap)

    def cnv_correct(self):
        # run cnv normalization
        self.cnvButton.setEnabled(False)
        res = int(self.reresCombox.currentText())
        hic_file = self.cnv_norm_hic_file
        genome = self.norm_param1Input.currentText()
        enzyme = self.norm_param2Input.currentText()
        nproc = int(self.norm_param3Input.text())
        outfile = correct.run(hic_file, genome, enzyme, res, nproc, self.rechrom_list, self.path)
        self.dialog1.accept()
        self.cnv_norm_hic_file = outfile
        self.page3SubmitButton.setEnabled(True)

    def pre_cnv_correct(self):
        # run cnv normalization
        self.pre_cnvButton.setEnabled(False)
        res = int(self.resCombox.currentText())
        hic_file = self.pre_cnv_norm_hic_file
        genome = self.pre_param1Input.currentText()
        enzyme = self.pre_param2Input.currentText()
        nproc = int(self.pre_param3Input.text())
        outfile = correct.run(hic_file, genome, enzyme, res, nproc, self.chrom_list, self.path)
        self.pre_cnv_norm_hic_file = outfile
        self.dialog2.accept()
        self.page1SubmitButton.setEnabled(True)

    def step3(self):
        result_file = self.reAssemblyFile.text()
        res = int(self.reresCombox.currentText())
        hic_file = self.cnv_norm_hic_file
        path, orient, chrom, start, end = read_one_assembly(result_file, res)
        clr = cooler.Cooler(hic_file) 
        dict_norm = {}
        chrom_list = self.rechrom_list
        # load matrix
        for i in range(len(chrom_list)):
            for j in range(i, len(chrom_list)):
                name = str(chrom_list[i]) + '_' + str(chrom_list[j])
                chr_pair = [chrom_list[i], chrom_list[j]]
                M = load_norm_matrix(clr, chr_pair)
                dict_norm[name] = M

        glob_exp = global_expect(dict_norm, chrom_list)
        SV_matrix = joint_SV(path, orient, chrom, start, end, dict_norm)
        matrix, interval = adjust_matrix(SV_matrix, glob_exp, path, start, end)
        canvas = plot_assembly_matrix(SV_matrix, interval, path, orient)
        self.reconstructMap = canvas
        self.plot3()
        self.SV_matrix = SV_matrix
        self.exportResult.setEnabled(True)
        self.page3SubmitButton.setEnabled(False)
        self.reAssemblyFileButton.setEnabled(False)
        self.rehicFolderButton.setEnabled(False)
        self.reLoadButton.setEnabled(False)
        self.cnvButton.setEnabled(False)

    def plot3(self):
        self.horizontalLayout_3.addWidget(self.reconstructMap)

    def saveResult(self):
        filepath, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', '', 'Text Files (*.txt);;All Files (*)')
        if filepath:
            np.savetxt(filepath, self.SV_matrix, delimiter='\t', fmt='%.2f')

    def resetPage(self):
        self.hicFolder.clear()
        self.HiSV_result.clear()
        self.position1Input.clear()
        self.position2Input.clear()
        self.Interaction_show.clear()
        self.svresCombox.setCurrentIndex(-1)
        while self.verticalLayout.count():
            item = self.verticalLayout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()
        self.verticalLayout.update()
        self.runHiSVButton.setEnabled(False)
        self.pageSubmitButton.setEnabled(False)
        self.pageResetButton.setEnabled(False)
        self.hicFolderButton.setEnabled(True)
        self.hicsvLoadButton.setEnabled(True)

    def resetPage1(self):
        self.tumorFolder.clear()
        self.svFile.clear()
        self.chromlist_show.clear()
        self.CGR_result.clear()
        while self.verticalLayout_1.count():
            item = self.verticalLayout_1.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()
        self.verticalLayout_1.update()
        self.resCombox.setCurrentIndex(-1)
        self.tumorFolderButton.setEnabled(True)
        self.LoadButton.setEnabled(True)
        self.svFileButton.setEnabled(False)
        self.pre_cnvButton.setEnabled(False)
        self.page1SubmitButton.setEnabled(False)

    def visualPage2_reset(self):
        while self.verticalLayout_2.count():
            item = self.verticalLayout_2.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()
        self.verticalLayout_2.update()
        self.assemblyFileButton.setEnabled(True)
        self.page2visualResetButton.setEnabled(False)

    def resetPage2(self):
        self.simutumorFolder.clear()
        self.conFolder.clear()
        self.assemblyFile.clear()
        self.simuchromlist_show.clear()
        self.simuInteraction_show.clear()
        self.simuresCombox.setCurrentIndex(-1)
        self.conFolderButton.setEnabled(True)
        self.SimutumorFolderButton.setEnabled(True)
        self.simuLoadButton.setEnabled(True)
        self.assemblyFileButton.setEnabled(False)
        self.page2SubmitButton.setEnabled(False)
        self.page2visualResetButton.setEnabled(False)

    def resetPage3(self):
        self.rehicFolder.clear()
        self.reAssemblyFile.clear() 
        self.fragment_show.clear()
        self.rehicFolderButton.setEnabled(True)
        self.reLoadButton.setEnabled(True)
        self.reresCombox.setCurrentIndex(-1)
        while self.horizontalLayout_3.count():
            item = self.horizontalLayout_3.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()
        self.horizontalLayout_3.update()
        self.reAssemblyFileButton.setEnabled(False)
        self.cnvButton.setEnabled(False)
        self.page2SubmitButton.setEnabled(False)
        self.exportResult.setEnabled(False)

    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(1395, 940)
        # mainPage: Introduction to DeCGR
        self.mainPage = QtWidgets.QWidget()
        self.mainPage.setObjectName("mainPage")
        left_layout = QVBoxLayout()
        button1 = QPushButton("User manual")
        button1.setStyleSheet(
            "QPushButton { background-color: rgb(245, 195, 66); color: white; border-radius: 15px; border: 2px solid transparent; } QPushButton:hover { border: 2px solid white; }")
        button1.clicked.connect(self.open_user_manual)
        button2 = QPushButton("Download test data")
        button2.setStyleSheet(
            "QPushButton { background-color: rgb(245, 195, 66); color: white; border-radius: 15px; border: 2px solid transparent; } QPushButton:hover { border: 2px solid white; }")
        button2.clicked.connect(self.download_test_data)
        button1.setFixedSize(200, 40)
        button2.setFixedSize(200, 40)
        buttons_layout = QHBoxLayout()
        buttons_layout.addWidget(button1)
        buttons_layout.addWidget(button2)
        left_layout.setContentsMargins(0, 20, 0, 0)
        text_edit = QTextEdit()
        text_edit.setHtml(
            "<b>DeCGR</b> is an interaction toolkit with a set of modules specifically designed for deciphering "
            "CGRs in chromatin contact maps."
            "<ul><li>The <b>Breakpoint Filtering</b> module identifies and filters rearrangement breakpoints based on the input Hi-C sample. </li><br/>"
            "<li>The <b>Fragment Assembly</b> module automatically reconstructs candidate CGRs based on the input breakpoints and the Hi-C sample. </li><br/>"
            "<li>The <b>Validation CGRs</b> module reproduces the Hi-C map containing CGRs through a simulation process to validate CGRs.</li><br/>"
            "<li>The <b>Reconstruct Hi-C Map </b> module provides reconstructed Hi-C contact maps</li></ul>")
        text_edit.setFrameStyle(QFrame.NoFrame)
        font = QFont()
        font.setPointSize(18)
        text_edit.setFont(font)
        text_edit.setFixedWidth(450)
        text_edit.setMaximumHeight(350)
        text_edit.setAlignment(QtCore.Qt.AlignJustify)
        text_edit.setWordWrapMode(QTextOption.WrapAtWordBoundaryOrAnywhere)
        left_layout.addLayout(buttons_layout)
        left_layout.addWidget(text_edit)
        right_layout = QVBoxLayout()
        image_label = QLabel()
        image_path = pkg_resources.resource_filename('decgr', 'mainplot.jpg')
        print(image_path)
        pixmap = QPixmap(image_path)
        image_label.setPixmap(pixmap)
        image_label.setAlignment(QtCore.Qt.AlignCenter)
        right_layout.addWidget(image_label)
        main_layout = QHBoxLayout(self.mainPage)
        main_layout.addLayout(left_layout)
        main_layout.addLayout(right_layout)
        self.mainPage.setLayout(main_layout)
        # inputPage: calling breakpoint and visual SV
        self.inputPage = QtWidgets.QWidget()
        self.inputPage.setObjectName("inputPage")
        self.gridLayout = QtWidgets.QGridLayout(self.inputPage)
        self.gridLayout.setObjectName("gridLayout")
        ## input hicfile 
        self.hicFolder = QtWidgets.QLineEdit(self.inputPage)
        self.hicFolder.setEnabled(True)
        self.hicFolder.setReadOnly(True)
        self.hicFolder.setObjectName("hicFolder")
        self.hicFolder.setPlaceholderText("Select tumor sample (cool format).")
        self.gridLayout.addWidget(self.hicFolder, 0, 0, 1, 2)
        self.hicFolderButton = QtWidgets.QPushButton(self.inputPage)
        self.hicFolderButton.setObjectName("hicFolderButton")
        self.hicFolderButton.setFixedSize(150, 30)
        self.gridLayout.addWidget(self.hicFolderButton, 0, 2, 1, 1)
        self.hicsvLoadButton = QtWidgets.QPushButton(self.inputPage)
        self.hicsvLoadButton.setObjectName("LoadButton")
        self.hicsvLoadButton.setFixedSize(150, 30)
        self.gridLayout.addWidget(self.hicsvLoadButton, 1, 2, 1, 1)
        ## resolution setting
        self.resolutionLabel = QtWidgets.QLabel("Resolution:")
        self.gridLayout.addWidget(self.resolutionLabel, 2, 0, 1, 1)
        self.svresCombox = QtWidgets.QComboBox(self.inputPage)
        self.gridLayout.addWidget(self.svresCombox, 2, 1, 1, 2) 
        ## running HiSV
        self.runHiSVButton = QtWidgets.QPushButton(self.inputPage)
        self.runHiSVButton.setObjectName("runHiSVButton")
        self.runHiSVButton.setFixedSize(150, 30)
        self.runHiSVButton.setEnabled(False)
        self.gridLayout.addWidget(self.runHiSVButton, 3, 0, 1, 3)
        ## HiSV result shown
        self.HiSV_result = QtWidgets.QTextEdit(self.inputPage)
        self.HiSV_result.setReadOnly(True)
        self.HiSV_result.setObjectName("HiSV_result")
        self.HiSV_result.setPlaceholderText("The result of HiSV...")
        font = QFont()
        font.setPointSize(12)
        self.HiSV_result.setFont(font)
        self.HiSV_result.setFixedWidth(450)
        self.HiSV_result.setMaximumHeight(100)
        self.HiSV_result.setAlignment(QtCore.Qt.AlignJustify)
        self.HiSV_result.setWordWrapMode(QTextOption.WrapAtWordBoundaryOrAnywhere)
        self.gridLayout.addWidget(self.HiSV_result, 4, 0, 1, 3)
        ## input SV event
        ### pos1
        self.position1Label = QtWidgets.QLabel("Position 1:")
        self.position1Input = QtWidgets.QLineEdit()
        self.position1Input.setPlaceholderText("e.g. chr1:10000-50000.")
        self.gridLayout.addWidget(self.position1Label, 5, 0, 1, 1)
        self.gridLayout.addWidget(self.position1Input, 5, 1, 1, 2) 
        ### pos2
        self.position2Label = QtWidgets.QLabel("Position 2:")
        self.position2Input = QtWidgets.QLineEdit()
        self.position2Input.setPlaceholderText("e.g. chr1:10000-50000.")
        self.gridLayout.addWidget(self.position2Label, 6, 0, 1, 1)
        self.gridLayout.addWidget(self.position2Input, 6, 1, 1, 2)
        ## Interaction count shown
        self.Interaction_show =  QtWidgets.QTextEdit(self.inputPage)
        self.Interaction_show.setReadOnly(True)
        self.Interaction_show.setPlaceholderText("Interaction counts in the clicked region")
        font = QFont()
        font.setPointSize(12)
        self.Interaction_show.setFont(font)
        self.Interaction_show.setFixedWidth(450)
        self.Interaction_show.setMaximumHeight(100)
        self.Interaction_show.setAlignment(QtCore.Qt.AlignJustify)
        self.Interaction_show.setWordWrapMode(QTextOption.WrapAtWordBoundaryOrAnywhere)
        self.gridLayout.addWidget(self.Interaction_show, 7, 0, 1, 3)
        ## submit
        self.pageSubmitButton = QtWidgets.QPushButton(self.inputPage)
        self.pageSubmitButton.setObjectName("pageSubmitButton")
        self.pageSubmitButton.setFixedSize(150, 30)
        self.pageSubmitButton.setEnabled(False)
        self.gridLayout.addWidget(self.pageSubmitButton, 8, 0, 1, 1)
        ## reset
        self.pageResetButton = QtWidgets.QPushButton(self.inputPage)
        self.pageResetButton.setObjectName("pageResetButton")
        self.pageResetButton.setEnabled(False)
        self.pageResetButton.setFixedSize(150, 30)
        self.gridLayout.addWidget(self.pageResetButton, 8, 2, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem1, 9, 0, 1, 2)
        self.resetButton = QtWidgets.QPushButton(self.inputPage)
        self.resetButton.setObjectName("resetButton")
        self.resetButton.setFixedSize(150, 30)
        self.resetButton.setText("Reset")
        self.gridLayout.addWidget(self.resetButton, 10, 2, 1, 1)
        self.resetButton.clicked.connect(self.resetPage)
        # inputPage1: assemble CSV based brealpoints
        self.inputPage1 = QtWidgets.QWidget()
        self.inputPage1.setObjectName("inputPage1")
        self.gridLayout1 = QtWidgets.QGridLayout(self.inputPage1)
        self.gridLayout1.setObjectName("gridLayout1")
        ## tumor sample input
        self.tumorFolder = QtWidgets.QLineEdit(self.inputPage1)
        self.tumorFolder.setEnabled(True)
        self.tumorFolder.setReadOnly(True)
        self.tumorFolder.setObjectName("tumorFolder")
        self.tumorFolder.setPlaceholderText("Select tumor sample (cool format).")
        self.gridLayout1.addWidget(self.tumorFolder, 0, 0, 1, 2)
        self.tumorFolderButton = QtWidgets.QPushButton(self.inputPage1)
        self.tumorFolderButton.setObjectName("tumorFolderButton")
        self.tumorFolderButton.setFixedSize(150, 30)
        self.gridLayout1.addWidget(self.tumorFolderButton, 0, 2, 1, 1)
        # load
        self.LoadButton = QtWidgets.QPushButton(self.inputPage1)
        self.LoadButton.setObjectName("LoadButton")
        self.LoadButton.setFixedSize(150, 30)
        self.gridLayout1.addWidget(self.LoadButton, 1, 2, 1, 1)
        ## rearrangement input
        self.svFile = QtWidgets.QLineEdit(self.inputPage1)
        self.svFile.setEnabled(True)
        self.svFile.setReadOnly(True)
        self.svFile.setObjectName("svFile")
        self.svFile.setPlaceholderText("Select a breakpoints file.")
        self.gridLayout1.addWidget(self.svFile, 3, 0, 1, 2)
        self.svFileButton = QtWidgets.QPushButton(self.inputPage1)
        self.svFileButton.setObjectName("svFileButton")
        self.svFileButton.setFixedSize(150, 30)
        self.svFileButton.setEnabled(False)
        self.gridLayout1.addWidget(self.svFileButton, 3, 2, 1, 1)
        ## resolution setting
        self.label = QtWidgets.QLabel("Resolution:")
        self.gridLayout1.addWidget(self.label, 2, 0, 1, 1)
        self.resCombox = QtWidgets.QComboBox(self.inputPage1)
        self.gridLayout1.addWidget(self.resCombox, 2, 1, 1, 2) 
        self.chromlist_show = QtWidgets.QTextEdit(self.inputPage1)
        self.chromlist_show.setReadOnly(True)
        self.chromlist_show.setObjectName("chromlist_show")
        self.chromlist_show.setPlaceholderText("The chromosomes of rearrangement file...")
        font = QFont()
        font.setPointSize(12)
        self.chromlist_show.setFont(font)
        self.chromlist_show.setFixedWidth(450)
        self.chromlist_show.setMaximumHeight(50)
        self.chromlist_show.setAlignment(QtCore.Qt.AlignJustify)
        self.chromlist_show.setWordWrapMode(QTextOption.WrapAtWordBoundaryOrAnywhere)
        self.gridLayout1.addWidget(self.chromlist_show, 4, 0, 1, 3)
        ## run cnv normaliztion
        self.pre_cnvButton = QtWidgets.QPushButton(self.inputPage1)
        self.pre_cnvButton.setObjectName("cnvButton")
        self.pre_cnvButton.setEnabled(False)
        self.pre_cnvButton.setFixedSize(150, 30)
        self.gridLayout1.addWidget(self.pre_cnvButton, 5, 2, 1, 1)
        # run assembly
        self.page1SubmitButton = QtWidgets.QPushButton(self.inputPage1)
        self.page1SubmitButton.setObjectName("page1SubmitButton")
        self.page1SubmitButton.setFixedSize(150, 30)
        self.page1SubmitButton.setEnabled(False)
        self.gridLayout1.addWidget(self.page1SubmitButton, 6, 2, 1, 1)
        # assemled CGR result
        self.CGR_result = QtWidgets.QTextEdit(self.inputPage)
        self.CGR_result.setReadOnly(True)
        self.CGR_result.setObjectName("HiSV_result")
        self.CGR_result.setPlaceholderText("The result of assembled CGRs.")
        font = QFont()
        font.setPointSize(12)
        self.CGR_result.setFont(font)
        self.CGR_result.setFixedWidth(450)
        self.CGR_result.setMaximumHeight(100)
        self.CGR_result.setAlignment(QtCore.Qt.AlignJustify)
        self.CGR_result.setWordWrapMode(QTextOption.WrapAtWordBoundaryOrAnywhere)
        self.gridLayout1.addWidget(self.CGR_result, 7, 0, 1, 3)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout1.addItem(spacerItem, 8, 0, 1, 3)
        self.page1ResetButton = QtWidgets.QPushButton(self.inputPage1)
        self.page1ResetButton.setObjectName("page1ResetButton")
        # self.page1ResetButton.setEnabled(False)
        self.page1ResetButton.setFixedSize(150, 30)
        self.page1ResetButton.clicked.connect(self.resetPage1)
        self.gridLayout1.addWidget(self.page1ResetButton, 9, 2, 1, 1)
        # inputPage2
        self.inputPage2 = QtWidgets.QWidget()
        self.inputPage2.setObjectName("inputPage2")
        self.gridLayout2 = QtWidgets.QGridLayout(self.inputPage2)
        self.gridLayout2.setObjectName("gridLayout2")
        ## tumor folder
        self.simutumorFolder = QtWidgets.QLineEdit(self.inputPage2)
        self.simutumorFolder.setEnabled(True)
        self.simutumorFolder.setReadOnly(True)
        self.simutumorFolder.setObjectName("simutumorFolder")
        self.simutumorFolder.setPlaceholderText("Select tumor sample (cool format).")
        self.gridLayout2.addWidget(self.simutumorFolder, 0, 0, 1, 2)
        self.SimutumorFolderButton = QtWidgets.QPushButton(self.inputPage2)
        self.SimutumorFolderButton.setObjectName("SimutumorFolderButton")
        self.SimutumorFolderButton.setFixedSize(150, 30)
        self.gridLayout2.addWidget(self.SimutumorFolderButton, 0, 2, 1, 1)
        ## control folder
        self.conFolder = QtWidgets.QLineEdit(self.inputPage2)
        self.conFolder.setEnabled(True)
        self.conFolder.setReadOnly(True)
        self.conFolder.setObjectName("conFolder")
        self.conFolder.setPlaceholderText("Select control sample (cool format).")
        self.gridLayout2.addWidget(self.conFolder, 1, 0, 1, 2)
        self.conFolderButton = QtWidgets.QPushButton(self.inputPage2)
        self.conFolderButton.setObjectName("conFolderButton")
        self.conFolderButton.setFixedSize(150, 30)
        self.gridLayout2.addWidget(self.conFolderButton, 1, 2, 1, 1)
        self.simuLoadButton = QtWidgets.QPushButton(self.inputPage2)
        self.simuLoadButton.setObjectName("simuLoadButton")
        self.simuLoadButton.setFixedSize(150, 30)
        self.gridLayout2.addWidget(self.simuLoadButton, 2, 2, 1, 1)
        ## resolution
        self.Simulabel = QtWidgets.QLabel("Resolution:")
        self.gridLayout2.addWidget(self.Simulabel, 3, 0, 1, 1)
        self.simuresCombox = QtWidgets.QComboBox(self.inputPage2)
        self.gridLayout2.addWidget(self.simuresCombox, 3, 1, 1, 2) 
        ## chrom list
        self.simuchromlist_show = QtWidgets.QTextEdit(self.inputPage1)
        self.simuchromlist_show.setReadOnly(True)
        self.simuchromlist_show.setObjectName("simuchromlist_show")
        self.simuchromlist_show.setPlaceholderText("List of CGR events.")
        font = QFont()
        font.setPointSize(12)
        self.simuchromlist_show.setFont(font)
        self.simuchromlist_show.setFixedWidth(450)
        self.simuchromlist_show.setMaximumHeight(50)
        self.simuchromlist_show.setAlignment(QtCore.Qt.AlignJustify)
        self.simuchromlist_show.setWordWrapMode(QTextOption.WrapAtWordBoundaryOrAnywhere)
        self.gridLayout2.addWidget(self.simuchromlist_show, 5, 0, 1, 3)
        ## Assembled CGR
        self.assemblyFile = QtWidgets.QLineEdit(self.inputPage2)
        self.assemblyFile.setEnabled(True)
        self.assemblyFile.setReadOnly(True)
        self.assemblyFile.setObjectName("assemblyFile")
        self.assemblyFile.setPlaceholderText("Select a assembled CGR file.")
        self.gridLayout2.addWidget(self.assemblyFile, 4, 0, 1, 2)
        self.assemblyFileButton = QtWidgets.QPushButton(self.inputPage2)
        self.assemblyFileButton.setObjectName("assemblyFileButton")
        self.assemblyFileButton.setEnabled(False)
        self.assemblyFileButton.setFixedSize(150, 30)
        self.gridLayout2.addWidget(self.assemblyFileButton, 4, 2, 1, 1)
        ## simulation process
        self.page2SubmitButton = QtWidgets.QPushButton(self.inputPage2)
        self.page2SubmitButton.setObjectName("page2SubmitButton")
        self.page2SubmitButton.setEnabled(False)
        self.gridLayout2.addWidget(self.page2SubmitButton, 9, 0, 1, 1)
        self.page2visualResetButton = QtWidgets.QPushButton(self.inputPage2)
        self.page2visualResetButton.setObjectName("page2SubmitButton")
        self.page2visualResetButton.setEnabled(False)
        self.page2visualResetButton.setFixedSize(150, 30)
        self.gridLayout2.addWidget(self.page2visualResetButton, 9, 2, 1, 1)
        self.page2visualResetButton.clicked.connect(self.visualPage2_reset)
        self.page2SubmitButton.setFixedSize(150, 30)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout2.addItem(spacerItem1, 10, 0, 1, 2)
        self.simuInteraction_show =  QtWidgets.QLineEdit()
        self.simuInteraction_show.setReadOnly(True)
        self.simuInteraction_show.setPlaceholderText("Interaction counts in the clicked region")
        font = QFont()
        font.setPointSize(12)
        self.simuInteraction_show.setFont(font)
        self.gridLayout2.addWidget(self.simuInteraction_show,11, 0, 1, 3)
        self.page2ResetButton = QtWidgets.QPushButton(self.inputPage2)
        self.page2ResetButton.setObjectName("page2ResetButton")
        self.page2ResetButton.clicked.connect(self.resetPage2)
        self.page2ResetButton.setFixedSize(150, 30)
        self.gridLayout2.addWidget(self.page2ResetButton, 12, 2, 1, 1)
        # inputPage3
        self.inputPage3 = QtWidgets.QWidget()
        self.inputPage3.setObjectName("inputPage3")
        self.gridLayout3 = QtWidgets.QGridLayout(self.inputPage3)
        self.gridLayout3.setObjectName("gridLayout3")
        self.rehicFolder = QtWidgets.QLineEdit(self.inputPage3)
        self.rehicFolder.setEnabled(True)
        self.rehicFolder.setReadOnly(True)
        self.rehicFolder.setObjectName("hicFolder")
        self.rehicFolder.setPlaceholderText("Select tumor sample (cool format).")
        self.gridLayout3.addWidget(self.rehicFolder, 0, 0, 1, 2)
        self.rehicFolderButton = QtWidgets.QPushButton(self.inputPage3)
        self.rehicFolderButton.setObjectName("rehicFolderButton")
        self.rehicFolderButton.setFixedSize(150, 30)
        self.gridLayout3.addWidget(self.rehicFolderButton, 0, 2, 1, 1)
        self.reLoadButton = QtWidgets.QPushButton(self.inputPage3)
        self.reLoadButton.setObjectName("reLoadButton")
        self.reLoadButton.setFixedSize(150, 30)
        self.gridLayout3.addWidget(self.reLoadButton, 1, 2, 1, 1)
        ## resolution setting
        self.reresolutionLabel = QtWidgets.QLabel("Resolution:")
        self.gridLayout3.addWidget(self.reresolutionLabel, 2, 0, 1, 1)
        self.reresCombox = QtWidgets.QComboBox(self.inputPage3)
        self.gridLayout3.addWidget(self.reresCombox, 2, 1, 1, 2) 
        ##runCNV_normalization
        self.cnvButton = QtWidgets.QPushButton(self.inputPage3)
        self.cnvButton.setObjectName("cnvButton")
        self.cnvButton.setEnabled(False)
        self.cnvButton.setFixedSize(150, 30)
        self.gridLayout3.addWidget(self.cnvButton, 5, 0, 1, 1)
        ## reconstruct result
        self.reAssemblyFileButton = QtWidgets.QPushButton(self.inputPage3)
        self.reAssemblyFileButton.setObjectName("reAssemblyFileButton")
        self.reAssemblyFileButton.setFixedSize(150, 30)
        self.reAssemblyFileButton.setEnabled(False)
        self.gridLayout3.addWidget(self.reAssemblyFileButton, 3, 2, 1, 1)
        self.reAssemblyFile = QtWidgets.QLineEdit(self.inputPage3)
        self.reAssemblyFile.setPlaceholderText("Select a CGR event file.")
        self.reAssemblyFile.setReadOnly(True)
        self.reAssemblyFile.setObjectName("reAssemblyFile")
        self.gridLayout3.addWidget(self.reAssemblyFile, 3, 0, 1, 2)
        self.fragment_show = QtWidgets.QTextEdit(self.inputPage3)
        self.fragment_show.setReadOnly(True)
        self.fragment_show.setObjectName("chromlist_show")
        self.fragment_show.setPlaceholderText("The fragment of the CGR event.")
        font = QFont()
        font.setPointSize(12)
        self.fragment_show.setFont(font)
        self.fragment_show.setFixedWidth(450)
        self.fragment_show.setMaximumHeight(80)
        self.fragment_show.setAlignment(QtCore.Qt.AlignJustify)
        self.fragment_show.setWordWrapMode(QTextOption.WrapAtWordBoundaryOrAnywhere)
        self.gridLayout3.addWidget(self.fragment_show, 4, 0, 1, 3)
        self.page3SubmitButton = QtWidgets.QPushButton(self.inputPage3)
        self.page3SubmitButton.setObjectName("page3SubmitButton")
        self.page3SubmitButton.setFixedSize(150, 30)
        self.page3SubmitButton.setEnabled(False)
        self.gridLayout3.addWidget(self.page3SubmitButton, 5, 2, 1, 1)
        self.exportResult = QtWidgets.QPushButton(self.inputPage3)
        self.exportResult.setObjectName("exportResult")
        self.exportResult.setFixedSize(150, 30)
        self.exportResult.setEnabled(False)
        self.gridLayout3.addWidget(self.exportResult, 7, 0, 1, 1)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout3.addItem(spacerItem2, 6, 0, 1, 2)
        self.page3ResetButton = QtWidgets.QPushButton(self.inputPage3)
        self.page3ResetButton.setObjectName("page3ResetButton")
        self.page3ResetButton.setFixedSize(150, 30)
        self.page3ResetButton.clicked.connect(self.resetPage3)
        self.gridLayout3.addWidget(self.page3ResetButton, 7, 2, 1, 1)
        # showPage
        self.page = QtWidgets.QWidget()
        self.page.setObjectName("page")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.page)
        self.verticalLayout.setObjectName("verticalLayout")
        # showPage1
        self.page_1 = QtWidgets.QWidget()
        self.page_1.setObjectName("page")
        self.verticalLayout_1 = QtWidgets.QVBoxLayout(self.page_1)
        self.verticalLayout_1.setObjectName("verticalLayout")
        # showPage2
        self.page_2 = QtWidgets.QWidget()
        self.page_2.setObjectName("page_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.page_2)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        # showPage3
        self.page_3 = QtWidgets.QWidget()
        self.page_3.setObjectName("page_3")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.page_3)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        # page
        self.hicFolderButton.setText(_translate("Form", "Tumor sample"))
        self.runHiSVButton.setText(_translate("Form", "Run HiSV"))
        self.pageSubmitButton.setText(_translate("Form", "Visualization"))
        self.pageResetButton.setText(_translate("Form", "Reset visualization"))
        self.hicsvLoadButton.setText(_translate("Form", "Load File"))
        # page1
        self.tumorFolderButton.setText(_translate("Form", "Tumor sample"))
        self.svFileButton.setText(_translate("Form", "Breakpoint File"))
        self.page1SubmitButton.setText(_translate("Form", "Fragment Assembly"))
        self.LoadButton.setText(_translate("Form", "Load file"))
        self.pre_cnvButton.setText(_translate("Form", "CNV normalization"))
        self.page1ResetButton.setText(_translate("Form", "Reset"))
        # page2
        self.SimutumorFolderButton.setText(_translate("Form", "Tumor Sample"))
        self.simuLoadButton.setText(_translate("Form", "Load file"))
        self.Simulabel.setText(_translate("Form", "Resolution"))
        self.assemblyFileButton.setText(_translate("Form", "Assembly CGRs File"))
        self.page2SubmitButton.setText(_translate("Form", "Validation"))
        self.conFolderButton.setText(_translate("Form", "Control Sample"))
        self.page2ResetButton.setText(_translate("Form", "Reset"))
        self.page2visualResetButton.setText(_translate("Form", "Reset visualization"))
        # page3
        self.rehicFolderButton.setText(_translate("Form", "Tumor Sample"))
        self.reLoadButton.setText(_translate("Form", "Load file"))
        self.reAssemblyFileButton.setText(_translate("Form", "CGR File"))
        self.cnvButton.setText(_translate("Form", "CNV normalization"))
        self.page3SubmitButton.setText(_translate("Form", "Reconstruct"))
        self.exportResult.setText(_translate("Form", "Export result"))
        self.page3ResetButton.setText(_translate("Form", "Reset"))

    def open_user_manual(self):
        QDesktopServices.openUrl(QUrl("https://decgr.readthedocs.io/en/latest/"))

    def download_test_data(self):
        QDesktopServices.openUrl(QUrl("https://github.com/GaoLabXDU/DeCGR/tree/main/TestData"))


if __name__ == "__main__":
    QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)
    app = QApplication(sys.argv)
    w = Interface()
    w.show()
    sys.exit(app.exec_())
