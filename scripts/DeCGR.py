#!/usr/bin/env python

from PyQt5.QtWidgets import (
    QMainWindow, QLabel, QToolBar, QAction, QWidget, QVBoxLayout, QStackedWidget, QPushButton, QToolButton, QSplitter,
    QHBoxLayout, QTextEdit, QFrame, QApplication, QSizePolicy
)
from PyQt5.QtGui import QPixmap, QTextOption, QFont, QPainter, QColor, QDesktopServices
from PyQt5.QtCore import QUrl, Qt
import sys
from PyQt5 import QtCore, QtGui, QtWidgets
from iced import normalization
from decgr.DetermineFragment import FragmentAssembly
from decgr.simuProcess import simulation, get_cut_pos, plot_simu_matrix
from decgr.plot_reconstruct_map import global_expect, joint_SV, adjust_matrix, plot_assembly_matrix, read_one_assembly
from decgr.utils import *
import cooler
from decgr import correct


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
        self.path = os.path.dirname(self.path)
        if not os.path.exists(self.path + '/result'):
            os.mkdir(self.path + '/result')
        if not os.path.exists(self.path + '/tmp'):
            os.mkdir(self.path + '/tmp')

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

        self.setStyleSheet('QAbstractButton {font-size: 20px;} .QLabel {font-size: 20px;} \
                     .QLineEdit {font-size: 14px;} #toLabel {font-size: 48px;} .QComboBox {font-size: 14px;} \
                     .QComboBox:disabled {color: rgb(0,0,0)} .QListWidget {font-size: 14px;}')
        # module 1
        self.mainWidget = QWidget()
        self.splitter = QSplitter(QtCore.Qt.Horizontal)
        self.splitter.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.splitter.addWidget(self.inputPage1)
        self.splitter.addWidget(self.page)
        self.splitter.setSizes([self.width // 3, self.width * 2 // 3])
        mainLayout = QHBoxLayout(self.mainWidget)
        mainLayout.addWidget(self.splitter)

        # module 2
        self.mainWidget_2 = QWidget()
        self.splitter_2 = QSplitter(QtCore.Qt.Horizontal)
        self.splitter_2.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.splitter_2.addWidget(self.inputPage2)
        self.splitter_2.addWidget(self.page_2)
        self.splitter_2.setSizes([self.width // 3, self.width * 2 // 3])
        mainLayout_2 = QHBoxLayout(self.mainWidget_2)
        mainLayout_2.addWidget(self.splitter_2)
        # module 3
        self.mainWidget_3 = QWidget()
        self.splitter_3 = QSplitter(QtCore.Qt.Horizontal)
        self.splitter_3.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.splitter_3.addWidget(self.inputPage3)
        self.splitter_3.addWidget(self.page_3)
        self.splitter_3.setSizes([self.width // 3, self.width * 2 // 3])
        mainLayout_3 = QHBoxLayout(self.mainWidget_3)
        mainLayout_3.addWidget(self.splitter_3)

        self.stackedWidget.addWidget(self.inputPage)
        self.stackedWidget.addWidget(self.mainWidget)
        self.stackedWidget.addWidget(self.mainWidget_2)
        self.stackedWidget.addWidget(self.mainWidget_3)

        toolbar = CustomToolBar("My main toolbar")
        toolbar.setFixedHeight(50)
        toolbar.setStyleSheet("background-color: rgb(53, 58, 64); color: white; font-size: 20px;")
        self.addToolBar(toolbar)
        button_action0 = QAction("DeCGR", self)
        button_action0.setStatusTip("Go to MainPage")
        button_action0.triggered.connect(lambda: self.stackedWidget.setCurrentIndex(0))
        tool_button = QToolButton()
        tool_button.setDefaultAction(button_action0)
        toolbar.addWidget(tool_button)
        font = QFont()
        font.setBold(True)
        tool_button.setFont(font)
        tool_button.setStyleSheet("color: rgb(245, 195, 66);")

        self.button_action1 = QAction("Fragment Assembly", self)
        self.button_action1.setStatusTip("Go to Page 1")
        self.button_action1.triggered.connect(lambda:  self.stackedWidget.setCurrentIndex(1))
        toolbar.addAction(self.button_action1)

        self.button_action2 = QAction("Validation CGRs", self)
        self.button_action2.setStatusTip("Go to Page 2")
        self.button_action2.triggered.connect(lambda: self.stackedWidget.setCurrentIndex(2))
        toolbar.addAction(self.button_action2)

        self.button_action3 = QAction("Reconstruct Hi-C Contact Map", self)
        self.button_action3.setStatusTip("Go to Page 3")
        self.button_action3.triggered.connect(lambda: self.stackedWidget.setCurrentIndex(3))
        # self.button_action3.setEnabled(False)
        toolbar.addAction(self.button_action3)

        self.tumorFolderButton.clicked.connect(lambda: self.inputFilePath('tumor'))
        self.svFileButton.clicked.connect(lambda: self.inputFilePath('sv'))
        re = QtCore.QRegExp('^[1-9][0-9]*$')
        validator = QtGui.QRegExpValidator(re)
        self.resolution.setValidator(validator)
        self.addButton.clicked.connect(self.addItem)
        self.delButton.clicked.connect(self.delItem)
        self.page1SubmitButton.clicked.connect(self.step1)

        self.conFolderButton.clicked.connect(lambda: self.inputFilePath('control'))
        self.SimutumorFolderButton.clicked.connect(lambda: self.inputFilePath('simutumor'))
        self.assemblyFileButton.clicked.connect(lambda: self.inputFilePath('assembly'))
        self.simuresolution.setValidator(validator)
        self.simuaddButton.clicked.connect(self.simuaddItem)
        self.simudelButton.clicked.connect(self.simudelItem)
        self.useButton.toggled.connect(self.fithic)
        self.dis2cntButton.clicked.connect(lambda: self.inputFilePath('dis2cnt'))
        self.page2SubmitButton.clicked.connect(self.step2)

        self.cnvButton.toggled.connect(self.cnv)
        self.reAssemblyFileButton.clicked.connect(lambda: self.inputFilePath('re'))
        self.page3SubmitButton.clicked.connect(self.step3)

    def inputFilePath(self, tag):
        if tag == 'tumor':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', self.path)[0]
            if path == '':
                return
            self.tumorFolder.setText(path)
        elif tag == 'simutumor':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', self.path)[0]
            if path == '':
                return
            self.simutumorFolder.setText(path)
        elif tag == 'sv':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', self.path)[0]
            if path == '':
                return
            self.svFile.setText(path)
        elif tag == 'control':
            # path = QtWidgets.QFileDialog.getExistingDirectory(self, 'Input Directory Path', self.path)
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', self.path)[0]
            if path == '':
                return
            self.conFolder.setText(path)
        elif tag == 'assembly':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', self.path)[0]
            if path == '':
                return
            self.assemblyFile.setText(path)
        elif tag == 'dis2cnt':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', self.path)[0]
            if path == '':
                return
            self.dis2cntFile.setText(path)
        elif tag == 're':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', self.path)[0]
            if path == '':
                return
            self.reAssemblyFile.setText(path)

    def addItem(self):
        chro = self.chrCombox.currentText()
        if chro == 'Chromosome List':
            return
        for i in range(self.chrList.count()):
            item = self.chrList.item(i).text()
            if chro == item:
                return
        self.chrList.addItem(chro)

    def delItem(self):
        curRow = self.chrList.currentRow()
        self.chrList.takeItem(curRow)

    def simuaddItem(self):
        chro = self.simuchrCombox.currentText()
        if chro == 'Chromosome List':
            return
        for i in range(self.simuchrList.count()):
            item = self.simuchrList.item(i).text()
            if chro == item:
                return
        self.simuchrList.addItem(chro)

    def simudelItem(self):
        curRow = self.simuchrList.currentRow()
        self.simuchrList.takeItem(curRow)

    def fithic(self, checked):
        if checked:
            self.dis2cntButton.setEnabled(True)
            self.dis2cntFile.setEnabled(True)
        else:
            self.dis2cntButton.setEnabled(False)
            self.dis2cntFile.setEnabled(False)

    def cnv(self, checked):
        if checked:
            self.genome.setEnabled(True)
            self.enzyme.setEnabled(True)
        else:
            self.genome.setEnabled(False)
            self.enzyme.setEnabled(False)

    def step1(self):
        SV_file = self.svFile.text()
        self.res = int(self.resolution.text())
        hic_path = self.tumorFolder.text()
        if 'mcool' in hic_path:
            hic_file = hic_path + '::/resolutions/' + str(self.res)
        else:
            hic_file = hic_path
        clr = cooler.Cooler(hic_file)
        self.chrom_list = []
        for i in range(self.chrList.count()):
            item = self.chrList.item(i).text()
            self.chrom_list.append(item)
        self.chrom_list = sorted(self.chrom_list, key=lambda x: int(x[3:]))
        self.dict_matrix = {}
        self.dict_count = {}
        self.dict_norm = {}
        for i in range(len(self.chrom_list)):
            for j in range(i, len(self.chrom_list)):
                name = str(self.chrom_list[i]) + '_' + str(self.chrom_list[j])
                # cur_matrix = get_MatrixFile(hic_path, name)
                chr_pair = [self.chrom_list[i], self.chrom_list[j]]
                cur_matrix = load_matrix(clr, chr_pair)
                self.dict_count[self.chrom_list[i]] = np.shape(cur_matrix)[0]
                self.dict_matrix[name] = cur_matrix
                if i == j:
                    normed = normalization.ICE_normalization(cur_matrix)
                    self.dict_norm[name] = normed
                else:
                    self.dict_norm[name] = cur_matrix

        all_path, all_orient, chrom, start, end, canvas_matrix, canvas_link, canvas_assembly = FragmentAssembly(SV_file,
                                            self.res, self.dict_norm, self.chrom_list, self.dict_count, self.path)
        result_file = self.path + "/result/assembly_result.txt"
        with open(result_file, 'w') as out:
            out.write("num" + '\t' + 'chrom' + '\t' + "start" + '\t' + "end" + '\t' + "node" + '\t' + "orient" + '\n')
            for i in range(len(all_path)):
                cur_path = all_path[i]
                for j in range(len(cur_path)):
                    index = ord(cur_path[j]) - 65
                    out.write(
                        str(i + 1) + '\t' + chrom[index] + '\t' + str(start[index] * self.res) + '\t' +
                        str(end[index] * self.res) + '\t' + cur_path[j] + '\t' + all_orient[i][j] + '\n')
        self.hicmap1 = canvas_matrix
        self.link = canvas_link
        self.assembly = canvas_assembly
        self.plot1()

    def plot1(self):
        parent_width = self.page.width()
        parent_height = self.page.height()
        assembly_width = int(parent_width * 0.2)
        assembly_height = int(parent_height * 0.2)
        self.assembly.setGeometry(0, 0, assembly_width, assembly_height)
        self.verticalLayout.addWidget(self.assembly)
        # self.assembly.setFixedSize(assembly_width, assembly_height)
        hicmap1_height = int(parent_height * 0.6)
        link_height = int(parent_height * 0.2)
        self.verticalLayout.addWidget(self.hicmap1)
        self.verticalLayout.addWidget(self.link)
        self.hicmap1.setFixedSize(parent_width, hicmap1_height)
        self.link.setFixedSize(parent_width, link_height)

    def step2(self):
        result_file = self.assemblyFile.text()
        self.res = int(self.simuresolution.text())
        hic_path = self.simutumorFolder.text()
        if 'mcool' in hic_path:
            hic_file = hic_path + '::/resolutions/' + str(self.res)
        else:
            hic_file = hic_path
        clr = cooler.Cooler(hic_file)
        self.simuchrom_list = []
        for i in range(self.simuchrList.count()):
            item = self.simuchrList.item(i).text()
            self.simuchrom_list.append(item)
        self.simuchrom_list = sorted(self.simuchrom_list, key=lambda x: int(x[3:]))
        self.dict_matrix = {}
        self.dict_count = {}
        self.dict_norm = {}
        control_path = self.conFolder.text()
        control_file = control_path + '::/resolutions/' + str(self.res)
        control_clr = cooler.Cooler(control_file)
        control_dict = {}
        control_len = []
        for i in range(len(self.simuchrom_list)):
            for j in range(i, len(self.simuchrom_list)):
                name = str(self.simuchrom_list[i]) + '_' + str(self.simuchrom_list[j])
                # cur_matrix = get_MatrixFile(control_path, name)
                chr_pair = [self.simuchrom_list[i], self.simuchrom_list[j]]
                cur_control_matrix = load_matrix(control_clr, chr_pair)
                control_dict[name] = cur_control_matrix
                cur_tumor_matrix = load_matrix(clr, chr_pair)
                self.dict_matrix[name] = cur_tumor_matrix
                self.dict_count[self.simuchrom_list[i]] = np.shape(cur_tumor_matrix)[0]
                if i == j:
                    control_len.append(cur_control_matrix.shape[0])
                    normed = normalization.ICE_normalization(cur_tumor_matrix)
                    self.dict_norm[name] = normed
                else:
                    self.dict_norm[name] = cur_control_matrix
        if self.runButton.isChecked():
            dis2cnt = runFitHiC(self.path, control_dict, self.simuchrom_list, control_len, self.res)
        else:
            dis2cnt = self.dis2cntFile.text()
        all_path, all_orient, chrom, start, end = read_assembly_result(result_file, self.res)
        simu_dict, rows_dict, cols_dict = simulation(all_path, all_orient, chrom, start, end, dis2cnt, control_dict,
                                                     self.simuchrom_list,
                                                     self.res)
        for key in simu_dict:
            SV_matrix = simu_dict[key]
            filename = self.path + "/result/simu_" + key + '_matrix.txt'
            np.savetxt(filename, SV_matrix, fmt='%.2f', delimiter='\t')
        dict_start, dict_end = get_cut_pos(self.simuchrom_list, chrom, start, end, self.res, self.dict_count)
        print(dict_start, dict_end, self.simuchrom_list)
        canvas = plot_simu_matrix(self.dict_norm, simu_dict, dict_start, dict_end, self.simuchrom_list, rows_dict,
                                  cols_dict)
        self.simuMap = canvas
        self.plot2()
        self.button_action3.setEnabled(True)

    def plot2(self):
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.addWidget(self.simuMap)
        print(self.simuMap.width(), self.simuMap.height())
        simuLabel = QtWidgets.QLabel('Simulated', self.simuMap)
        simuLabel.move(int(self.simuMap.width() * 0.6), int(self.simuMap.height() * 0.15))
        font = QtGui.QFont('Arial', 26, 75)
        simuLabel.setFont(font)
        simuLabel.show()
        varLabel = QtWidgets.QLabel('Original', self.simuMap)
        varLabel.move(int(self.simuMap.width() * 0.6), int(self.simuMap.height() * 1.5))
        varLabel.setFont(font)
        varLabel.show()

    def step3(self):
        result_file = self.reAssemblyFile.text()
        path, orient, chrom, start, end = read_one_assembly(result_file, self.res)
        hic_path = self.simutumorFolder.text()
        if 'mcool' in hic_path:
            hic_file = hic_path + '::/resolutions/' + str(self.res)
        else:
            hic_file = hic_path
        clr = cooler.Cooler(hic_file)
        # ICE or sweight
        if self.cnvButton.isChecked():
            self.dict_norm = {}
            # running correct-CNV
            genome = self.genome.text()
            enzyme = self.enzyme.text()
            correct.run(hic_file, genome, enzyme, self.res)
            # load matrix
            for i in range(len(self.simuchrom_list)):
                for j in range(i, len(self.simuchrom_list)):
                    name = str(self.simuchrom_list[i]) + '_' + str(self.simuchrom_list[j])
                    chr_pair = [self.simuchrom_list[i], self.simuchrom_list[j]]
                    if i == j:
                        bias = load_bias(clr, chr_pair)
                        M = load_matrix(clr, chr_pair)
                        print(M)
                        Mbias = bias * bias
                        normed = Mbias * M
                        self.dict_norm[name] = normed
                    else:
                        M = load_matrix(clr, chr_pair)
                        self.dict_norm[name] = M

        glob_exp = global_expect(self.dict_norm, self.simuchrom_list)
        SV_matrix = joint_SV(path, orient, chrom, start, end, self.dict_norm)
        matrix, interval = adjust_matrix(SV_matrix, glob_exp, path, start, end)
        canvas = plot_assembly_matrix(SV_matrix, interval, path, orient)
        self.reconstructMap = canvas
        self.plot3()
        self.SV_matrix = SV_matrix
        self.exportResult.clicked.connect(self.saveResult)

    def plot3(self):
        self.horizontalLayout_3.addWidget(self.reconstructMap)

    def saveResult(self):
        filepath, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', '', 'Text Files (*.txt);;All Files (*)')
        if filepath:
            np.savetxt(filepath, self.SV_matrix, delimiter='\t', fmt='%.2f')

    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(1395, 940)
        # instruction
        self.inputPage = QtWidgets.QWidget()
        self.inputPage.setObjectName("inputPage")
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
            "<ul><li>The <b>Fragment Assembly</b> module automatically reconstructs candidate CGRs based on the input breakpoints and the Hi-C map. </li><br/>"
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
        path = os.path.split(os.path.abspath(__file__))[0]
        path = os.path.dirname(path)
        pixmap = QPixmap(path+"/figure/mainplot.jpg")
        image_label.setPixmap(pixmap)
        image_label.setAlignment(QtCore.Qt.AlignCenter)
        right_layout.addWidget(image_label)
        main_layout = QHBoxLayout(self.inputPage)
        main_layout.addLayout(left_layout)
        main_layout.addLayout(right_layout)
        self.inputPage.setLayout(main_layout)
        # inputPage1
        self.inputPage1 = QtWidgets.QWidget()
        self.inputPage1.setObjectName("inputPage1")
        self.gridLayout = QtWidgets.QGridLayout(self.inputPage1)
        self.gridLayout.setObjectName("gridLayout")
        self.addButton = QtWidgets.QPushButton(self.inputPage1)
        self.addButton.setObjectName("addButton")
        self.gridLayout.addWidget(self.addButton, 7, 1, 1, 1)
        self.svFile = QtWidgets.QLineEdit(self.inputPage1)
        self.svFile.setEnabled(True)
        self.svFile.setReadOnly(True)
        self.svFile.setObjectName("svFile")
        self.gridLayout.addWidget(self.svFile, 2, 0, 1, 3)
        self.chrCombox = QtWidgets.QComboBox(self.inputPage1)
        self.chrCombox.setObjectName("chrCombox")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.chrCombox.addItem("")
        self.gridLayout.addWidget(self.chrCombox, 7, 0, 1, 1)
        self.delButton = QtWidgets.QPushButton(self.inputPage1)
        self.delButton.setObjectName("delButton")
        self.gridLayout.addWidget(self.delButton, 7, 2, 1, 1)
        self.resolution = QtWidgets.QLineEdit(self.inputPage1)
        self.resolution.setObjectName("resolution")
        self.gridLayout.addWidget(self.resolution, 5, 1, 1, 2)
        self.tumorFolderButton = QtWidgets.QPushButton(self.inputPage1)
        self.tumorFolderButton.setObjectName("tumorFolderButton")
        self.gridLayout.addWidget(self.tumorFolderButton, 1, 0, 1, 3)
        self.label = QtWidgets.QLabel(self.inputPage1)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 5, 0, 1, 1)
        self.svFileButton = QtWidgets.QPushButton(self.inputPage1)
        self.svFileButton.setObjectName("svFileButton")
        self.gridLayout.addWidget(self.svFileButton, 3, 0, 1, 3)
        self.chrList = QtWidgets.QListWidget(self.inputPage1)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.chrList.sizePolicy().hasHeightForWidth())
        self.chrList.setSizePolicy(sizePolicy)
        self.chrList.setObjectName("chrList")
        self.gridLayout.addWidget(self.chrList, 8, 0, 1, 3)
        self.page1SubmitButton = QtWidgets.QPushButton(self.inputPage1)
        self.page1SubmitButton.setObjectName("page1SubmitButton")
        self.gridLayout.addWidget(self.page1SubmitButton, 9, 1, 1, 2)
        self.tumorFolder = QtWidgets.QLineEdit(self.inputPage1)
        self.tumorFolder.setEnabled(True)
        self.tumorFolder.setReadOnly(True)
        self.tumorFolder.setObjectName("tumorFolder")
        self.gridLayout.addWidget(self.tumorFolder, 0, 0, 1, 3)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem, 10, 0, 1, 3)
        # inputPage2
        self.inputPage2 = QtWidgets.QWidget()
        self.inputPage2.setObjectName("inputPage2")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.inputPage2)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.simutumorFolder = QtWidgets.QLineEdit(self.inputPage2)
        self.simutumorFolder.setEnabled(True)
        self.simutumorFolder.setReadOnly(True)
        self.simutumorFolder.setObjectName("simutumorFolder")
        self.gridLayout_3.addWidget(self.simutumorFolder, 0, 0, 1, 3)
        self.SimutumorFolderButton = QtWidgets.QPushButton(self.inputPage2)
        self.SimutumorFolderButton.setObjectName("SimutumorFolderButton")
        self.gridLayout_3.addWidget(self.SimutumorFolderButton, 1, 0, 1, 3)
        self.conFolder = QtWidgets.QLineEdit(self.inputPage2)
        self.conFolder.setEnabled(True)
        self.conFolder.setReadOnly(True)
        self.conFolder.setObjectName("conFolder")
        self.gridLayout_3.addWidget(self.conFolder, 2, 0, 1, 3)
        self.conFolderButton = QtWidgets.QPushButton(self.inputPage2)
        self.conFolderButton.setObjectName("conFolderButton")
        self.gridLayout_3.addWidget(self.conFolderButton, 3, 0, 1, 3)
        self.Simulabel = QtWidgets.QLabel(self.inputPage2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Simulabel.sizePolicy().hasHeightForWidth())
        self.Simulabel.setSizePolicy(sizePolicy)
        self.Simulabel.setObjectName("Simulabel")
        self.gridLayout_3.addWidget(self.Simulabel, 4, 0, 1, 1)
        self.simuresolution = QtWidgets.QLineEdit(self.inputPage2)
        self.simuresolution.setObjectName("resolution")
        self.gridLayout_3.addWidget(self.simuresolution, 4, 1, 1, 2)
        self.simuchrCombox = QtWidgets.QComboBox(self.inputPage2)
        self.simuchrCombox.setObjectName("simuchrCombox")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.simuchrCombox.addItem("")
        self.gridLayout_3.addWidget(self.simuchrCombox, 5, 0, 1, 1)
        self.simuaddButton = QtWidgets.QPushButton(self.inputPage2)
        self.simuaddButton.setObjectName("addButton")
        self.gridLayout_3.addWidget(self.simuaddButton, 5, 1, 1, 1)
        self.simudelButton = QtWidgets.QPushButton(self.inputPage2)
        self.simudelButton.setObjectName("delButton")
        self.gridLayout_3.addWidget(self.simudelButton, 5, 2, 1, 1)
        self.simuchrList = QtWidgets.QListWidget(self.inputPage2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.simuchrList.sizePolicy().hasHeightForWidth())
        self.simuchrList.setSizePolicy(sizePolicy)
        self.simuchrList.setObjectName("chrList")
        self.gridLayout_3.addWidget(self.simuchrList, 6, 0, 1, 3)
        self.assemblyFileButton = QtWidgets.QPushButton(self.inputPage2)
        self.assemblyFileButton.setObjectName("assemblyFileButton")
        self.gridLayout_3.addWidget(self.assemblyFileButton, 8, 0, 1, 3)
        self.assemblyFile = QtWidgets.QLineEdit(self.inputPage2)
        self.assemblyFile.setEnabled(True)
        self.assemblyFile.setReadOnly(True)
        self.assemblyFile.setObjectName("assemblyFile")
        self.gridLayout_3.addWidget(self.assemblyFile, 7, 0, 1, 3)
        self.page2SubmitButton = QtWidgets.QPushButton(self.inputPage2)
        self.page2SubmitButton.setObjectName("page2SubmitButton")
        self.gridLayout_3.addWidget(self.page2SubmitButton, 11, 2, 1, 1)
        self.useButton = QtWidgets.QRadioButton(self.inputPage2)
        self.useButton.setObjectName("useButton")
        self.gridLayout_3.addWidget(self.useButton, 9, 1, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_3.addItem(spacerItem1, 12, 0, 1, 2)
        self.runButton = QtWidgets.QRadioButton(self.inputPage2)
        self.runButton.setChecked(True)
        self.runButton.setObjectName("runButton")
        self.gridLayout_3.addWidget(self.runButton, 9, 0, 1, 1)
        self.dis2cntButton = QtWidgets.QPushButton(self.inputPage2)
        self.dis2cntButton.setEnabled(False)
        self.dis2cntButton.setObjectName("dis2cntButton")
        self.gridLayout_3.addWidget(self.dis2cntButton, 10, 2, 1, 1)
        self.dis2cntFile = QtWidgets.QLineEdit(self.inputPage2)
        self.dis2cntFile.setEnabled(True)
        self.dis2cntFile.setReadOnly(True)
        self.dis2cntFile.setObjectName("dis2cntFile")
        self.gridLayout_3.addWidget(self.dis2cntFile, 10, 0, 1, 2)
        # inputPage3
        self.inputPage3 = QtWidgets.QWidget()
        self.inputPage3.setObjectName("inputPage3")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.inputPage3)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.iceButton = QtWidgets.QRadioButton(self.inputPage3)
        self.iceButton.setObjectName("iceButton")
        self.gridLayout_4.addWidget(self.iceButton, 0, 1, 1, 1)
        self.cnvButton = QtWidgets.QRadioButton(self.inputPage3)
        self.cnvButton.setChecked(True)
        self.cnvButton.setObjectName("cnvButton")
        self.gridLayout_4.addWidget(self.cnvButton, 0, 0, 1, 1)
        self.genomelabel = QtWidgets.QLabel(self.inputPage3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.genomelabel.sizePolicy().hasHeightForWidth())
        self.genomelabel.setSizePolicy(sizePolicy)
        self.genomelabel.setObjectName("genomelabel")
        self.gridLayout_4.addWidget(self.genomelabel, 1, 0, 1, 1)
        self.genome = QtWidgets.QLineEdit(self.inputPage3)
        self.genome.setObjectName("genome")
        self.gridLayout_4.addWidget(self.genome, 1, 1, 1, 1)
        self.enzymelabel = QtWidgets.QLabel(self.inputPage3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.enzymelabel.sizePolicy().hasHeightForWidth())
        self.enzymelabel.setSizePolicy(sizePolicy)
        self.enzymelabel.setObjectName("enzymelabel")
        self.gridLayout_4.addWidget(self.enzymelabel, 2, 0, 1, 1)
        self.enzyme = QtWidgets.QLineEdit(self.inputPage3)
        self.enzyme.setObjectName("enzyme")
        self.gridLayout_4.addWidget(self.enzyme, 2, 1, 1, 1)
        self.reAssemblyFileButton = QtWidgets.QPushButton(self.inputPage3)
        self.reAssemblyFileButton.setObjectName("reAssemblyFileButton")
        self.gridLayout_4.addWidget(self.reAssemblyFileButton, 4, 0, 1, 2)
        self.label_2 = QtWidgets.QLabel(self.inputPage3)
        self.label_2.setText("")
        self.label_2.setObjectName("label_2")
        self.gridLayout_4.addWidget(self.label_2, 5, 0, 1, 1)
        self.reAssemblyFile = QtWidgets.QLineEdit(self.inputPage3)
        self.reAssemblyFile.setEnabled(True)
        self.reAssemblyFile.setReadOnly(True)
        self.reAssemblyFile.setObjectName("reAssemblyFile")
        self.gridLayout_4.addWidget(self.reAssemblyFile, 3, 0, 1, 2)
        self.page3SubmitButton = QtWidgets.QPushButton(self.inputPage3)
        self.page3SubmitButton.setObjectName("page3SubmitButton")
        self.gridLayout_4.addWidget(self.page3SubmitButton, 5, 1, 1, 1)
        self.exportResult = QtWidgets.QPushButton(self.inputPage3)
        self.exportResult.setObjectName("exportResult")
        self.gridLayout_4.addWidget(self.exportResult, 8, 1, 1, 1)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_4.addItem(spacerItem2, 6, 0, 1, 2)
        # showPage1
        self.page = QtWidgets.QWidget()
        self.page.setObjectName("page")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.page)
        self.verticalLayout.setObjectName("verticalLayout")
        # self.hicMap = QtWidgets.QLabel(self.page)
        # self.hicMap.setText(" ")
        # self.hicMap.setObjectName("hicMap")
        # self.verticalLayout.addWidget(self.hicMap)
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
        # page1
        self.addButton.setText(_translate("Form", "+"))
        self.chrCombox.setItemText(0, _translate("Form", "Chromosome List"))
        self.chrCombox.setItemText(1, _translate("Form", "chr1"))
        self.chrCombox.setItemText(2, _translate("Form", "chr2"))
        self.chrCombox.setItemText(3, _translate("Form", "chr3"))
        self.chrCombox.setItemText(4, _translate("Form", "chr4"))
        self.chrCombox.setItemText(5, _translate("Form", "chr5"))
        self.chrCombox.setItemText(6, _translate("Form", "chr6"))
        self.chrCombox.setItemText(7, _translate("Form", "chr7"))
        self.chrCombox.setItemText(8, _translate("Form", "chr8"))
        self.chrCombox.setItemText(9, _translate("Form", "chr9"))
        self.chrCombox.setItemText(10, _translate("Form", "chr10"))
        self.chrCombox.setItemText(11, _translate("Form", "chr11"))
        self.chrCombox.setItemText(12, _translate("Form", "chr12"))
        self.chrCombox.setItemText(13, _translate("Form", "chr13"))
        self.chrCombox.setItemText(14, _translate("Form", "chr14"))
        self.chrCombox.setItemText(15, _translate("Form", "chr15"))
        self.chrCombox.setItemText(16, _translate("Form", "chr16"))
        self.chrCombox.setItemText(17, _translate("Form", "chr17"))
        self.chrCombox.setItemText(18, _translate("Form", "chr18"))
        self.chrCombox.setItemText(19, _translate("Form", "chr19"))
        self.chrCombox.setItemText(20, _translate("Form", "chr20"))
        self.chrCombox.setItemText(21, _translate("Form", "chr21"))
        self.chrCombox.setItemText(22, _translate("Form", "chr22"))
        self.chrCombox.setItemText(23, _translate("Form", "chr23"))
        self.delButton.setText(_translate("Form", "-"))
        self.tumorFolderButton.setText(_translate("Form", "Tumor Hi-C Sample"))
        self.label.setText(_translate("Form", "Resolution"))
        self.svFileButton.setText(_translate("Form", "Breakpoints File"))
        self.page1SubmitButton.setText(_translate("Form", "Submit -->"))
        # page2
        self.simuaddButton.setText(_translate("Form", "+"))
        self.simuchrCombox.setItemText(0, _translate("Form", "Chromosome List"))
        self.simuchrCombox.setItemText(1, _translate("Form", "chr1"))
        self.simuchrCombox.setItemText(2, _translate("Form", "chr2"))
        self.simuchrCombox.setItemText(3, _translate("Form", "chr3"))
        self.simuchrCombox.setItemText(4, _translate("Form", "chr4"))
        self.simuchrCombox.setItemText(5, _translate("Form", "chr5"))
        self.simuchrCombox.setItemText(6, _translate("Form", "chr6"))
        self.simuchrCombox.setItemText(7, _translate("Form", "chr7"))
        self.simuchrCombox.setItemText(8, _translate("Form", "chr8"))
        self.simuchrCombox.setItemText(9, _translate("Form", "chr9"))
        self.simuchrCombox.setItemText(10, _translate("Form", "chr10"))
        self.simuchrCombox.setItemText(11, _translate("Form", "chr11"))
        self.simuchrCombox.setItemText(12, _translate("Form", "chr12"))
        self.simuchrCombox.setItemText(13, _translate("Form", "chr13"))
        self.simuchrCombox.setItemText(14, _translate("Form", "chr14"))
        self.simuchrCombox.setItemText(15, _translate("Form", "chr15"))
        self.simuchrCombox.setItemText(16, _translate("Form", "chr16"))
        self.simuchrCombox.setItemText(17, _translate("Form", "chr17"))
        self.simuchrCombox.setItemText(18, _translate("Form", "chr18"))
        self.simuchrCombox.setItemText(19, _translate("Form", "chr19"))
        self.simuchrCombox.setItemText(20, _translate("Form", "chr20"))
        self.simuchrCombox.setItemText(21, _translate("Form", "chr21"))
        self.simuchrCombox.setItemText(22, _translate("Form", "chr22"))
        self.simuchrCombox.setItemText(23, _translate("Form", "chr23"))
        self.simudelButton.setText(_translate("Form", "-"))
        self.SimutumorFolderButton.setText(_translate("Form", "Tumor Hi-C Sample"))
        self.Simulabel.setText(_translate("Form", "Resolution"))
        self.assemblyFileButton.setText(_translate("Form", "Assembly CGRs File"))
        self.page2SubmitButton.setText(_translate("Form", "Submit-->"))
        self.useButton.setText(_translate("Form", "Use Local File"))
        self.runButton.setText(_translate("Form", "Run Fit-Hi-C"))
        self.conFolderButton.setText(_translate("Form", "Control Hi-C Sample"))
        self.dis2cntButton.setText(_translate("Form", "Dis2cnt File"))
        # page3
        self.reAssemblyFileButton.setText(_translate("Form", "The assembly CGR File"))
        self.iceButton.setText(_translate("Form", "ICE normalization"))
        self.cnvButton.setText(_translate("Form", "CNV normalization"))
        self.genomelabel.setText(_translate("Form", "Genome"))
        self.enzymelabel.setText(_translate("Form", "Enzyme"))
        self.page3SubmitButton.setText(_translate("Form", "Submit-->"))
        self.exportResult.setText(_translate("Form", "Export result"))

    def open_user_manual(self):
        QDesktopServices.openUrl(QUrl("https://github.com/GaoLabXDU/DeCGR"))

    def download_test_data(self):
        QDesktopServices.openUrl(QUrl("https://github.com/GaoLabXDU/DeCGR/tree/main/TestData"))


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = Interface()
    w.show()
    sys.exit(app.exec_())
