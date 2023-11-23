import sys
import os

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication
import numpy as np
from PIL import Image
from iced import normalization

from DetermineFragment import FragmentAssembly
from simuProcess import simulation, get_cut_pos, plot_simu_matrix
from plot_reconstruct_map import global_expect, joint_SV, adjust_matrix, plot_assembly_matrix, read_one_assembly
from utils import *



class Interface(QtWidgets.QWidget):

    def __init__(self, parent=None):
        super(Interface, self).__init__(parent)
        self.setupUi(self)

        self.path = os.path.split(os.path.abspath(__file__))[0]
        if not os.path.exists(self.path + '/result'):
            os.mkdir(self.path + '/result')
        if not os.path.exists(self.path + '/tmp'):
            os.mkdir(self.path + '/tmp')

        self.setWindowTitle('DeCGRs')
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

        self.displayWidget.setCurrentIndex(0)
        self.inputWidget.setCurrentIndex(0)

        # self.verticalLayout.setStretch(0, 7)
        # self.verticalLayout.setStretch(1, 3)

        self.gridLayout.setColumnStretch(0, 2)
        self.gridLayout.setColumnStretch(1, 1)
        self.gridLayout.setColumnStretch(2, 1)
        self.gridLayout_2.setColumnStretch(0, 1)
        self.gridLayout_2.setColumnStretch(1, 4)
        self.gridLayout_2.setRowStretch(0, 10)
        self.gridLayout_2.setRowStretch(1, 1)

        self.setStyleSheet('QAbstractButton {font-size: 18px;} .QLabel {font-size: 18px;} \
             .QLineEdit {font-size: 16px;} #toLabel {font-size: 48px;} .QComboBox {font-size: 16px;} \
             .QComboBox:disabled {color: rgb(0,0,0)} .QListWidget {font-size: 18px;}')

        self.displayWidget.setFrameShape(QtWidgets.QFrame.Box)

        self.run_1 = False
        self.run_2 = False
        self.run_3 = False

        self.grp = QtWidgets.QButtonGroup(self)
        self.grp.addButton(self.page1Button, 0)
        self.grp.addButton(self.page2Button, 1)
        self.grp.addButton(self.page3Button, 2)
        self.grp.idClicked.connect(self.changePage)

        self.tumorFolderButton.clicked.connect(lambda: self.inputFilePath('tumor'))
        self.svFileButton.clicked.connect(lambda: self.inputFilePath('sv'))
        re = QtCore.QRegExp('^[1-9][0-9]*$')
        validator = QtGui.QRegExpValidator(re)
        self.resolution.setValidator(validator)
        self.addButton.clicked.connect(self.addItem)
        self.delButton.clicked.connect(self.delItem)
        self.page1SubmitButton.clicked.connect(self.step1)

        self.conFolderButton.clicked.connect(lambda: self.inputFilePath('control'))
        self.assemblyFileButton.clicked.connect(lambda: self.inputFilePath('assembly'))
        self.useButton.toggled.connect(self.fithic)
        self.dis2cntButton.clicked.connect(lambda: self.inputFilePath('dis2cnt'))
        self.page2SubmitButton.clicked.connect(self.step2)

        self.reAssemblyFileButton.clicked.connect(lambda: self.inputFilePath('re'))
        self.page3SubmitButton.clicked.connect(self.step3)

    def changePage(self, id):
        self.displayWidget.setCurrentIndex(id)
        self.inputWidget.setCurrentIndex(id)

        if id == 0 and self.run_1:
            self.plot1()
        elif id == 1 and self.run_2:
            self.plot2()
        elif id == 2 and self.run_3:
            self.plot3()

    def inputFilePath(self, tag):
        if tag == 'tumor':
            path = QtWidgets.QFileDialog.getExistingDirectory(self, 'Input Directory Path', self.path)
            if path == '':
                return
            self.tumorFolder.setText(path)
        elif tag == 'sv':
            path = QtWidgets.QFileDialog.getOpenFileName(self, 'Input File Path', self.path)[0]
            if path == '':
                return
            self.svFile.setText(path)
        elif tag == 'control':
            path = QtWidgets.QFileDialog.getExistingDirectory(self, 'Input Directory Path', self.path)
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

    def fithic(self, checked):
        if checked:
            self.dis2cntButton.setEnabled(True)
            self.dis2cntFile.setEnabled(True)
        else:
            self.dis2cntButton.setEnabled(False)
            self.dis2cntFile.setEnabled(False)

    def step1(self):
        SV_file = self.svFile.text()
        self.res = int(self.resolution.text())
        hic_path = self.tumorFolder.text()
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
                cur_matrix = get_MatrixFile(hic_path, name)
                self.dict_count[self.chrom_list[i]] = np.shape(cur_matrix)[0]
                self.dict_matrix[name] = cur_matrix
                if i == j:
                    normed = normalization.ICE_normalization(cur_matrix)
                    self.dict_norm[name] = normed
                else:
                    self.dict_norm[name] = cur_matrix

        all_path, all_orient, chrom, start, end = FragmentAssembly(SV_file, self.res, self.dict_norm, 
                                                                   self.chrom_list, self.dict_count, self.path)

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

        self.plot1()

        self.run_1 = True
        self.page2Button.setEnabled(True)

    def plot1(self):
        svImage = QtGui.QPixmap(self.path + '/result/CGR_event.png')
        aspect_ratio = svImage.width() / svImage.height()
        new_width = self.page.width() * 0.4
        # new_height = int(new_width / aspect_ratio)
        new_height = self.page.height() * 0.2
        svImage = svImage.scaled(new_width, new_height, QtCore.Qt.IgnoreAspectRatio)
        self.svEvent = QtWidgets.QLabel(self.page)
        self.svEvent.setGeometry(self.page.width() - new_width, 0, new_width, new_height)
        self.svEvent.setPixmap(svImage)
        self.svEvent.show()

        # pixmap = QtGui.QPixmap(self.path + '/tmp/page1_hicmap.png').scaled(self.hicMap.size())
        # self.hicMap.setScaledContents(True)
        # self.hicMap.setPixmap(pixmap)

        # pixmap = QtGui.QPixmap(self.path + '/result/fragment_link.png')#.scaled(self.width(), self.page.height() * 0.3, QtCore.Qt.IgnoreAspectRatio)
        # self.fragmentLink.setScaledContents(True)
        # self.fragmentLink.setPixmap(pixmap)

        image1=Image.open(self.path + '/tmp/page1_hicmap.png')
        image2=Image.open(self.path + '/result/fragment_link.png')
        image1=image1.resize((int(self.page.width()),int(self.page.height()*0.7)))
        image2=image2.resize((int(self.page.width()),int(self.page.height()*0.3)))
        new=Image.new('RGB',(int(self.page.width()),int(self.page.height())))
        new.paste(image1,(0,0))
        new.paste(image2,(0,image1.height))
        new.save(self.path + '/tmp/page1_all.png')

        pixmap=QtGui.QPixmap(self.path + '/tmp/page1_all.png')
        self.hicMap.setScaledContents(True)
        self.hicMap.setPixmap(pixmap)

    def step2(self):
        result_file = self.assemblyFile.text()
        all_path, all_orient, chrom, start, end = read_assembly_result(result_file, self.res)
        control_path = self.conFolder.text()

        control_dict = {}
        control_len = []
        for i in range(len(self.chrom_list)):
            for j in range(i, len(self.chrom_list)):
                name = str(self.chrom_list[i]) + '_' + str(self.chrom_list[j])
                cur_matrix = get_MatrixFile(control_path, name)
                control_dict[name] = cur_matrix
                if i == j:
                    control_len.append(cur_matrix.shape[0])

        if self.runButton.isChecked():
            dis2cnt = runFitHiC(self.path, control_dict, self.chrom_list, control_len, self.res)
        else:
            dis2cnt = self.dis2cntFile.text()

        simu_dict = simulation(all_path, all_orient, chrom, start, end, dis2cnt, control_dict, self.chrom_list,
                               self.res)
        for key in simu_dict:
            SV_matrix = simu_dict[key]
            filename = self.path + "/result/simu_" + key + '_matrix.txt'
            np.savetxt(filename, SV_matrix, fmt='%.2f', delimiter='\t')

        dict_start, dict_end = get_cut_pos(self.chrom_list, chrom, start, end, self.res, self.dict_count)
    
        plot_simu_matrix(self.dict_norm, dict_start, dict_end, self.chrom_list, self.path + '/tmp/page2_varmap.png', simu=False)
        plot_simu_matrix(simu_dict, dict_start, dict_end, self.chrom_list, self.path + '/tmp/page2_simumap.png', simu=True)

        self.plot2()

        self.run_2 = True
        self.page3Button.setEnabled(True)

    def plot2(self):
        mirror = QtGui.QTransform()
        mirror.scale(1, -1)
        pixmap = QtGui.QPixmap(self.path + '/tmp/page2_simumap.png').scaled(self.simuMap.size()).transformed(mirror)
        self.simuMap.setScaledContents(True)
        self.simuMap.setPixmap(pixmap)

        simuLabel = QtWidgets.QLabel('Simulation Hi-C Map', self.simuMap)
        simuLabel.move(self.simuMap.width() * 0.1, self.simuMap.height() * 0.7)
        font = QtGui.QFont('Arial', 32, 75)
        simuLabel.setFont(font)
        simuLabel.show()

        pixmap = QtGui.QPixmap(self.path + '/tmp/page2_varmap.png').scaled(self.varMap.size())
        self.varMap.setScaledContents(True)
        self.varMap.setPixmap(pixmap)

        varLabel = QtWidgets.QLabel('Tumor Hi-C Map', self.varMap)
        varLabel.move(self.varMap.width() * 0.1, self.varMap.height() * 0.3)
        varLabel.setFont(font)
        varLabel.show()

    def step3(self):
        result_file = self.reAssemblyFile.text()
        path, orient, chrom, start, end = read_one_assembly(result_file, self.res)
        glob_exp = global_expect(self.dict_matrix, self.chrom_list)
        SV_matrix = joint_SV(path, orient, chrom, start, end, self.dict_matrix)
        matrix, interval = adjust_matrix(SV_matrix, glob_exp, path, start, end)
        plot_assembly_matrix(SV_matrix, interval, path, orient, self.path)

        self.plot3()
        self.run_3 = True

    def plot3(self):
        # rotation = QtGui.QTransform().rotate(-45)
        pixmap = QtGui.QPixmap(self.path + '/result/restruct_map.png').scaled(
            self.reconstructMap.height(), self.reconstructMap.height())

        # self.reconstructMap.setAlignment(QtCore.Qt.AlignCenter)
        self.reconstructMap.setScaledContents(True)
        self.reconstructMap.setPixmap(pixmap)

    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(1395, 940)
        self.gridLayout_2 = QtWidgets.QGridLayout(Form)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.widget = QtWidgets.QWidget(Form)
        self.widget.setObjectName("widget")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.widget)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.inputWidget = QtWidgets.QStackedWidget(self.widget)
        self.inputWidget.setObjectName("inputWidget")
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
        self.inputWidget.addWidget(self.inputPage1)
        self.inputPage2 = QtWidgets.QWidget()
        self.inputPage2.setObjectName("inputPage2")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.inputPage2)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.assemblyFileButton = QtWidgets.QPushButton(self.inputPage2)
        self.assemblyFileButton.setObjectName("assemblyFileButton")
        self.gridLayout_3.addWidget(self.assemblyFileButton, 3, 0, 1, 2)
        self.assemblyFile = QtWidgets.QLineEdit(self.inputPage2)
        self.assemblyFile.setEnabled(True)
        self.assemblyFile.setReadOnly(True)
        self.assemblyFile.setObjectName("assemblyFile")
        self.gridLayout_3.addWidget(self.assemblyFile, 2, 0, 1, 2)
        self.page2SubmitButton = QtWidgets.QPushButton(self.inputPage2)
        self.page2SubmitButton.setObjectName("page2SubmitButton")
        self.gridLayout_3.addWidget(self.page2SubmitButton, 7, 1, 1, 1)
        self.useButton = QtWidgets.QRadioButton(self.inputPage2)
        self.useButton.setObjectName("useButton")
        self.gridLayout_3.addWidget(self.useButton, 4, 1, 1, 1)
        self.conFolder = QtWidgets.QLineEdit(self.inputPage2)
        self.conFolder.setEnabled(True)
        self.conFolder.setReadOnly(True)
        self.conFolder.setObjectName("conFolder")
        self.gridLayout_3.addWidget(self.conFolder, 0, 0, 1, 2)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_3.addItem(spacerItem1, 8, 0, 1, 2)
        self.runButton = QtWidgets.QRadioButton(self.inputPage2)
        self.runButton.setChecked(True)
        self.runButton.setObjectName("runButton")
        self.gridLayout_3.addWidget(self.runButton, 4, 0, 1, 1)
        self.conFolderButton = QtWidgets.QPushButton(self.inputPage2)
        self.conFolderButton.setObjectName("conFolderButton")
        self.gridLayout_3.addWidget(self.conFolderButton, 1, 0, 1, 2)
        self.dis2cntButton = QtWidgets.QPushButton(self.inputPage2)
        self.dis2cntButton.setEnabled(False)
        self.dis2cntButton.setObjectName("dis2cntButton")
        self.gridLayout_3.addWidget(self.dis2cntButton, 6, 1, 1, 1)
        self.dis2cntFile = QtWidgets.QLineEdit(self.inputPage2)
        self.dis2cntFile.setEnabled(True)
        self.dis2cntFile.setReadOnly(True)
        self.dis2cntFile.setObjectName("dis2cntFile")
        self.gridLayout_3.addWidget(self.dis2cntFile, 6, 0, 1, 1)
        self.inputWidget.addWidget(self.inputPage2)
        self.inputPage3 = QtWidgets.QWidget()
        self.inputPage3.setObjectName("inputPage3")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.inputPage3)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.reAssemblyFileButton = QtWidgets.QPushButton(self.inputPage3)
        self.reAssemblyFileButton.setObjectName("reAssemblyFileButton")
        self.gridLayout_4.addWidget(self.reAssemblyFileButton, 1, 0, 1, 2)
        self.label_2 = QtWidgets.QLabel(self.inputPage3)
        self.label_2.setText("")
        self.label_2.setObjectName("label_2")
        self.gridLayout_4.addWidget(self.label_2, 2, 0, 1, 1)
        self.reAssemblyFile = QtWidgets.QLineEdit(self.inputPage3)
        self.reAssemblyFile.setEnabled(True)
        self.reAssemblyFile.setReadOnly(True)
        self.reAssemblyFile.setObjectName("reAssemblyFile")
        self.gridLayout_4.addWidget(self.reAssemblyFile, 0, 0, 1, 2)
        self.page3SubmitButton = QtWidgets.QPushButton(self.inputPage3)
        self.page3SubmitButton.setObjectName("page3SubmitButton")
        self.gridLayout_4.addWidget(self.page3SubmitButton, 2, 1, 1, 1)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_4.addItem(spacerItem2, 3, 0, 1, 2)
        self.inputWidget.addWidget(self.inputPage3)
        self.verticalLayout_3.addWidget(self.inputWidget)
        self.gridLayout_2.addWidget(self.widget, 0, 0, 2, 1)
        self.widget_2 = QtWidgets.QWidget(Form)
        self.widget_2.setObjectName("widget_2")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.widget_2)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.displayWidget = QtWidgets.QStackedWidget(self.widget_2)
        self.displayWidget.setObjectName("displayWidget")
        self.page = QtWidgets.QWidget()
        self.page.setObjectName("page")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.page)
        self.verticalLayout.setObjectName("verticalLayout")
        self.hicMap = QtWidgets.QLabel(self.page)
        self.hicMap.setText("")
        self.hicMap.setObjectName("hicMap")
        self.verticalLayout.addWidget(self.hicMap)
        self.displayWidget.addWidget(self.page)
        self.page_2 = QtWidgets.QWidget()
        self.page_2.setObjectName("page_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.page_2)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.varMap = QtWidgets.QLabel(self.page_2)
        self.varMap.setText("")
        self.varMap.setObjectName("varMap")
        self.verticalLayout_2.addWidget(self.varMap)
        self.simuMap = QtWidgets.QLabel(self.page_2)
        self.simuMap.setText("")
        self.simuMap.setObjectName("simuMap")
        self.verticalLayout_2.addWidget(self.simuMap)
        self.displayWidget.addWidget(self.page_2)
        self.page_3 = QtWidgets.QWidget()
        self.page_3.setObjectName("page_3")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.page_3)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.reconstructMap = QtWidgets.QLabel(self.page_3)
        self.reconstructMap.setText("")
        self.reconstructMap.setObjectName("reconstructMap")
        self.horizontalLayout_3.addWidget(self.reconstructMap)
        self.displayWidget.addWidget(self.page_3)
        self.verticalLayout_4.addWidget(self.displayWidget)
        self.gridLayout_2.addWidget(self.widget_2, 0, 1, 1, 1)
        self.functionWidget = QtWidgets.QWidget(Form)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.functionWidget.sizePolicy().hasHeightForWidth())
        self.functionWidget.setSizePolicy(sizePolicy)
        self.functionWidget.setObjectName("functionWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.functionWidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem3)
        self.page1Button = QtWidgets.QPushButton(self.functionWidget)
        self.page1Button.setObjectName("page1Button")
        self.horizontalLayout.addWidget(self.page1Button)
        spacerItem4 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem4)
        self.page2Button = QtWidgets.QPushButton(self.functionWidget)
        self.page2Button.setEnabled(False)
        self.page2Button.setObjectName("page2Button")
        self.horizontalLayout.addWidget(self.page2Button)
        spacerItem5 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem5)
        self.page3Button = QtWidgets.QPushButton(self.functionWidget)
        self.page3Button.setEnabled(False)
        self.page3Button.setObjectName("page3Button")
        self.horizontalLayout.addWidget(self.page3Button)
        spacerItem6 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem6)
        self.gridLayout_2.addWidget(self.functionWidget, 1, 1, 1, 1)

        self.retranslateUi(Form)
        self.inputWidget.setCurrentIndex(0)
        self.displayWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
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
        self.tumorFolderButton.setText(_translate("Form", "Tumor Hi-C Folder"))
        self.label.setText(_translate("Form", "Resolution"))
        self.svFileButton.setText(_translate("Form", "Structure Variance File"))
        self.page1SubmitButton.setText(_translate("Form", "Submit -->"))
        self.assemblyFileButton.setText(_translate("Form", "Assembly File"))
        self.page2SubmitButton.setText(_translate("Form", "Submit-->"))
        self.useButton.setText(_translate("Form", "Use Local File"))
        self.runButton.setText(_translate("Form", "Run Fit-Hi-C"))
        self.conFolderButton.setText(_translate("Form", "Control Hi-C Folder"))
        self.dis2cntButton.setText(_translate("Form", "Dis2cnt File"))
        self.reAssemblyFileButton.setText(_translate("Form", "Assembly File"))
        self.page3SubmitButton.setText(_translate("Form", "Submit-->"))
        self.page1Button.setText(_translate("Form", "CGR Assembly"))
        self.page2Button.setText(_translate("Form", "Simulation"))
        self.page3Button.setText(_translate("Form", "Reconstruct Hi-C Map"))

if __name__ == '__main__':

    app = QApplication(sys.argv)
    interface = Interface()
    interface.show()
    sys.exit(app.exec_())
