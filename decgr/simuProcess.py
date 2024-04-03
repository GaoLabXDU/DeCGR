#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import operator
from iced import normalization
from collections import  defaultdict
from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QPushButton, QVBoxLayout, QHBoxLayout, QMessageBox, QWidget, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas


class MyFigureCanvas(FigureCanvas):
    def __init__(self, width, height, dpi, matrix, interval, rows, cols, parent=None):
        self.fig = plt.figure(figsize=(width, height), dpi=dpi)
        super(MyFigureCanvas, self).__init__(self.fig)
        self.setParent(parent)
        self.matrix = matrix
        self.vmax = np.percentile(matrix, 98)
        self.vmin = matrix.min()
        cdict = {
            'red': ((0.0, 0.0, 1.0), (1.0, 1.0, 1.0)),
            'green': ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
            'blue': ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
        }
        self.cmap = LinearSegmentedColormap('Rd_Bl_Rd', cdict, 256)

        self.axes = self.fig.add_subplot(111)
        self.image = self.axes.imshow(matrix, cmap=self.cmap, vmin=self.vmin, vmax=self.vmax)
        self.axes.axis('off')
        row = matrix.shape[0]
        # plot chrom boundary
        for i in range(len(interval)):
            if interval[i] != 0:
                self.axes.plot([0, row], [interval[i], interval[i]],
                        color="black", linestyle='--', linewidth=2)
                self.axes.plot([interval[i], interval[i]], [0, row],
                        color="black", linestyle='--', linewidth=2)
        # plot fragment boundary
        for i in range(len(rows)):
            if i % 2 == 0:
                print(rows[i], rows[i+1], cols[i], cols[i+1])
                rows[i] += 1
                cols[i] += 1
                if rows[i] < cols[i]:
                    self.axes.plot([rows[i], cols[i+1]], [rows[i], rows[i]], color="black", linestyle='--', linewidth=0.5)
                    self.axes.plot([rows[i+1], cols[i+1]], [rows[i+1], rows[i+1]], color="black", linestyle='--', linewidth=0.5)
                    self.axes.plot([cols[i], cols[i]], [rows[i], cols[i]], color="black", linestyle='--', linewidth=0.5)
                    self.axes.plot([cols[i+1], cols[i+1]], [rows[i], cols[i+1]], color="black", linestyle='--', linewidth=0.5)

                    self.axes.plot([rows[i], rows[i]], [rows[i], cols[i+1]], color="black", linestyle='--', linewidth=0.5)
                    self.axes.plot([rows[i+1], rows[i+1]], [rows[i+1], cols[i+1]], color="black", linestyle='--', linewidth=0.5)
                    self.axes.plot([rows[i], cols[i]], [cols[i], cols[i]], color="black", linestyle='--', linewidth=0.5)
                    self.axes.plot([rows[i], cols[i+1]], [cols[i+1], cols[i+1]], color="black", linestyle='--', linewidth=0.5)

                else:
                    self.axes.plot([cols[i], rows[i+1]], [cols[i], cols[i]], color="black", linestyle='--', linewidth=0.5)
                    self.axes.plot([cols[i+1], rows[i+1]], [cols[i+1], cols[i+1]], color="black", linestyle='--', linewidth=0.5)
                    self.axes.plot([rows[i], rows[i]], [cols[i], rows[i]], color="black", linestyle='--', linewidth=0.5)
                    self.axes.plot([rows[i + 1], rows[i + 1]], [cols[i], rows[i+1]], color="black", linestyle='--', linewidth=0.5)
                    self.axes.plot([cols[i], cols[i]], [cols[i], rows[i+1]], color="black", linestyle='--',
                                   linewidth=0.5)
                    self.axes.plot([cols[i+1], cols[i+1]], [cols[i+1], rows[i+1]], color="black",
                                   linestyle='--', linewidth=0.5)
                    self.axes.plot([cols[i], rows[i]], [rows[i], rows[i]], color="black", linestyle='--', linewidth=0.5)
                    self.axes.plot([cols[i], rows[i+1]], [rows[i + 1], rows[i + 1]], color="black", linestyle='--',
                                   linewidth=0.5)

        self.axes.plot([0, row], [0, row], color='black', linewidth=2)
        self.axes.plot([0, 0], [0, row], color='black', linewidth=2)
        self.axes.plot([0, row], [0, 0], color='black', linewidth=2)
        self.axes.plot([row, row], [row, 0], color='black', linewidth=2)
        self.axes.plot([row, 0], [row, row], color='black', linewidth=2)
        self.zoom_factor = 1.1
        self.move_timer = QTimer(self)
        self.move_timer.timeout.connect(self.move_heatmap)
        self.info_label = None
        self.cid = self.figure.canvas.mpl_connect('button_press_event', self.on_click)
        self.add_move_buttons()
        self.isshowing=False
    # scaling
    def wheelEvent(self, event):
        if event.angleDelta().y() > 0:
            self.axes.set_xlim(self.axes.get_xlim()[0] * self.zoom_factor, self.axes.get_xlim()[1] * self.zoom_factor)
            self.axes.set_ylim(self.axes.get_ylim()[0] * self.zoom_factor, self.axes.get_ylim()[1] * self.zoom_factor)
        else:
            self.axes.set_xlim(self.axes.get_xlim()[0] / self.zoom_factor, self.axes.get_xlim()[1] / self.zoom_factor)
            self.axes.set_ylim(self.axes.get_ylim()[0] / self.zoom_factor, self.axes.get_ylim()[1] / self.zoom_factor)
        self.draw()
        event.accept()
    # move

    def move_heatmap(self):
        direction = self.move_direction
        if direction == "up":
            self.axes.set_ylim(self.axes.get_ylim()[0] + 5, self.axes.get_ylim()[1] + 5)
        elif direction == "down":
            self.axes.set_ylim(self.axes.get_ylim()[0] - 5, self.axes.get_ylim()[1] - 5)
        elif direction == "left":
            self.axes.set_xlim(self.axes.get_xlim()[0] - 5, self.axes.get_xlim()[1] - 5)
        elif direction == "right":
            self.axes.set_xlim(self.axes.get_xlim()[0] + 5, self.axes.get_xlim()[1] + 5)
        self.draw()

    def start_move(self, direction):
        self.move_direction = direction
        self.move_timer.start(100)

    def stop_move(self):
        self.move_timer.stop()

    def on_click(self, event):
        if event.inaxes == self.axes:
            x, y = int(event.xdata), int(event.ydata)
            value = self.matrix[y, x]
            QMessageBox.information(self, "Clicked Position", f"Position: ({x}, {y})\nValue: {value}")

    def add_move_buttons(self):
        outer_layout = QVBoxLayout()
        layout = QHBoxLayout()
        button_up = QPushButton('↑')
        button_up.setStyleSheet("font: bold; font-size: 20px")
        button_up.setFixedSize(60, 40)
        button_up.move(400, 10)
        button_up.pressed.connect(lambda: self.start_move("down"))
        button_up.released.connect(self.stop_move)
        layout.addWidget(button_up)
        button_down = QPushButton('↓')
        button_down.setStyleSheet("font: bold; font-size: 20px;")
        button_down.setFixedSize(60, 40)
        button_down.pressed.connect(lambda: self.start_move("up"))
        button_down.released.connect(self.stop_move)
        layout.addWidget(button_down)

        button_left = QPushButton('←')
        button_left.setStyleSheet("font: bold; font-size: 20px;")
        button_left.setFixedSize(60, 40)
        button_left.pressed.connect(lambda: self.start_move("left"))
        button_left.released.connect(self.stop_move)
        layout.addWidget(button_left)

        button_right = QPushButton('→')
        button_right.setStyleSheet("font: bold; font-size: 20px;")
        button_right.pressed.connect(lambda: self.start_move("right"))
        button_right.setFixedSize(60, 40)
        button_right.released.connect(self.stop_move)
        layout.addWidget(button_right)

        layout.setContentsMargins(200, 600, 200, 0)
        layout.setSpacing(20)
        outer_layout.addLayout(layout)
        self.setLayout(outer_layout)


def simuCount(matrix_dict, dis2cnt, res, genome, SV_chrom):
    d2c = {}
    with open(dis2cnt, 'r') as f:
        for line in f:
            dat = line.split()
            if dat[0][0] == '-':
                d2c['inter'] = float(dat[1])
            else:
                d2c[int(int(dat[0]) / res)] = float(dat[1])
    d2cLen = len(d2c) - 2
    SV_matrix = {}
    for key in matrix_dict:
        row, col = np.shape(matrix_dict[key])
        SV_matrix[key] = np.full((row, col), 0.0)
    for i in range(len(genome)):
        for j in range(i+1, len(genome)):
            # new_dist = abs(genome.index(j) - genome.index(i))
            if SV_chrom[i] == SV_chrom[j]:
                order_dist = abs(genome[i] - genome[j])
            else:
                order_dist = d2cLen
            new_dist = abs(i - j)
            scale = d2c[min(new_dist, d2cLen)] / d2c[min(order_dist, d2cLen)]
            pos1 = min(genome[i], genome[j])
            pos2 = max(genome[i], genome[j])
            chrom_pair = SV_chrom[i] + '_' + SV_chrom[j]
            SV_matrix[chrom_pair][pos1, pos2] += matrix_dict[chrom_pair][pos1, pos2] * scale
    return SV_matrix


def adjust_fragment(chrom, start, end, dict_count, all_path, all_orient):
    for i in range(len(all_path)):
        for j in range(len(all_path[i])):
            index = ord(all_path[i][j]) - 65
            if j == 0:
                if all_orient[i][j] == '+':
                    start[index] = 0
                else:
                    end[index] = dict_count[chrom[index]]-1
            elif j == len(all_path[i]) - 1:
                if all_orient[i][j] == '+':
                    end[index] = dict_count[chrom[index]]-1
                else:
                    start[index] = 0
    return chrom, start, end


def simulation(all_path, all_orient, chrom, start, end, dis2cnt, control_dict, chrom_list, res):
    simu_dict = {}
    rows_dict = defaultdict(list)
    cols_dict = defaultdict(list)
    for i in range(len(chrom_list)):
        for j in range(i, len(chrom_list)):
            name = str(chrom_list[i]) + '_' + str(chrom_list[j])
            simu_dict[name] = np.full(np.shape(control_dict[name]), 0.0)
    for i in range(len(all_path)):
        path = all_path[i]
        orient = all_orient[i]
        genome = []
        SV_chrom = []
        for j in range(len(path)):
            index = ord(path[j]) - 65
            if j < len(path)-1:
                next_index = ord(path[j+1]) - 65
                rows_dict[chrom[index]].append(start[index]-1)
                rows_dict[chrom[index]].append(end[index]-1)
                cols_dict[chrom[next_index]].append(start[next_index]-1)
                cols_dict[chrom[next_index]].append(end[next_index]-1)

            if orient[j] == '+':
                fragment = list(range(start[index], end[index]))
                cur_chrom = [chrom[index]] * (end[index] - start[index])
                genome.extend(fragment)
                SV_chrom.extend(cur_chrom)
            else:
                fragment = list(range(end[index], start[index], -1))
                cur_chrom = [chrom[index]] * (end[index] - start[index])
                SV_chrom.extend(cur_chrom)
                genome.extend(fragment)
            cur_simu_dict = simuCount(control_dict, dis2cnt, res, genome, SV_chrom)
            for key in cur_simu_dict:
                simu_dict[key] += cur_simu_dict[key]

    for key in simu_dict:
        simu_matrix = simu_dict[key]
        for i in range(len(simu_matrix[0])):
            for j in range(len(simu_matrix[0])):
                if simu_matrix[i, j] == 0:
                    simu_dict[key][i, j] = control_dict[key][i, j]

    return simu_dict, rows_dict, cols_dict


def combine_matrix(dict_matrix, dict_start, dict_end, chrom_list):
    # init_matrix
    sum_bin = 0
    for i in range(len(chrom_list)):
        sum_bin += (dict_end[chrom_list[i]] - dict_start[chrom_list[i]])
    all_matrix = np.full((sum_bin, sum_bin), 0)
    chrom_1 = 0
    interval = []
    for i in range(len(chrom_list)):
        start1 = dict_start[chrom_list[i]]
        end1 = dict_end[chrom_list[i]]
        length1 = end1 - start1
        chrom_2 = chrom_1
        interval.append(chrom_1)
        for j in range(i, len(chrom_list)):
            name = str(chrom_list[i]) + '_' + str(chrom_list[j])
            start2 = dict_start[chrom_list[j]]
            end2 = dict_end[chrom_list[j]]
            length2 = end2 - start2
            if j == i:
                cur_matrix = dict_matrix[name]
                cur_matrix = np.triu(cur_matrix)
                cur_matrix += cur_matrix.T - np.diag(cur_matrix.diagonal())
                dict_matrix[name] = normalization.ICE_normalization(cur_matrix)
            all_matrix[chrom_1:chrom_1+length1, chrom_2:chrom_2+length2] = dict_matrix[name][start1:end1, start2:end2]
            chrom_2 += length2
        chrom_1 += length1
    return all_matrix, interval


def get_cut_pos(chrom_list, chrom, start, end, resolution, dict_count):
    dict_start = {}
    dict_end = {}
    for i in range(len(chrom_list)):
        cur_pos = []
        for j in range(len(chrom)):
            if chrom[j] == chrom_list[i]:
                cur_pos.append(start[j])
                cur_pos.append(end[j])
        dict_start[chrom_list[i]] = max(0, min(cur_pos)-5000000//resolution)
        # print(max(0, min(cur_pos)-5000000//resolution))
        dict_end[chrom_list[i]] = min(max(cur_pos) + 5000000 // resolution, dict_count[chrom_list[i]])
        # print(min(max(cur_pos) + 5000000 // resolution, dict_count[chrom_list[i]]))
    return dict_start, dict_end


def plot_simu_matrix(dict_matrix, simu_dict, dict_start, dict_end, chrom_list, rows_dict, cols_dict):
    # rows and cols
    rows_pos = []
    cols_pos = []
    for cur_chrom in chrom_list:
        cur_rows = rows_dict[cur_chrom]
        cur_cols = cols_dict[cur_chrom]
        cur_start = dict_start[cur_chrom]
        for i in range(len(cur_rows)):
            rows_pos.append(cur_rows[i] - cur_start)
            cols_pos.append(cur_cols[i] - cur_start)
    print(rows_pos, cols_pos)
    original_matrix, interval = combine_matrix(dict_matrix, dict_start, dict_end, chrom_list)
    simu_matrix, interval = combine_matrix(simu_dict, dict_start, dict_end, chrom_list)
    # simu_matrix = (simu_matrix - np.min(simu_matrix)) / (np.max(simu_matrix) - np.min(simu_matrix))
    # original_matrix = (original_matrix - np.min(original_matrix)) / (np.max(original_matrix) - np.min(original_matrix))
    upper = np.triu(simu_matrix)
    upper_indices = np.triu_indices(simu_matrix.shape[0])
    original_matrix[upper_indices] = upper[upper_indices]
    np.fill_diagonal(original_matrix, 0)
    canvas = MyFigureCanvas(width=5, height=4, dpi=100, matrix=original_matrix, interval=interval, rows=rows_pos, cols=cols_pos)
    return canvas


def read_assembly_result(filaname, resolution):
    path = {}
    orient = {}
    chrom = {}
    start = {}
    end = {}
    with open(filaname) as f:
        f.readline()
        for line in f:
            line_data = line.strip().split('\t')
            index = ord(line_data[4]) - 65
            chrom[index] = line_data[1]
            start[index] = int(line_data[2]) // resolution
            end[index] = int(line_data[3]) // resolution
            path.setdefault(int(line_data[0]), []).append(line_data[4])
            orient.setdefault(int(line_data[0]), []).append(line_data[5])
    sort_chrom = dict(sorted(chrom.items(), key=operator.itemgetter(0)))
    sort_start = dict(sorted(start.items(), key=operator.itemgetter(0)))
    sort_end = dict(sorted(end.items(), key=operator.itemgetter(0)))
    return list(path.values()), list(orient.values()), list(sort_chrom.values()), list(sort_start.values()), list(
        sort_end.values())

