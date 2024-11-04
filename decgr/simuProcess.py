#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from collections import  defaultdict
from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QPushButton, QVBoxLayout, QHBoxLayout, QMessageBox, QWidget, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import bisect
import matplotlib.patches as patches
from matplotlib import gridspec

class MyFigureCanvas(FigureCanvas):
    def __init__(self, width, height, dpi, matrix, interval, assembly_dict, chrom_list, dict_start, dict_end, resolution, output_widget):
        self.fig = plt.figure(figsize=(width, height), dpi=dpi)
        super(MyFigureCanvas, self).__init__(self.fig)
        self.ax = self.fig.add_subplot(111)
        self.data = matrix
        self.vmax = np.percentile(matrix, 98)
        self.vmin = matrix.min()
        cdict = {
            'red': ((0.0, 0.0, 1.0), (1.0, 1.0, 1.0)),
            'green': ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
            'blue': ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
        }
        self.cmap = LinearSegmentedColormap('Rd_Bl_Rd', cdict, 256)
        self.interval = interval
        self.chrom_list = chrom_list
        self.dict_start = dict_start
        self.resolution = resolution
        self.output_widget = output_widget
        self.image = self.ax.imshow(matrix, cmap=self.cmap, vmin=self.vmin, vmax=self.vmax)
        row = matrix.shape[0]
        self.ax.text(row / 2, row + 30, "Original", 
                horizontalalignment='center', verticalalignment='bottom', 
                fontsize=14, color="black", fontweight="bold")
        self.ax.text(row / 2, -30, "Simulated", 
                horizontalalignment='center', verticalalignment='top', 
                fontsize=14, color="black", fontweight="bold")
        # plot chrom boundary
        for chrom, cur_chrom_start in interval.items():
            info = str(chrom) + ":" + str(dict_start[chrom] * self.resolution // 1000000) + '-' + str(dict_end[chrom] * self.resolution // 1000000) + 'M'
            self.ax.text(cur_chrom_start, -10, f'{info}', horizontalalignment='center', color='black', fontsize=12)
            if cur_chrom_start != 0:
                self.ax.plot([0, row], [cur_chrom_start, cur_chrom_start],
                        color="black", linestyle='--', linewidth=2)
                self.ax.plot([cur_chrom_start, cur_chrom_start], [0, row],
                        color="black", linestyle='--', linewidth=2)
        # plot matrix boundary       
        self.ax.plot([0, row], [0, row], color='black', linewidth=2)
        self.ax.plot([0, 0], [0, row], color='black', linewidth=2)
        self.ax.plot([0, row], [0, 0], color='black', linewidth=2)
        self.ax.plot([row, row], [row, 0], color='black', linewidth=2)
        self.ax.plot([row, 0], [row, row], color='black', linewidth=2)
        # plot fragment boundary
        for num, group_df in assembly_dict.items():
            row_pair = []
            for i in range(len(group_df) - 1):
                row1 = group_df.iloc[i]
                row2 = group_df.iloc[i + 1]
                if row1['start'] <= row2['start']:
                    x1 = row1['start'] + interval[row1['chrom']]
                    x2 = row1['end'] + interval[row1['chrom']]
                    y1 = row2['start'] + interval[row2['chrom']]
                    y2 = row2['end'] + interval[row2['chrom']]
                else:
                    y1 = row1['start'] + interval[row1['chrom']]
                    y2 = row1['end'] + interval[row1['chrom']]
                    x1 = row2['start'] + interval[row2['chrom']]
                    x2 = row2['end'] + interval[row2['chrom']]
                info = str(row1['node']) + '-' + str(row2['node'])
                self.ax.text(x1, y1, f'{info}', fontsize=10, verticalalignment='bottom', horizontalalignment='right', color='black')
                self.ax.plot([y1, x1], [y1, y1], color="black", linestyle='--', linewidth=0.5)
                self.ax.plot([x1, x1], [x1, y2], color="black", linestyle='--', linewidth=0.5)
                self.ax.plot([y2, x1], [y2, y2], color="black", linestyle='--', linewidth=0.5)
                self.ax.plot([x2, x2], [x2, y2], color="black", linestyle='--', linewidth=0.5)

                self.ax.plot([y1, y1], [y1, x1], color="black", linestyle='--', linewidth=0.5)
                self.ax.plot([y2, y2], [y2, x1], color="black", linestyle='--', linewidth=0.5)
                self.ax.plot([x1, y2], [x1, x1], color="black", linestyle='--', linewidth=0.5)
                self.ax.plot([x2, y2], [x2, x2], color="black", linestyle='--', linewidth=0.5)
                self.ax.scatter([x1, x2, x1, x2], [y1, y1, y2, y2], color="black", s=5)
                self.ax.scatter([y1, y1, y2, y2], [x1, x2, x1, x2], color="black", s=5)
        self.ax.axis('off')
        # mouse_connect
        self.mpl_connect('scroll_event', self.on_scroll)
        self.mpl_connect('button_press_event', self.on_press)
        self.mpl_connect('motion_notify_event', self.on_move)
        self.mpl_connect('button_release_event', self.on_release)
        self.mpl_connect('button_press_event', self.on_click)
        # set press
        self.press = None
        self.draw()

    def on_scroll(self, event):
        x_left, x_right = self.ax.get_xlim()
        y_bottom, y_top = self.ax.get_ylim()
        # set zoom factor
        zoom_factor = 1.2

        if event.button == 'up':
            scale_factor = 1 / zoom_factor
        elif event.button == 'down':
            scale_factor = zoom_factor
        else:
            scale_factor = 1

        xdata, ydata = event.xdata, event.ydata
        new_width = (x_right - x_left) * scale_factor
        new_height = (y_top - y_bottom) * scale_factor

        x_left = xdata - (xdata - x_left) * scale_factor
        x_right = xdata + (x_right - xdata) * scale_factor
        y_bottom = ydata - (ydata - y_bottom) * scale_factor
        y_top = ydata + (y_top - ydata) * scale_factor

        self.ax.set_xlim([x_left, x_right])
        self.ax.set_ylim([y_bottom, y_top])
        self.draw()

    def on_press(self, event):
        if event.button == 1:
            if event.inaxes != self.ax:
                return
            self.press = (event.xdata, event.ydata)
            if event.dblclick:
                self.on_click(event) 

    def on_move(self, event):
        if self.press is None:
            return
        if event.xdata is None or event.ydata is None:
            return
        x_left, x_right = self.ax.get_xlim()
        y_bottom, y_top = self.ax.get_ylim()
        dx = event.xdata - self.press[0]
        dy = event.ydata - self.press[1]
        self.ax.set_xlim(x_left - dx, x_right - dx)
        self.ax.set_ylim(y_bottom - dy, y_top - dy)
        self.press = (event.xdata, event.ydata)
        self.draw()

    def on_release(self, event):
        self.press = None

    def on_click(self, event):
        if event.inaxes != self.ax:
            return
        x_click = int(event.xdata)
        y_click = int(event.ydata)
        edges = sorted(self.interval.values())
        chromosomes = sorted(self.interval, key=self.interval.get)
        x_chr_idx = bisect.bisect_right(edges, x_click) - 1
        y_chr_idx = bisect.bisect_right(edges, y_click) - 1
        if x_chr_idx == y_chr_idx and 0 <= x_chr_idx < len(chromosomes):
            chr_clicked = chromosomes[x_chr_idx]
            start_pos = self.dict_start[chr_clicked]
            pixel_value = self.data[y_click, x_click]
            final_position_x = (start_pos + (x_click - self.interval[chr_clicked])) * self.resolution
            final_position_y = (start_pos + (y_click - self.interval[chr_clicked])) * self.resolution
            output_text = str(chr_clicked) + ":" + str(final_position_x) + "-" + str(final_position_y)
        elif 0 <= x_chr_idx < len(chromosomes) and 0 <= y_chr_idx < len(chromosomes):
            chr_x = chromosomes[x_chr_idx]
            chr_y = chromosomes[y_chr_idx]
            pixel_value = self.data[y_click, x_click]
            final_position_x_in_chr_x = (self.dict_start[chr_x] + (x_click - self.interval[chr_x])) * self.resolution
            final_position_y_in_chr_y = (self.dict_start[chr_y] + (y_click - self.interval[chr_y])) * self.resolution
            output_text = str(chr_x) + ":" + str(final_position_x_in_chr_x) + "," + str(chr_y) + ":" + str(final_position_y_in_chr_y)
        else:
            output_text = "Error position!"
        self.output_widget.setText(f"Clicked at ({output_text}), Count: {pixel_value:.3f}")
       

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
        # for j in range(i+1, len(genome)):
        for j in range(len(genome)):
            # new_dist = abs(genome.index(j) - genome.index(i))
            if SV_chrom[i] == SV_chrom[j]:
                order_dist = abs(genome[i] - genome[j])
                pos1 = min(genome[i], genome[j])
                pos2 = max(genome[i], genome[j])
            else:
                order_dist = d2cLen
                pos1 = genome[i]
                pos2 = genome[j]

            new_dist = abs(i - j)
            scale = d2c[min(new_dist, d2cLen)] / d2c[min(order_dist, d2cLen)]
            if i < j:
                chrom_pair = SV_chrom[i] + '_' + SV_chrom[j]
            else:
                chrom_pair = SV_chrom[j] + '_' + SV_chrom[i]
            SV_matrix[chrom_pair][pos1, pos2] += matrix_dict[chrom_pair][pos1, pos2] * scale
    return SV_matrix



def simulation(assembly_dict, dis2cnt, control_dict, chrom_list, res):
    simu_dict = {}
    # initial matrix
    for i in range(len(chrom_list)):
        for j in range(i, len(chrom_list)):
            name = str(chrom_list[i]) + '_' + str(chrom_list[j])
            simu_dict[name] = np.full(np.shape(control_dict[name]), 0.0)
   
    for num, group_df in assembly_dict.items():
        genome = []
        SV_chrom = []
        for index, row in group_df.iterrows():            
            if row['orient'] == '+':
                fragment = list(range(row['start'], row['end']))
                cur_chrom = [row['chrom']] * (row['end'] - row['start'])
                genome.extend(fragment)
                SV_chrom.extend(cur_chrom)
            else:
                fragment = list(range(row['end'], row['start'], -1))
                cur_chrom = [row['chrom']] * (row['end'] - row['start'])
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

    return simu_dict


def combine_matrix(dict_matrix, dict_start, dict_end, chrom_list):
    # init_matrix
    sum_bin = 0
    for i in range(len(chrom_list)):
        sum_bin += (dict_end[chrom_list[i]] - dict_start[chrom_list[i]])
    all_matrix = np.full((sum_bin, sum_bin), 0)
    chrom_1 = 0
    interval = {}
    for i in range(len(chrom_list)):
        start1 = dict_start[chrom_list[i]]
        end1 = dict_end[chrom_list[i]]
        length1 = end1 - start1
        chrom_2 = chrom_1
        interval[chrom_list[i]] = chrom_1
        for j in range(i, len(chrom_list)):
            name = str(chrom_list[i]) + '_' + str(chrom_list[j])
            start2 = dict_start[chrom_list[j]]
            end2 = dict_end[chrom_list[j]]
            length2 = end2 - start2
            if j == i:
                cur_matrix = dict_matrix[name]
                cur_matrix = np.triu(cur_matrix)
                cur_matrix += cur_matrix.T - np.diag(cur_matrix.diagonal())
                dict_matrix[name] = cur_matrix
            all_matrix[chrom_1:chrom_1+length1, chrom_2:chrom_2+length2] = dict_matrix[name][start1:end1, start2:end2]
            chrom_2 += length2
        chrom_1 += length1
    return all_matrix, interval


def get_cut_pos(assembly_dict, resolution, dict_count):
    dict_start = {}
    dict_end = {}
    for num, group_df in assembly_dict.items():
        for index, row in group_df.iterrows():
            cur_start = max(0, row['start']-5000000//resolution)
            cur_end = min(row['end'] + 5000000 // resolution, dict_count[row['chrom']])
            if row['chrom'] in dict_start and dict_start[row['chrom']] is not None:
                dict_start[row['chrom']] = min(cur_start, dict_start[row['chrom']])
                dict_end[row['chrom']] = max(cur_end, dict_end[row['chrom']])
            else:
                dict_start[row['chrom']] = cur_start
                dict_end[row['chrom']] = cur_end
    return dict_start, dict_end


def plot_simu_matrix(dict_matrix, simu_dict, dict_start, dict_end, chrom_list, assembly_dict, resolution, output_widget):
    for num, group_df in assembly_dict.items():
        for index, row in group_df.iterrows():
            cur_chrom = row['chrom']
            assembly_dict[num].loc[index, 'start'] -= dict_start[cur_chrom]
            assembly_dict[num].loc[index, 'end'] -= dict_start[cur_chrom]
    original_matrix, interval = combine_matrix(dict_matrix, dict_start, dict_end, chrom_list)
    simu_matrix, interval = combine_matrix(simu_dict, dict_start, dict_end, chrom_list)
    upper = np.triu(simu_matrix)
    upper_indices = np.triu_indices(simu_matrix.shape[0])
    original_matrix[upper_indices] = upper[upper_indices]
    np.fill_diagonal(original_matrix, 0)
    canvas = MyFigureCanvas(width=5, height=4, dpi=100, matrix=original_matrix, interval=interval, resolution=resolution,
                            assembly_dict=assembly_dict, chrom_list=chrom_list, dict_start=dict_start, dict_end=dict_end, output_widget=output_widget)
    return canvas


def read_assembly_result(filaname, resolution):
    df = pd.read_csv(filaname, sep='\t')
    df[['start', 'end']] = df[['start', 'end']] // resolution
    assembly_dict = {num: group.drop(columns=['num']).reset_index(drop=True) for num, group in df.groupby('num')}
    return assembly_dict


