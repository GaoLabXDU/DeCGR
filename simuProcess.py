#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import itertools, operator, os
from iced import normalization
from collections import  defaultdict


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


    return simu_dict


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


def matrix_plot(matrix, interval, figure_name, text):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    cdict = {
        'red': ((0.0, 0.0, 1.0), (1.0, 1.0, 1.0)),
        'green': ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
        'blue': ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
    }
    cmap = LinearSegmentedColormap('Rd_Bl_Rd', cdict, 256)
    M = matrix
    n = M.shape[0]
    vmax = np.percentile(M, 98)
    vmin = M.min()
    # Create the rotation matrix
    t = np.array([[1, 0.5], [-1, 0.5]])
    A = np.dot(np.array([(i[1], i[0]) for i in itertools.product(range(n, -1, -1), range(0, n + 1, 1))]), t)

    # Plot the Heatmap ...
    x = A[:, 1].reshape(n + 1, n + 1)
    y = A[:, 0].reshape(n + 1, n + 1)
    y[y < 0] = -y[y < 0]

    ax.pcolormesh(x, y, np.flipud(M), vmin=vmin, vmax=vmax, cmap=cmap,
                             edgecolor='none', snap=True, linewidth=.01, rasterized=True)
    # plot boundary
    # for i in range(len(interval)):
    #     if interval[i] != 0:
    #         ax.plot([x[interval[i], interval[i]], x[n, interval[i]]], [y[interval[i], interval[i]], y[n, interval[i]]],
    #                  color="black", linestyle='--', linewidth=2)
    #         ax.plot([x[interval[i], interval[i]], x[interval[i], n]], [y[interval[i], interval[i]], y[interval[i], n]],
    #                 color="black", linestyle='--', linewidth=2)
    for i in range(len(interval)):
        if interval[i] != 0:
            # print(x[interval[i], interval[i]], x[n, interval[i]], y[interval[i], interval[i]], y[n, interval[i]])
            ax.plot([interval[i], x[n, interval[i]]], [0, y[n, interval[i]]],
                     color="black", linestyle='--', linewidth=2)
            # print(x[interval[i], interval[i]], x[interval[i], n], y[interval[i], interval[i]], y[interval[i], n])
            ax.plot([interval[i], x[interval[i], n]], [0, y[interval[i], n]],
                    color="black", linestyle='--', linewidth=2)

    ax.plot([0, x[0, 0]], [0, y[0, 0]], color="black", linewidth=2)
    ax.plot([x[0, 0], n], [y[0, 0], 0], color="black", linewidth=2)
    # plt.text(0, n, text, fontsize=16, weight='bold')
    ax.axis('off')

    fig.subplots_adjust(0.15, 0.01, 0.85, 0.90)

    plt.savefig(figure_name)
    # plt.show()


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


def plot_simu_matrix(dict_matrix, dict_start, dict_end, chrom_list, figure_name, simu):
    if simu:
        text = "Simulated Hi-C map"
    else:
        text = "Original Hi-C map"
    chrom_matrix, interval = combine_matrix(dict_matrix, dict_start, dict_end, chrom_list)
    matrix_plot(chrom_matrix, interval, figure_name, text)


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

