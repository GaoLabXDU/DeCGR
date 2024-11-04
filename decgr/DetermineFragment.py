#!/usr/bin/env python
from scipy.ndimage import gaussian_filter
from sklearn.decomposition import PCA
import pandas as pd
import random, itertools
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas



def calculate_orient(matrix, bp1, bp2, bpLink, chrom1, chrom2):
    '''
    Calculate the orient
    '''
    # Location peak
    index = np.unravel_index(np.argmax(matrix), np.array(matrix).shape)
    cur_row = np.shape(matrix)[0] // 2
    cur_col = np.shape(matrix)[1] // 2
    if 0 <= index[0] < cur_row and 0 <= index[1] < cur_col:
        row_peak = -1
        col_peak = 1
    elif 0 <= index[0] < cur_row and index[1] >= cur_col:
        row_peak = -1
        col_peak = -1
    elif index[0] >= cur_row and 0 <= index[1] < cur_col:
        row_peak = 1
        col_peak = 1
    else:
        row_peak = 1
        col_peak = -1
    # row and col correlation
    row_pos = np.arange(1, np.shape(matrix)[0] + 1)
    col_pos = np.arange(1, np.shape(matrix)[1] + 1)
    cols = np.mean(matrix, axis=0)
    rows = np.mean(matrix, axis=1)
    rows_coff = np.corrcoef(row_pos, rows)[0, 1]
    cols_coff = -np.corrcoef(col_pos, cols)[0, 1]
    # diagonal correlation
    d = np.diag(matrix)
    diag_rows_coff = np.corrcoef(rows, d)[0, 1]
    diag_cols_coff = np.corrcoef(cols, d)[0, 1]
    sum_row = row_peak + rows_coff + diag_rows_coff
    sum_col = col_peak + cols_coff + diag_cols_coff
    if sum_row > 0 and sum_col > 0:
        bpLink.loc[len(bpLink.index)] = [chrom1, bp1, '+', chrom2, bp2, '+']
    elif sum_row > 0 > sum_col:
        bpLink.loc[len(bpLink.index)] = [chrom1, bp1, '+', chrom2, bp2, '-']
    elif sum_row < 0 < sum_col:
        bpLink.loc[len(bpLink.index)] = [chrom1, bp1, '-', chrom2, bp2, '+']
    else:
        bpLink.loc[len(bpLink.index)] = [chrom1, bp1, '-', chrom2, bp2, '-']
    return bpLink


def calculate_std(matrix):
    count = []
    row, col = np.shape(matrix)
    for i in range(row):
        for j in range(col):
            if matrix[i][j] != 0:
                count.append(matrix[i][j])
    # remove outlier
    count = np.array(count)
    q1 = np.quantile(count, 0.25)
    q3 = np.quantile(count, 0.75)
    iqr = q3 - q1
    upper_bound = q3 + (3 * iqr)
    index = count < upper_bound

    return np.std(count[index])


def read_SV_info(bp_file):
    '''
        read SV breakpoint information
        :param breakpoint file name :
        :return: {'chrom1', 'pos1', 'chrom2', 'pos2'}
    '''
    result = pd.DataFrame(columns=['chrom1', 'pos1', 'chrom2', 'pos2'])
    with open(bp_file) as f:
        for line in f:
            line_data = line.strip().split('\t')
            result.loc[len(result.index)] = [line_data[0], int(line_data[1]), line_data[2], int(line_data[3])]
    return result


def extend_stretches(arr, min_seed_len=2, max_gap=1, min_stripe_len=5):
    arr.sort()
    pieces = np.split(arr, np.where(np.diff(arr) != 1)[0] + 1)
    filtered = [p for p in pieces if len(p) >= min_seed_len]  # can't be empty
    stripes = []
    seed = filtered[0]
    for p in filtered[1:]:
        if p[0] - seed[-1] < (max_gap + 2):
            seed = np.r_[seed, p]
        else:
            if seed[-1] - seed[0] + 1 >= min_stripe_len:
                stripes.append([seed[0], seed[-1] + 1])
            seed = p

    if seed[-1] - seed[0] + 1 >= min_stripe_len:
        stripes.append([seed[0], seed[-1] + 1])
    return stripes


def locate(pca1, left_most=True):

    loci_1 = extend_stretches(np.where(pca1 >= 0)[0])
    loci_2 = extend_stretches(np.where(pca1 <= 0)[0])
    loci = loci_1 + loci_2
    loci.sort()
    if left_most:
        b = loci[0][1] - 1
    else:
        b = loci[-1][0]

    return b


def recalculate_breakpoint(matrix, orient1, orient2):
    '''

    :param SV_region matrix:
    :return: breakpoint position
    '''
    rowmask = matrix.sum(axis=1) != 0
    colmask = matrix.sum(axis=0) != 0
    new = matrix[rowmask][:, colmask]
    # row
    corr = gaussian_filter(np.corrcoef(new, rowvar=True), sigma=2)
    corr[np.isnan(corr)] = 0
    try:
        pca = PCA(n_components=3, whiten=True)
        pca1 = pca.fit_transform(corr)[:, 0]
        if orient1 == '+':
            loc = locate(pca1, left_most=False)
        else:
            loc = locate(pca1, left_most=True)
        up_i = np.where(rowmask)[0][loc]  # included
    except:
        up_i = 0
    # col
    corr = gaussian_filter(np.corrcoef(new, rowvar=False), sigma=2)
    try:
        pca = PCA(n_components=3, whiten=True)
        pca1 = pca.fit_transform(corr)[:, 0]
        if orient2 == '+':
            loc = locate(pca1, left_most=True)
        else:
            loc = locate(pca1, left_most=False)
        down_i = np.where(colmask)[0][loc]
    except:
        down_i = 0
    return up_i, down_i


def get_MatrixFile(hic_path, name):
    matrix = 0
    for filewalks in os.walk(hic_path):
        for files in filewalks[2]:
            if name in files:
                matrix = np.loadtxt(os.path.join(filewalks[0], files))
    return matrix


def get_SV_orient(bp_result, dict_matrix, res, merge):
    # init orient
    bpLink = pd.DataFrame(columns=['chrom1', 'pos1', 'orient1', 'chrom2', 'pos2', 'orient2'])
    overlap = False
    for k in bp_result.index:
        row = bp_result.loc[k]
        cur_chr1 = row['chrom1']
        cur_chr2 = row['chrom2']
        name = str(cur_chr1) + '_' + str(cur_chr2)
        matrix = dict_matrix[name]
        pos1 = row['pos1']
        pos2 = row['pos2']
        if pos1 < merge or pos2 < merge or (pos1 + merge) // res > np.shape(matrix)[0] or (pos2 + merge) // res > np.shape(matrix)[1]:
            merge = min(pos1, pos2, np.shape(matrix)[0]*res-pos1,  np.shape(matrix)[1]*res-pos2)
        pos1_start = (pos1 - merge) // res
        pos1_end = (pos1 + merge) // res
        pos2_start = (pos2 - merge) // res
        pos2_end = (pos2 + merge) // res
        cur_matrix = matrix[pos1_start:pos1_end, pos2_start:pos2_end]
        # instead of zero bin
        for m in range(np.shape(cur_matrix)[0]):
            for n in range(np.shape(cur_matrix)[1]):
                if cur_matrix[m][n] == 0:
                        cur_matrix[m][n] = round(random.uniform(0, 1), 2)

        row, col = np.shape(cur_matrix)
        # calculate std cut-off value
        cut_off = calculate_std(cur_matrix)
        cut_off *= 2
        # calculate orient
        cut_1 = merge // res
        cut_2 = merge // res
        cut_matrix = {}
        cut_std = []
        cur_cut_matrix = cur_matrix[0:cut_1, 0:cut_2]
        cut_std.append(calculate_std(cur_cut_matrix))
        cut_matrix[0] = cur_cut_matrix
        cur_cut_matrix = cur_matrix[cut_1:row, 0:cut_2]
        cut_std.append(calculate_std(cur_cut_matrix))
        cut_matrix[1] = cur_cut_matrix
        cur_cut_matrix = cur_matrix[0:cut_1, cut_2:col]
        cut_std.append(calculate_std(cur_cut_matrix))
        cut_matrix[2] = cur_cut_matrix
        cur_cut_matrix = cur_matrix[cut_1:row, cut_2:col]
        cut_std.append(calculate_std(cur_cut_matrix))
        cut_matrix[3] = cur_cut_matrix
        count = sum(i > cut_off for i in cut_std)
        if count == 2:
            filtered = filter(lambda x: x > cut_off, cut_std)
            pos = list(map(cut_std.index, filtered))
            if pos[0] == 0 and pos[1] == 3:
                bpLink = calculate_orient(cut_matrix[0], pos1 // res, pos2 // res, bpLink, cur_chr1, cur_chr2)
                bpLink = calculate_orient(cut_matrix[3], pos1 // res, pos2 // res, bpLink, cur_chr1, cur_chr2)
            elif pos[0] == 1 and pos[1] == 2:
                bpLink = calculate_orient(cut_matrix[1], pos1 // res, pos2 // res, bpLink, cur_chr1, cur_chr2)
                bpLink = calculate_orient(cut_matrix[2], pos1 // res, pos2 // res, bpLink, cur_chr1, cur_chr2)
            else:
                overlap = True
                index = cut_std.index(max(cut_std))
                bpLink = calculate_orient(cut_matrix[index], pos1 // res, pos2 // res, bpLink, cur_chr1, cur_chr2)

        else:
            index = cut_std.index(max(cut_std))
            bpLink = calculate_orient(cut_matrix[index], pos1 // res, pos2 // res, bpLink, cur_chr1, cur_chr2)

    return bpLink, overlap


def fragment_boundary(chrom, start, end, orient, index, SV_fragment):
    cur_start = start
    cur_end = end
    num = [i for i in range(SV_fragment.shape[0])]
    del num[index]
    for i in num:
        row = SV_fragment.loc[i]
        if row['chrom1'] == chrom:
            fragment1 = [m for m in range(start, end+1)]
            fragment2 = [n for n in range(row['start1'], row['end1']+1)]
            ret1 = list(set(fragment1).intersection(set(fragment2)))
            if len(ret1) > (len(fragment1) * 0.6) and len(ret1) > (len(fragment2) * 0.6):
                if orient == '+' and row['orient1'] == '-':
                    cur_start = start
                    cur_end = row['end1']
                elif orient == '+' and row['orient1'] == '+':
                    cur_start = min(row['start1'], start)
                    cur_end = max(row['end1'], end)
                elif orient == '-' and row['orient1'] == '-':
                    cur_start = min(row['start1'], start)
                    cur_end = max(row['end1'], end)
                else:
                    cur_start = row['start1']
                    cur_end = end

        if row['chrom2'] == chrom:
            fragment1 = [m for m in range(start, end + 1)]
            fragment2 = [n for n in range(row['start2'], row['end2'] + 1)]
            ret1 = list(set(fragment1).intersection(set(fragment2)))
            if len(ret1) > (len(fragment1) * 0.6) and len(ret1) > (len(fragment2) * 0.6):
                # print(chrom, start, end)
                # print(row['chrom1'], row['start1'], row['end1'] + 1)
                if orient == '+' and row['orient2'] == '-':
                    cur_start = start
                    cur_end = row['end2']
                elif orient == '+' and row['orient2'] == '+':
                    cur_start = min(row['start2'], start)
                    cur_end = max(row['end2'], end)
                elif orient == '-' and row['orient2'] == '-':
                    cur_start = min(row['start2'], start)
                    cur_end = max(row['end2'], end)
                else:
                    cur_start = row['start2']
                    cur_end = end

    return cur_start, cur_end


def coverFragment(SV_fragment):
    for i in SV_fragment.index:
        row = SV_fragment.loc[i]
        start1, end1 = fragment_boundary(row['chrom1'], row['start1'], row['end1'], row['orient1'], i, SV_fragment)
        start2, end2 = fragment_boundary(row['chrom2'], row['start2'], row['end2'], row['orient2'], i, SV_fragment)
        SV_fragment.loc[i, 'start1'] = start1
        SV_fragment.loc[i, 'end1'] = end1
        SV_fragment.loc[i, 'start2'] = start2
        SV_fragment.loc[i, 'end2'] = end2
    return SV_fragment


def genomeFragment(SV_fragment):
    chrom = []
    start = []
    end = []
    node_list = []
    for i in SV_fragment.index:
        row = SV_fragment.loc[i]
        start1 = row['start1']
        end1 = row['end1']
        start2 = row['start2']
        end2 = row['end2']
        if start1 not in start and end1 not in end:
            start.append(start1)
            end.append(end1)
            chrom.append(row['chrom1'])
            node_list.append(row['fragment1'])
        if start2 not in start and end2 not in end:
            start.append(start2)
            end.append(end2)
            chrom.append(row['chrom2'])
            node_list.append(row['fragment2'])
        if start1 in start and end1 in end:
            index = start.index(start1)
            if chrom[index] != row['chrom1']:
                start.append(start1)
                end.append(end1)
                chrom.append(row['chrom1'])
                node_list.append(row['fragment1'])
        if start2 in start and end2 in end:
            index = start.index(start2)
            if chrom[index] != row['chrom2']:
                start.append(start2)
                end.append(end2)
                chrom.append(row['chrom2'])
                node_list.append(row['fragment2'])

    return chrom, start, end, node_list

segment_to_node = {}
node_labels = itertools.chain(
    (chr(i) for i in range(65, 91)),  # Single letters A-Z
    (''.join(pair) for pair in itertools.product("ABCDEFGHIJKLMNOPQRSTUVWXYZ", repeat=2))  # Double letters AA, AB, etc.
)
def get_node_name(chrom, start, end):
    segment_key = (chrom, start, end)
    if segment_key not in segment_to_node:
        segment_to_node[segment_key] = next(node_labels)
    return segment_to_node[segment_key]

def constructPath(SV_fragment):
    all_path = []
    all_orient = []
    SV_fragment['flag'] = 0
    SV_fragment.sort_values(["chrom1", "start1"], inplace=True)
    SV_fragment.reset_index(drop=True, inplace=True)
    SV_fragment['fragment1'] = SV_fragment.apply(lambda row: get_node_name(row['chrom1'], row['start1'], row['end1']), axis=1)
    SV_fragment['fragment2'] = SV_fragment.apply(lambda row: get_node_name(row['chrom2'], row['start2'], row['end2']), axis=1)
    # path
    while SV_fragment['flag'].sum() != SV_fragment.shape[0]:
        index = SV_fragment[SV_fragment['flag'] == 0].index[0]
        init_link = SV_fragment.loc[index]
        node = init_link['fragment1']
        next_orient = init_link['orient2']
        next_node = init_link['fragment2']
        path_list = []
        orient_list = []
        path_list.append(node)
        path_list.append(next_node)
        orient_list.append(init_link['orient1'])
        orient_list.append(next_orient)
        SV_fragment.loc[index, 'flag'] = 1
        count = 1
        while count < SV_fragment.shape[0]:
            row = SV_fragment.loc[count]
            if row['flag'] == 0 and row['orient1'] == next_orient and next_node == row['fragment1']:
                next_node = row['fragment2']
                next_orient = row['orient2']
                path_list.append(next_node)
                orient_list.append(next_orient)
                SV_fragment.loc[count, 'flag'] = 1
                count = 1
            elif row['flag'] == 0 and row['orient2'] != next_orient and next_node == row['fragment2'] :
                next_node = row['fragment1']
                if row['orient1'] == '+':
                    next_orient = '-'
                else:
                    next_orient = '+'
                path_list.append(next_node)
                orient_list.append(next_orient)
                SV_fragment.loc[count, 'flag'] = 1
                count = 1
            else:
                count += 1
        all_path.append(path_list)
        all_orient.append(orient_list)
    return all_path, all_orient


def plot_SV_genome(path, all_orient, sv_fragment):
    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot(111)
    h = 0
    length_list = []
    CGR_info = {'num': [], 'chrom': [], 'start': [], 'end': [], 'node': [], 'orient': []}
    for i in range(len(path)):
        cur_start = 0
        count = 0
        colors = list(mcolors.TABLEAU_COLORS.keys())
        for j in range(len(path[i])):
            first_match = sv_fragment[(sv_fragment['fragment1'] == path[i][j]) | (sv_fragment['fragment2'] == path[i][j])].iloc[0]
            if first_match['fragment1'] == path[i][j]:
                fragment_chrom = first_match['chrom1']
                fragment_start = first_match['start1']
                fragment_end = first_match['end1']
            if first_match['fragment2'] == path[i][j]:
                fragment_chrom = first_match['chrom2']
                fragment_start = first_match['start2']
                fragment_end = first_match['end2']
            CGR_info['num'].append(i)
            CGR_info['chrom'].append(fragment_chrom)
            CGR_info['start'].append(fragment_start)
            CGR_info['end'].append(fragment_end)
            CGR_info['node'].append(path[i][j])
            CGR_info['orient'].append(all_orient[i][j])
            length = fragment_end - fragment_start + 1
            count += length
            textPos = cur_start + length // 2
            plt.text(textPos, h+2, path[i][j], fontsize=10, weight='bold')
            if all_orient[i][j] == '+':
                pos1 = cur_start
                pos2 = pos1 + length
            else:
                pos2 = cur_start
                pos1 = pos2 + length
            ax.annotate("",
                        xy=[pos2, h],
                        xytext=[pos1, h],
                        size=4,
                        arrowprops=dict(color=mcolors.TABLEAU_COLORS[colors[i % len(colors)]], headwidth=8, width=4, headlength=4), zorder=0)

            cur_start = cur_start + length + 1
        length_list.append(count)
        h -= 10

    plt.xlim([0, max(length_list)+5])
    plt.ylim([h, 3])
    plt.axis('off')
    canvas = FigureCanvas(fig)
    return canvas, CGR_info



def combine_matrix(dict_matrix, dict_start, dict_end, chrom_list, resolution):
    # init_matrix
    sum_bin = 0
    for i in range(len(chrom_list)):
        sum_bin += (dict_end[chrom_list[i]] - dict_start[chrom_list[i]])
    all_matrix = np.full((sum_bin, sum_bin), 0)
    chrom_1 = 0
    interval = []
    text_pos = []
    text = []
    for i in range(len(chrom_list)):
        start1 = dict_start[chrom_list[i]]
        end1 = dict_end[chrom_list[i]]
        length1 = end1 - start1
        chrom_2 = chrom_1
        interval.append(chrom_1)
        text_pos.append(chrom_1)
        text_pos.append(chrom_1+length1)
        text.append(str(start1*resolution // 1000000) + 'M')
        text.append(str(end1*resolution // 1000000) + 'M')
        for j in range(i, len(chrom_list)):
            name = str(chrom_list[i]) + '_' + str(chrom_list[j])
            start2 = dict_start[chrom_list[j]]
            end2 = dict_end[chrom_list[j]]
            length2 = end2 - start2
            all_matrix[chrom_1:chrom_1+length1, chrom_2:chrom_2+length2] = dict_matrix[name][start1:end1, start2:end2]
            chrom_2 += length2
        chrom_1 += length1
    return all_matrix, interval, text_pos, text


def matrix_plot(matrix, interval, chrom_text_pos, chrom_text,
                SV_fragment, chrom, start, end, node_list, chrom_list, plot_start, plot_end):
    
    fig, (ax_heatmap, ax_genome) = plt.subplots(2, 1, 
                                                gridspec_kw={'height_ratios': [8, 3]}, figsize=(10, 10), sharex=True)
    fig.subplots_adjust(hspace=0.05) 
    cdict = {
        'red': ((0.0, 0.0, 1.0), (1.0, 1.0, 1.0)),
        'green': ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
        'blue': ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
    }

    cmap = LinearSegmentedColormap('Rd_Bl_Rd', cdict, 256)
    M = matrix
    n = M.shape[0]
    vmax = np.percentile(M, 96)
    vmin = M.min()
    # Create the rotation matrix
    t = np.array([[1, 0.5], [-1, 0.5]])
    A = np.dot(np.array([(i[1], i[0]) for i in itertools.product(range(n, -1, -1), range(0, n + 1, 1))]), t)
    # Plot the Heatmap ...
    x = A[:, 1].reshape(n + 1, n + 1)
    y = A[:, 0].reshape(n + 1, n + 1)
    y[y < 0] = -y[y < 0]
    ax_heatmap.pcolormesh(x, y, np.flipud(M), vmin=vmin, vmax=vmax, cmap=cmap,
                  edgecolor='none', snap=True, linewidth=.01, rasterized=True)
    # plot boundary
    for i in range(len(interval)):
        if interval[i] != 0:
            ax_heatmap.plot([interval[i], x[n, interval[i]]], [0, y[n, interval[i]]],
                     color="black", linestyle='--', linewidth=2)
            ax_heatmap.plot([interval[i], x[(n-interval[i]), n]], [0, y[(n-interval[i]), n]],
                     color="black", linestyle='--', linewidth=2)
            
    ax_heatmap.plot([0, x[0, 0]], [0, y[0, 0]], color="black", linewidth=2)
    ax_heatmap.plot([x[0, 0], n], [y[0, 0], 0], color="black", linewidth=2)
    ax_heatmap.axis('off')
    ## plot genome
    color_list = []
    chrom_start = 0
    dict_start = {}
    for i in range(len(chrom_list)):
        count = 0
        dict_start[chrom_list[i]] = chrom_start
        chrom_length = plot_end[chrom_list[i]]-plot_start[chrom_list[i]]
        ax_genome.add_patch(patches.Rectangle((chrom_start, 0), width=chrom_length, height=0.8, color='gray', zorder=1, alpha=.4, linewidth=1, edgecolor='black'))
        rect_mid = chrom_start + chrom_length / 2
        ax_genome.text(rect_mid, 0.8, chrom_list[i], ha='center', va='center', fontsize=10, color='black')
        ax_genome.text(chrom_start, 0.8, f"{chrom_text[count]}", ha='right', va='center', fontsize=8, color='black')
        ax_genome.text(chrom_start + chrom_length, 0.8, f"{chrom_text[count+1]}", ha='left', va='center', fontsize=8, color='black')
        chrom_start += chrom_length
        count += 1
    ## plot fragment
    for i in range(len(start)):
        cur_start = dict_start[chrom[i]]
        rect_start = cur_start + (start[i] - plot_start[chrom[i]])
        rect_width = end[i] - start[i]
        ax_genome.add_patch(patches.Rectangle((rect_start, 0), width=rect_width, height=0.8, color='black', zorder=1, alpha=.4, linewidth=0))
        rect_mid = rect_start + rect_width / 2
        ax_genome.text(rect_mid, 0.3, node_list[i], ha='center', va='center', fontsize=10, color='black')
    ax_genome.axis('off')
    ## plot link
    h = -0.5
    for i in SV_fragment.index:
        row = SV_fragment.loc[i]
        start1 = (row['start1']-plot_start[row['chrom1']]) + dict_start[row['chrom1']]
        end1 = (row['end1']-plot_start[row['chrom1']]) + dict_start[row['chrom1']]
        start2 = (row['start2']-plot_start[row['chrom2']]) + dict_start[row['chrom2']]
        end2 = (row['end2']-plot_start[row['chrom2']]) + dict_start[row['chrom2']]
        colors = list(mcolors.TABLEAU_COLORS.keys())
        color_list.append(mcolors.TABLEAU_COLORS[colors[i % len(colors)]])
        ax_genome.add_patch(patches.Rectangle((start1, h), width=end1-start1, height=0.2, facecolor=mcolors.TABLEAU_COLORS[colors[i % len(colors)]], zorder=1, linewidth=2, clip_on=False))
        ax_genome.add_patch(patches.Rectangle((start2, h), width=end2-start2, height=0.2, facecolor=mcolors.TABLEAU_COLORS[colors[i % len(colors)]], zorder=1, linewidth=2, clip_on=False))
        # plot link
        if row['orient1'] == '+':
            pos1 = (row['end1']-plot_start[row['chrom1']]) + dict_start[row['chrom1']]
        else:
            pos1 = (row['start1']-plot_start[row['chrom1']]) + dict_start[row['chrom1']]
        if row['orient2'] == '+':
            pos2 = (row['start2']-plot_start[row['chrom2']]) + dict_start[row['chrom2']]
        else:
            pos2 = (row['end2']-plot_start[row['chrom2']]) + dict_start[row['chrom2']]
        # height = -0.2
        ax_genome.annotate("",
                    xy=[pos1, h],
                    xytext=[pos2, h],
                    arrowprops=dict(color=mcolors.TABLEAU_COLORS[colors[i % len(colors)]],
                                    arrowstyle="-",
                                    connectionstyle="arc3,rad=-0.3",
                                    linewidth=2
                                    ), zorder=0, clip_on=False)
        h -= 0.3
    ylim = 0 - SV_fragment.shape[0] * 1.5
    ax_genome.set_ylim(ylim, 2)
    
    canvas = FigureCanvas(fig)
    return canvas

def FragmentAssembly(SV_file, resolution, dict_matrix, chrom_list, dict_count, savepath):
    merge = 1000000
    win = 10000000
    # get breakpoint
    bp_result = read_SV_info(SV_file)
    num = bp_result.shape[0]
    SV_fragment = pd.DataFrame(columns=['chrom1', 'start1', 'end1', 'orient1', 'chrom2', 'start2', 'end2', 'orient2'])
    for i in range(num):
        index = [i]
        cur_bp_result = bp_result.loc[index]
        bp_Link, overlap = get_SV_orient(cur_bp_result, dict_matrix, resolution, merge)
        for j in bp_Link.index:
            row = bp_Link.loc[j]
            # smaller that 1Mb region
            if row['pos1'] < (1000000 // resolution) or overlap:
                row_win = 1000000
            else:
                row_win = win
            if row['pos2'] < (1000000 // resolution) or overlap:
                col_win = 1000000
            else:
                col_win = win
            if row['orient1'] == '+':
                pos1_start = row['pos1'] - row_win // resolution
                pos1_end = row['pos1']
            else:
                pos1_start = row['pos1']
                pos1_end = row['pos1'] + row_win // resolution
            if row['orient2'] == '+':
                pos2_start = row['pos2']
                pos2_end = row['pos2'] + col_win // resolution
            else:
                pos2_start = row['pos2'] - col_win // resolution
                pos2_end = row['pos2']
            if pos1_start < 0:
                pos1_start = 0
            if pos2_start < 0:
                pos2_start = 0
            name = row['chrom1'] + '_' + row['chrom2']
            matrix = dict_matrix[name]
            cur_matrix = matrix[pos1_start:pos1_end, pos2_start:pos2_end]
            cut_1, cut_2 = recalculate_breakpoint(cur_matrix, row['orient1'], row['orient2'])
            if row['orient1'] == '-' and row['orient2'] == '-':
                frag1start = pos1_start
                frag1end = pos1_start + cut_1 - 1
                frag2start = pos2_start + cut_2 - 1
                frag2end = pos2_end
            elif row['orient1'] == '+' and row['orient2'] == '+':
                frag1start = pos1_start + cut_1 - 1
                frag1end = pos1_end
                frag2start = pos2_start
                frag2end = pos2_start + cut_2 - 1
            elif row['orient1'] == '+' and row['orient2'] == '-':
                frag1start = pos1_start + cut_1 - 1
                frag1end = pos1_end
                frag2start = pos2_start + cut_2 - 1
                frag2end = pos2_end
            else:
                frag1start = pos1_start
                frag1end = pos1_start + cut_1 - 1
                frag2start = pos2_start
                frag2end = pos2_start + cut_2 - 1
            SV_fragment.loc[len(SV_fragment.index)] = [row['chrom1'], frag1start, frag1end, row['orient1'],
                                                       row['chrom2'], frag2start, frag2end, row['orient2']]
    SV_fragment = coverFragment(SV_fragment)
    all_path, all_orient = constructPath(SV_fragment)
    chrom, start, end, node_list = genomeFragment(SV_fragment)

    dict_start = {}
    dict_end = {}
    for i in range(len(chrom_list)):
        cur_pos = []
        for j in range(len(chrom)):
            if chrom[j] == chrom_list[i]:
                cur_pos.append(start[j])
                cur_pos.append(end[j])
        dict_start[chrom_list[i]] = max(0, min(cur_pos)-5000000//resolution)
        dict_end[chrom_list[i]] = min(max(cur_pos) + 5000000 // resolution, dict_count[chrom_list[i]])

    # plot Hi-C matrix
    chrom_matrix, interval, chrom_text_pos, chrom_text = combine_matrix(dict_matrix, dict_start, dict_end, chrom_list, resolution)
    canvas_matrix = matrix_plot(chrom_matrix, interval, chrom_text_pos, chrom_text, 
                                SV_fragment, chrom, start, end, node_list, chrom_list, dict_start, dict_end)
    
    # plot assembly
    canvas_assembly, CGR_info = plot_SV_genome(all_path, all_orient, SV_fragment)
    return canvas_matrix, canvas_assembly, CGR_info


