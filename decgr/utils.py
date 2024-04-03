import os
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import operator
import itertools
import matplotlib.pyplot as plt


def load_matrix(clr, chr_pair):
    M = clr.matrix(balance=False).fetch(chr_pair[0], chr_pair[1])
    M[np.isnan(M)] = 0
    return M

def load_bias(clr, chr_pair):
    bias = clr.bins().fetch(chr_pair[0])['sweight'].values
    # bias = clr.matrix(balance='sweight').fetch(chr_pair[0], chr_pair[1])
    bias[np.isnan(bias)] = 0
    print(bias)
    return bias

def get_MatrixFile(hic_path, name):
    matrix = 0
    for filewalks in os.walk(hic_path):
        for files in filewalks[2]:
            if name in files:
                matrix = np.loadtxt(os.path.join(filewalks[0], files))
    return matrix


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


def runFitHiC(path, control_dict, chrom_list, chro_lens, resolution):
    with open(path + '/tmp/fithic.bed', 'w+') as f:
        chr_first = {}
        index = 1
        for idx, chro_len in enumerate(chro_lens):
            chr = chrom_list[idx]
            chr_first[chr] = index
            start = 0
            for _ in range(chro_len):
                f.write(chr + '\t' + str(start) + '\t' + str(start + resolution) + '\t' + str(index) + '\n')
                start += resolution
                index += 1

    with open(path + '/tmp/fithic.matrix', 'w+') as f:
        for i in range(len(chrom_list)):
            for j in range(i, len(chrom_list)):
                chr1_first = chr_first[chrom_list[i]]
                chr2_first = chr_first[chrom_list[j]]
                name = str(chrom_list[i]) + '_' + str(chrom_list[j])
                cur_matrix = control_dict[name]
                for bin_i in range(cur_matrix.shape[0]):
                    for bin_j in range(bin_i, cur_matrix.shape[1]):
                        if cur_matrix[bin_i][bin_j] > 0:
                            f.write(str(chr1_first + bin_i) + '\t' + str(chr2_first + bin_j) + '\t' + str(cur_matrix[bin_i][bin_j]) + '\n')

    os.system('python {0}/fithic/HiCPro2FitHiC.py -i {1} -b {2} -o {0}/tmp -r {3}'.format(
        path, path + '/tmp/fithic.matrix', path + '/tmp/fithic.bed', str(resolution)))
    os.system(
        'python {0}/fithic/fithic_copy.py -i {0}/tmp/fithic.interactionCounts.gz -f {0}/tmp/fithic.fragmentMappability.gz -o {0}/tmp -r {1} -x All'
        .format(path, str(resolution)))

    return path + '/result/dis2cnt'


