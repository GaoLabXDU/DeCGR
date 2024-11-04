import os
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import operator
import matplotlib.pyplot as plt


def load_matrix(clr, chr_pair):
    M = clr.matrix(balance=False).fetch(chr_pair[0], chr_pair[1])
    M[np.isnan(M)] = 0
    return M

def load_norm_matrix(clr, chr_pair):
    M = clr.matrix(balance='sweight').fetch(chr_pair[0], chr_pair[1])
    M[np.isnan(M)] = 0
    return M

def load_bias(clr, chr_pair):
    bias = clr.bins().fetch(chr_pair[0])['sweight'].values
    # bias = clr.matrix(balance='sweight').fetch(chr_pair[0], chr_pair[1])
    bias[np.isnan(bias)] = 0
    return bias

def get_MatrixFile(hic_path, name):
    matrix = 0
    for filewalks in os.walk(hic_path):
        for files in filewalks[2]:
            if name in files:
                matrix = np.loadtxt(os.path.join(filewalks[0], files))
    return matrix



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
    current_path = os.path.dirname(os.path.abspath(__file__))
    fithic_path = os.path.join(current_path, '../decgr/fithic')

    os.system('python {0}/HiCPro2FitHiC.py -i {2} -b {3} -o {1}/tmp -r {4}'.format(
        fithic_path, path, path + '/tmp/fithic.matrix', path + '/tmp/fithic.bed', str(resolution)))
    os.system(
        'python {0}/fithic_copy.py -i {1}/tmp/fithic.interactionCounts.gz -f {1}/tmp/fithic.fragmentMappability.gz -o {1}/tmp -r {2} -x All'
        .format(fithic_path, path, str(resolution)))

    return path + '/result/dis2cnt'


