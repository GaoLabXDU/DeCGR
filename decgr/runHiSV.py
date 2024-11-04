# running HiSV
import cooler
import numpy as np
import numpy as np
import math
from itertools import groupby
from collections import defaultdict
from skimage.restoration import denoise_tv_chambolle

def combine(lst):
    # Merge consecutive identical values
    pos = (j - i for i, j in enumerate(lst))
    t = 0
    for i, els in groupby(pos):
        l = len(list(els))
        el = lst[t]
        t += l
        yield range(el, el+l)


def group(data):
    '''
    # Merge consecutive identical values
    :param data: pos index
    :return: start and end of consecutive same values
    '''
    last = data[0]
    start = end = 0
    for n in data[1:]:
        if n - last <= 5: # Part of the group, bump the end
            last = n
            end += 1
        else: # Not part of the group, yield current group and start a new
            yield range(data[start], data[end]+1)
            last = n
            start = end = end + 1
    # yield start, end
    yield range(data[start], data[end]+1)



def group_position(pos1list, pos2list, binsize):
    result = defaultdict(list)
    l1 = sorted(set(pos1list), key=pos1list.index)
    pos1group = list(group(l1))
    for i in range(len(pos1group)):
        start_pos = pos1list.index(pos1group[i][0])
        end_pos = len(pos1list) - 1 - pos1list[::-1].index(pos1group[i][-1])
        cur_pos2_list = pos2list[start_pos:end_pos + 1]
        cur_pos1_list = pos1list[start_pos:end_pos + 1]
        cur_pos2_list_sort = sorted(cur_pos2_list)
        l2 = sorted(set(cur_pos2_list_sort), key=cur_pos2_list_sort.index)
        pos2group = list(group(l2))
        if len(pos2group) > 1:
            for j in range(len(pos2group)):
                new_pos1 = []
                for k in range(len(pos2group[j])):
                    if pos2group[j][k] in cur_pos2_list:
                        pos = cur_pos2_list.index(pos2group[j][k])
                        new_pos1.append(cur_pos1_list[pos])        
                new_pos1 = sorted(new_pos1)
                result['pos1_start'].append(new_pos1[0] * binsize)
                result['pos1_end'].append(new_pos1[-1] * binsize + binsize)
                result['pos2_start'].append(pos2group[j][0] * binsize)
                result['pos2_end'].append(pos2group[j][-1] * binsize + binsize)

        else:
            result['pos1_start'].append(pos1group[i][0] * binsize)
            result['pos1_end'].append(pos1group[i][-1] * binsize + binsize)
            result['pos2_start'].append(pos2group[0][0] * binsize)
            result['pos2_end'].append(pos2group[0][-1] * binsize + binsize)
    return result


def load_matrix(clr, chrom1, chrom2):
    # load matrix from cool file (normalized matrix)
    if chrom1 == chrom2:
        M = clr.matrix(balance=False).fetch(chrom1)
    else:
        M = clr.matrix(balance=False).fetch(chrom1, chrom2)
    M[np.isnan(M)] = 0
    return M
    
def cool2mat(clr, chrom1, chrom2):
    chromsomes = clr.chromnames
    if chrom1 not in chromsomes or chrom2 not in chromsomes:
        raise ValueError(f"Chromosome number is invalidation.")
    else:
        M = load_matrix(clr, chrom1, chrom2)
    return M


def local_saliency(matrix, win):
    """
    # calculate local saliency map
    :param matrix: current Hi-C contact matrix
    :param win: 'Local region window.
    :return: saliency map
    """
    d = int(win / 2)
    N, M = len(matrix), len(matrix[0])
    nm_matrix = np.full((N, M), 0.0)
    Dist_matrix = np.full((win, win), 0.0)

    for m in range(win):
        for n in range(win):
            Dist_matrix[m][n] = np.sqrt(np.square(m) + np.square(n))

    it = np.nditer(matrix, flags=['multi_index'])
    while not it.finished:
        idx = it.multi_index
        data = matrix[idx]
        if data != 0:
            cur_matrix = matrix[idx[0] - d:idx[0] + d, idx[1] - d:idx[1] + d]
            cur_raw = cur_matrix.shape[0]
            cur_col = cur_matrix.shape[1]
            if cur_col == cur_raw == win and np.sum(cur_matrix) != 0:
                nm_matrix[idx] = 1 - math.exp(-np.mean(abs(data - cur_matrix) / (1 + Dist_matrix)))
            else:
                nm_matrix[idx] = 0
        it.iternext()
    return nm_matrix


def Calculating_diagonal_data(matrix):
    """
    # normalization matrix by diagonal to remove distance bias
    Calculating each diagonal mean and std
    """
    N, M = len(matrix), len(matrix[0])
    Diagonal_mean = np.full(M, 0.0)
    Diagonal_std = np.full(M, 0.0)
    std = []
    for d in range(N):
        intermediate = []
        c = d
        r = 0
        while r < N - d:
            intermediate.append(matrix[r][c])
            r += 1
            c += 1
        intermediate = np.array(intermediate)
        Diagonal_mean[d] = (np.mean(intermediate))
        Diagonal_std[d] = (np.std(intermediate))
    return Diagonal_mean, Diagonal_std


def Distance_normalization(matrix):
    """
    # normalization matrix by diagonal to remove distance bias
    norm_data = (data - mean_data) / mean_std
    """
    Diagonal_mean, Diagonal_std = Calculating_diagonal_data(matrix)
    N, M = len(matrix), len(matrix[0])
    for d in range(N):
        c = d
        r = 0
        while r < N - d:
            if Diagonal_std[d] == 0:
                matrix[r][c] = 0
            else:
                if matrix[r][c] - Diagonal_mean[d] < 0:
                    matrix[r][c] = 0
                else:
                    matrix[r][c] = (matrix[r][c] - Diagonal_mean[d]) / Diagonal_std[d]
            r += 1
            c += 1
    return matrix


def HiSV(clr, chromList, win, weight, cutoff):
    # load matrix from cool file
    combine_result = defaultdict(list)
    resolution = clr.info['bin-size']
    for i in range(len(chromList)):
        for j in range(i, len(chromList)):
            chrom1 = chromList[i]
            chrom2 = chromList[j]
            mat = cool2mat(clr, chrom1, chrom2)
            # running HiSV
            pos1 = []
            pos2 = []
            # filter low reads
            if chrom1 == chrom2:
            # intra-chromosomes SV
                # filter low reads
                non_zero_rows_and_cols = np.where(np.any(mat != 0, axis=1))[0]
                mat = mat[non_zero_rows_and_cols][:, non_zero_rows_and_cols]
                num = np.shape(mat)[0]
                il = np.tril_indices(num)
                mat[il] = 0
                # z-score for distance normalization
                Dist_norm_mat = Distance_normalization(mat)
                # calculating local saliency
                local_sali_matrix = local_saliency(Dist_norm_mat, win)
                # Segmentation
                seg_mat = denoise_tv_chambolle(local_sali_matrix, weight=weight)
                # Filter
                for m in range(num):
                    for n in range(num):
                        if seg_mat[m][n] > cutoff:
                            pos1.append(non_zero_rows_and_cols[m])
                            pos2.append(non_zero_rows_and_cols[n])
                result = group_position(pos1, pos2, resolution)
            else:
                # inter-chromosomes SV
                # filter
                non_zero_rows = np.where(np.any(mat != 0, axis=1))[0]
                non_zero_cols = np.where(np.any(mat != 0, axis=0))[0]
                mat = mat[non_zero_rows][:, non_zero_cols]
                raw_num, col_num = np.shape(mat)
                # calculating local saliency
                local_sali_matrix = local_saliency(mat, win)
                # Segmentation
                seg_mat = denoise_tv_chambolle(local_sali_matrix, weight=0.2)
                # Filter
                for m in range(raw_num):
                    for n in range(col_num):
                        if seg_mat[m][n] > cutoff:
                            pos1.append(non_zero_rows[m])
                            pos2.append(non_zero_cols[n])
                result = group_position(pos1, pos2, resolution)
            result['chrom1'] = [chrom1] * len(result['pos1_start'])
            result['chrom2'] = [chrom2] * len(result['pos2_start']) 
            for key, values in result.items():
                combine_result[key].extend(values)
    return dict(combine_result)
