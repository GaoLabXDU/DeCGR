import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LinearRegression
import itertools, cooler
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import decgr.correct

class MyFigureCanvas(FigureCanvas):
    def __init__(self, width, height, dpi, matrix, interval, orient, path):
        self.fig = plt.figure(figsize=(width, height), dpi=dpi)
        super(MyFigureCanvas, self).__init__(self.fig)
        self.vmax = np.percentile(matrix, 98)
        self.vmin = matrix.min()
        cdict = {
            'red': ((0.0, 0.0, 1.0), (1.0, 1.0, 1.0)),
            'green': ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
            'blue': ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
        }
        self.cmap = LinearSegmentedColormap('Rd_Bl_Rd', cdict, 256)
        self.axes = self.fig.add_subplot(111)
        self.axes.axis('off')
        row = matrix.shape[0]
        M = matrix
        n = M.shape[0]
        # Create the rotation matrix
        t = np.array([[1, 0.5], [-1, 0.5]])
        A = np.dot(np.array([(i[1], i[0]) for i in itertools.product(range(n, -1, -1), range(0, n + 1, 1))]), t)

        # Plot the Heatmap ...
        x = A[:, 1].reshape(n + 1, n + 1)
        y = A[:, 0].reshape(n + 1, n + 1)
        y[y < 0] = -y[y < 0]
        self.axes.pcolormesh(x, y, np.flipud(M), vmin=self.vmin, vmax=self.vmax, cmap=self.cmap,
                      edgecolor='none', snap=True, linewidth=.01, rasterized=True)
        # plot boundary
        for i in range(len(interval)):
            if interval[i] != 0:
                self.axes.plot([x[(n - interval[i]), 0], x[(n - interval[i]), interval[i]]], [y[(n - interval[i]), 0], 0],
                        color="black", linestyle='--', linewidth=2)
                self.axes.plot([x[(n - interval[i]), n], x[(n - interval[i]), interval[i]]], [y[(n - interval[i]), n], 0],
                        color="black", linestyle='--', linewidth=2)
        self.axes.plot([0, x[0, 0]], [0, y[0, 0]], color="black", linewidth=2)
        self.axes.plot([x[0, 0], n], [y[0, 0], 0], color="black", linewidth=2)

        # plot CSV
        h = -10
        interval.append(n)
        for i in range(len(interval) - 1):
            if orient[i] == '+':
                pos1 = interval[i + 1]
                pos2 = interval[i]
            else:
                pos1 = interval[i]
                pos2 = interval[i + 1]
            textPos = interval[i] + (interval[i + 1] - interval[i]) // 2
            plt.text(textPos, h + 4, path[i],
                     fontdict={'family': 'Arial', 'size': 16, 'color': 'black', 'weight': 'bold'})
            self.axes.annotate("",
                        xy=[pos1, h],
                        xytext=[pos2, h],
                        size=10,
                        arrowprops=dict(color='gray', headwidth=20, width=10, headlength=10, ec="black", lw=1),
                        zorder=0)

        self.axes.axis('off')
        plt.ylim([-10, n])
        plt.subplots_adjust(0.1, 0.1, 0.9, 0.8)
        self.isshowing=False


def plot_assembly_matrix(matrix, interval, path, orient):
    canvas = MyFigureCanvas(width=5, height=4, dpi=100, matrix=matrix, interval=interval, orient=orient, path=path)
    return canvas


def joint_SV(path, orient, chrom, start, end, dict_matrix):
    genome = []
    SV_chrom = []
    for i in range(len(path)):
        if orient[i] == '+':
            fragment = list(range(start[i], end[i]))
            cur_chrom = [chrom[i]] * (end[i] - start[i])
            genome.extend(fragment)
            SV_chrom.extend(cur_chrom)
        else:
            fragment = list(range(end[i], start[i], -1))
            cur_chrom = [chrom[i]] * (end[i] - start[i])
            SV_chrom.extend(cur_chrom)
            genome.extend(fragment)
    # init_SV_matrix
    rows = len(genome)
    SV_matrix = np.full((rows, rows), 0.0)
    for i in range(len(genome)):
        for j in range(i, len(genome)):
            chrom_pair = SV_chrom[i] + '_' + SV_chrom[j]
            cur_matrix = dict_matrix[chrom_pair]
            SV_matrix[i][j] = cur_matrix[genome[i]][genome[j]]

    return SV_matrix


def get_start(path, start, end):
    cur_start = 0
    if len(path) > 1:
        for i in range(len(path)-1):
            # index = ord(path[i]) - 65
            cur_start += (end[i] - start[i])
    return cur_start


def linear_regression(SV_exp, global_exp):
    X = []
    Y = []
    Xi = []
    for i in sorted(SV_exp):
        if i in global_exp:
            Xi.append(i)
            X.append(global_exp[i])
            Y.append(SV_exp[i])
    X = np.r_[X]
    Y = np.r_[Y]

    IR = IsotonicRegression(increasing=False, out_of_bounds='clip')
    IR.fit(Xi, Y)
    Y1 = IR.predict(Xi)
    filter_X = []
    filter_Y = []
    for i in range(len(Y1)):
        if Y1[i] not in filter_Y:
            filter_Y.append(Y1[i])
            filter_X.append([X[i]])

    reg = LinearRegression()
    reg.fit(filter_X, filter_Y)
    slope = reg.coef_[0]

    return slope


def scale_matrix(matrix, start1, end1, start2, end2, glob_exp):
    SV_exp = {}
    dist_count = {}
    for i in range(start1, end1):
        for j in range(start2, end2):
            dist = abs(j - i)
            if dist not in SV_exp:
                SV_exp[dist] = 0
                dist_count[dist] = 0
            SV_exp[dist] += matrix[i][j]
            dist_count[dist] += 1
    filter_exp = {}
    for key in SV_exp:
        if dist_count[key] > 2:
            filter_exp[key] = SV_exp[key] / dist_count[key]
        # SV_exp[key] /= dist_count[key]
    if filter_exp:
        slope = linear_regression(filter_exp, glob_exp)
        if slope > 0:
            for i in range(start1, end1):
                for j in range(start2, end2):
                    matrix[i][j] /= slope

    return matrix


def adjust_matrix(matrix, glob_exp, path, start, end):
    interval = []
    for i in range(len(path)):
        row_start = get_start(path[:i+1], start, end)
        if row_start not in interval:
            interval.append(row_start)
        for j in range(i+1, len(path)):
            col_start = get_start(path[:j+1], start, end)
            start1 = row_start
            end1 = row_start+(end[i]-start[i])
            start2 = col_start
            end2 = col_start+(end[j]-start[j])
            matrix = scale_matrix(matrix, start1, end1, start2, end2, glob_exp)

    return matrix, interval


def read_one_assembly(filaname, resolution):
    path = []
    orient = []
    chrom = []
    start = []
    end = []
    with open(filaname) as f:
        f.readline()
        for line in f:
            line_data = line.strip().split('\t')
            chrom.append(line_data[0])
            start.append(int(line_data[1]) // resolution)
            end.append(int(line_data[2]) // resolution)
            path.append(line_data[3])
            orient.append(line_data[4])

    return path, orient, chrom, start, end


def Calculating_diagonal_data(matrix):
    """
    Calculating each diagonal mean
    """
    exp = {}
    N, M = len(matrix), len(matrix[0])
    for d in range(N):
        intermediate = []
        c = d
        r = 0
        while r < N - d:
            intermediate.append(matrix[r][c])
            r += 1
            c += 1
        intermediate = np.array(intermediate)
        exp[d] = np.mean(intermediate)
    return exp


def global_expect(dict_matrix, chrom_list):
    chrom = chrom_list[0] + '_' + chrom_list[0]
    matrix = dict_matrix[chrom]
    global_exp = Calculating_diagonal_data(matrix)
    IR = IsotonicRegression(increasing=False, out_of_bounds='clip')
    _d = np.r_[sorted(global_exp)]
    _y = np.r_[[global_exp[i] for i in _d]]
    IR.fit(_d, _y)
    d = np.arange(500 + 1)
    exp = IR.predict(d)
    expected = dict(zip(d, exp))
    return expected