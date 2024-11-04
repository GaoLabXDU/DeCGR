import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle


class VisualSVCanvas(FigureCanvas):
    def __init__(self, matrix, interval, chrom_list, start_list, end_list, SV_list, resolution, output_widget):
        self.fig, (self.ax_genome, self.ax_heatmap) = plt.subplots(2, 1, 
                                                                 gridspec_kw={'height_ratios': [1, 10]}, 
                                                                 figsize=(10, 10), sharex=True)
        super(VisualSVCanvas, self).__init__(self.fig)
        self.output_widget = output_widget
        self.interval = interval
        self.chrom_list = chrom_list
        self.start_list = start_list
        self.end_list = end_list
        self.resolution = resolution
        cdict = {
            'red': ((0.0, 0.0, 1.0), (1.0, 1.0, 1.0)),
            'green': ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
            'blue': ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
        }
        cmap = LinearSegmentedColormap('Rd_Bl_Rd', cdict, 256)
        
        self.data = matrix
        vmax = np.percentile(matrix, 95)
        vmin = matrix.min()
        self.image = self.ax_heatmap.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax)
        #self.ax_heatmap.axis('off')
        self.ax_heatmap.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        row = matrix.shape[0]
        # plot chrom boundary
        for i in range(len(interval)):
            if interval[i] != 0:
                self.ax_heatmap.plot([0, row], [interval[i], interval[i]],
                        color="black", linestyle='--', linewidth=2)
                self.ax_heatmap.plot([interval[i], interval[i]], [0, row],
                        color="black", linestyle='--', linewidth=2)
        # plot SV block
        block_x_start, block_y_start = SV_list[0], SV_list[2]
        block_width, block_height = SV_list[3]-SV_list[2], SV_list[1]-SV_list[0]
        if SV_list[0] == SV_list[1]:
            self.ax_heatmap.add_patch(Circle((block_y_start, block_x_start), radius=0.5,
                                     edgecolor='black', facecolor='none',
                                     linestyle='--', linewidth=2))
        else:
            self.ax_heatmap.add_patch(Rectangle((block_y_start, block_x_start), block_width, block_height, 
                                    edgecolor='black', facecolor='none', 
                                    linestyle='--', linewidth=2))
        # plot genome region
        start_pos = 0
        for k in range(len(start_list)):
            width = (end_list[k] - start_list[k]) // resolution
            self.ax_genome.add_patch(FancyBboxPatch((start_pos, 0), width, 1,
                                               boxstyle="round,pad=0.3", 
                                            facecolor=(245/255, 195/255, 66/255), edgecolor='black', alpha=0.7))
            self.ax_genome.text(start_pos + (start_list[k] + end_list[k]) // (2 * resolution), 0.5, chrom_list[k], ha='center', va='center')   
            self.ax_genome.text(start_pos, 1.05, f'{round(start_list[k]/1000000, 2)}', 
                    ha='left', va='bottom', fontsize=10, color='black')
            self.ax_genome.text(start_pos + width, 1.05, f'{round(end_list[k]/1000000, 2)}M', 
                    ha='right', va='bottom', fontsize=10, color='black')
            start_pos += width

        self.ax_heatmap.set_xlim(0, matrix.shape[1])
        self.ax_genome.set_xlim(0, matrix.shape[1])
        self.ax_genome.axis('off')
        self.ax_genome.set_position([0.1, 0.9, 0.8, 0.03])
        self.ax_heatmap.set_position([0.1, 0.1, 0.8, 0.8])
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
        x_left, x_right = self.ax_heatmap.get_xlim()
        y_bottom, y_top = self.ax_heatmap.get_ylim()
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

        self.ax_heatmap.set_xlim([x_left, x_right])
        self.ax_heatmap.set_ylim([y_bottom, y_top])

        self.draw()

    def on_press(self, event):
        if event.button == 1:
            if event.inaxes != self.ax_heatmap:
                return
            self.press = (event.xdata, event.ydata)
            if event.dblclick:
                self.on_click(event) 

    def on_move(self, event):
        if self.press is None:
            return
        if event.xdata is None or event.ydata is None:
            return
        x_left, x_right = self.ax_heatmap.get_xlim()
        y_bottom, y_top = self.ax_heatmap.get_ylim()
        dx = event.xdata - self.press[0]
        dy = event.ydata - self.press[1]
        self.ax_heatmap.set_xlim(x_left - dx, x_right - dx)
        self.ax_heatmap.set_ylim(y_bottom - dy, y_top - dy)
        self.press = (event.xdata, event.ydata)
        self.draw()

    def on_release(self, event):
        self.press = None

    def on_click(self, event):
        if event.inaxes != self.ax_heatmap:
            return
        x_click = int(event.xdata)
        y_click = int(event.ydata)
        output_text = ""
        if 0 <= x_click < self.data.shape[1] and 0 <= y_click < self.data.shape[0]:
            pixel_value = self.data[y_click, x_click]
            if len(self.interval) == 0:
                # intra
                cur_chrom = self.chrom_list[0]
                cur_start = self.start_list[0] + y_click * self.resolution
                cur_end = self.start_list[0] + x_click * self.resolution
                output_text += "position1: " + "\t"  + str(cur_chrom) + "\t" + str(cur_start) + "\n"
                output_text += "position2: " + "\t"  + str(cur_chrom) + "\t" + str(cur_end) + "\n"
            else:
                # inter
                if x_click <= self.interval[0] and y_click <= self.interval[0]:
                    # region_A
                    cur_start = self.start_list[0] + y_click * self.resolution
                    cur_end = self.start_list[0] + x_click * self.resolution
                    output_text += "position1: " + "\t"  + str(self.chrom_list[0]) + "\t" + str(cur_start) + "\n"
                    output_text += "position2: " + "\t"  + str(self.chrom_list[0]) + "\t" + str(cur_end) + "\n"
                elif x_click > self.interval[0] and y_click > self.interval[0]:
                    # region_B
                    cur_start = self.start_list[1] + y_click * self.resolution
                    cur_end = self.start_list[1] + x_click * self.resolution
                    output_text += "position1: " + "\t"  + str(self.chrom_list[0]) + "\t" + str(cur_start) + "\n"
                    output_text += "position2: " + "\t"  + str(self.chrom_list[0]) + "\t" + str(cur_end) + "\n"
                elif x_click <= self.interval[0] and y_click > self.interval[0]:
                    # region_C
                    cur_start = self.start_list[1] + (y_click-self.interval[0]) * self.resolution
                    cur_end = self.start_list[0] + x_click * self.resolution
                    output_text += "position1: " + "\t"  + str(self.chrom_list[0]) + "\t" + str(cur_start) + "\n"
                    output_text += "position2: " + "\t"  + str(self.chrom_list[1]) + "\t" + str(cur_end) + "\n"
                else:
                    # region_D
                    cur_start = self.start_list[0] + y_click * self.resolution
                    cur_end = self.start_list[1] + (x_click-self.interval[0]) * self.resolution
                    output_text += "position1: " + "\t"  + str(self.chrom_list[0]) + "\t" + str(cur_start) + "\n"
                    output_text += "position2: " + "\t"  + str(self.chrom_list[1]) + "\t" + str(cur_end) + "\n"
            output_text += "Count: " + "\t" + str(pixel_value)
            self.output_widget.setText(output_text)


def plot_current_SV(matrix, interval, chrom_list, start_list, end_list, SV_list, resolution, output_widget):
    canvas = VisualSVCanvas(matrix, interval, chrom_list, start_list, end_list, SV_list, resolution, output_widget)
    return canvas


def extra_SV_region_intra(clr, region):
    matrix = clr.matrix(balance=False).fetch(region)
    return matrix, []

def extra_SV_region_inter(clr, region1, region2):
    matrix1 = clr.matrix(balance=False).fetch(region1)
    matrix2 = clr.matrix(balance=False).fetch(region2)
    matrix_inter = clr.matrix(balance=False).fetch(region1, region2)
    # combine these matrix into one matrix
    n1 = matrix1.shape[0]
    n2 = matrix2.shape[0]
    final_matrix = np.zeros((n1 + n2, n1 + n2))
    final_matrix[:n1, :n1] = matrix1
    final_matrix[n1:, n1:] = matrix2
    final_matrix[:n1, n1:] = matrix_inter
    final_matrix[n1:, :n1] = matrix_inter.T
    return final_matrix, [n1]
