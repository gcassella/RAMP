import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.widgets as wid
from mpl_toolkits.mplot3d import Axes3D

import numpy as np

mpl.rcParams['toolbar'] = 'None'
plt.style.use('dark_background')

class Visualisation():
    def __init__(self, inst, controls=True):
        self.inst = inst
        self.controls = controls

        self.fig = plt.figure()
        self.ax_zx = self.fig.add_subplot(2, 2, 1)
        self.ax_zx.set_xlabel('z [m]')
        self.ax_zx.set_ylabel('x [m]')

        self.ax_zy = self.fig.add_subplot(2, 2, 2)
        self.ax_zy.set_xlabel('z [m]')
        self.ax_zy.set_ylabel('y [m]')

        self.ax_xy = self.fig.add_subplot(2, 2, 3)
        self.ax_xy.set_xlabel('x [m]')
        self.ax_xy.set_ylabel('y [m]')

        self.ax_or = self.fig.add_subplot(2, 2, 4, projection='3d')
        self.ax_or.set_xlabel('x [m]')
        self.ax_or.set_ylabel('z [m]')
        self.ax_or.set_zlabel('y [m]')

        if self.controls:
            plt.subplots_adjust(left=0.05, right=0.7)

            rectprops = dict(facecolor='blue', alpha=0.5)
            self.xlim_span_ax = plt.axes([0.72, 0.3, 0.25, 0.03])
            self.xlim_span_ax.set_yticks([])
            self.xlim_span = wid.SpanSelector(self.xlim_span_ax, self._inst_xlim_change, 'horizontal', rectprops=rectprops)
            xlim_span_label = self._create_text("xlim", [0.72, 0.35, 0.25, 0.03])

            self.ylim_span_ax = plt.axes([0.72, 0.2, 0.25, 0.03])
            self.ylim_span_ax.set_yticks([])
            self.ylim_span = wid.SpanSelector(self.ylim_span_ax, self._inst_ylim_change, 'horizontal', rectprops=rectprops)
            ylim_span_label = self._create_text("ylim", [0.72, 0.25, 0.25, 0.03])

            self.zlim_span_ax = plt.axes([0.72, 0.1, 0.25, 0.03])
            self.zlim_span_ax.set_yticks([])
            self.zlim_span = wid.SpanSelector(self.zlim_span_ax, self._inst_zlim_change, 'horizontal', rectprops=rectprops)
            zlim_span_label = self._create_text("zlim", [0.72, 0.15, 0.25, 0.03])

    def _update_limits(self):
        if self.controls:
            def_xlim = self.ax_zx.get_ylim()
            xlim_range = def_xlim[1] - def_xlim[0]
            self.xlim_span_ax.set_xlim([def_xlim[0]-0.5*xlim_range, def_xlim[1]+0.5*xlim_range])

            def_ylim = self.ax_xy.get_ylim()
            ylim_range = def_ylim[1] - def_ylim[0]
            self.ylim_span_ax.set_xlim([def_ylim[0]-0.5*ylim_range, def_ylim[1]+0.5*ylim_range])

            def_zlim = self.ax_zx.get_xlim()
            zlim_range = def_zlim[1] - def_zlim[0]
            self.zlim_span_ax.set_xlim([def_zlim[0]-0.5*zlim_range, def_zlim[1]+0.5*zlim_range])

    def _inst_xlim_change(self, vmin, vmax):
        if not vmin == vmax:
            self.ax_zx.set_ylim(vmin, vmax, emit=False, auto=False)
            self.ax_xy.set_xlim(vmin, vmax, emit=False, auto=False)
            self.ax_or.set_xlim(vmin, vmax, emit=False, auto=False)

    def _inst_ylim_change(self, vmin, vmax):
        if not vmin == vmax:
            self.ax_zy.set_ylim(vmin, vmax, emit=False, auto=False)
            self.ax_xy.set_ylim(vmin, vmax, emit=False, auto=False)
            self.ax_or.set_zlim(vmin, vmax, emit=False, auto=False)

    def _inst_zlim_change(self, vmin, vmax):
        if not vmin == vmax:
            self.ax_zx.set_xlim(vmin, vmax, emit=False, auto=False)
            self.ax_zy.set_xlim(vmin, vmax, emit=False, auto=False)
            self.ax_or.set_ylim(vmin, vmax, emit=False, auto=False)

    def _create_text(self, text, position):
        text_axes = plt.axes(position)
        text_axes.text(0, 0, text)
        text_axes.axis('off')
        return text_axes

    def plot_component(self, comp):
        s_lines = self.inst.blocks[comp.block].components[comp.comp].scat_kernel.lines()
        g_lines = self.inst.blocks[comp.block].components[comp.comp].geom_kernel.lines()

        g_lines_xz = [
            np.add(g_lines[2], comp.pos[2]),
            np.add(g_lines[0], comp.pos[0])
        ]
        self.ax_zx.plot(g_lines_xz[0], g_lines_xz[1])

        g_lines_zy = [
            np.add(g_lines[2], comp.pos[2]),
            np.add(g_lines[1], comp.pos[1])
        ]
        self.ax_zy.plot(g_lines_zy[0], g_lines_zy[1])

        g_lines_xy = [
            np.add(s_lines[0], comp.pos[0]),
            np.add(s_lines[1], comp.pos[1])
        ]
        self.ax_xy.plot(g_lines_xy[0], g_lines_xy[1])

        g_lines_xyz = [
            np.add(g_lines[0], comp.pos[0]),
            np.add(g_lines[1], comp.pos[1]),
            np.add(g_lines[2], comp.pos[2])
        ]
        self.ax_or.plot(g_lines_xyz[0],
                   g_lines_xyz[2],
                   g_lines_xyz[1])

        s_lines_xz = [
            np.add(s_lines[2], comp.pos[2]),
            np.add(s_lines[0], comp.pos[0])
        ]
        self.ax_zx.plot(s_lines_xz[0], s_lines_xz[1])

        s_lines_zy = [
            np.add(s_lines[2], comp.pos[2]),
            np.add(s_lines[1], comp.pos[1])
        ]
        self.ax_zy.plot(s_lines_zy[0], s_lines_zy[1])

        s_lines_xy = [
            np.add(s_lines[0], comp.pos[0]),
            np.add(s_lines[1], comp.pos[1])
        ]
        self.ax_xy.plot(s_lines_xy[0], s_lines_xy[1])

        s_lines_xyz = [
            np.add(s_lines[0], comp.pos[0]),
            np.add(s_lines[1], comp.pos[1]),
            np.add(s_lines[2], comp.pos[2])
        ]
        self.ax_or.plot(s_lines_xyz[0],
                   s_lines_xyz[2],
                   s_lines_xyz[1])

        self._update_limits()

    def show(self):
        plt.show()
