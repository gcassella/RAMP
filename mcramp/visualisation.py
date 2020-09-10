import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.widgets as wid
from mpl_toolkits.mplot3d import Axes3D

import numpy as np

mpl.rcParams['toolbar'] = 'None'

class Visualisation():
    def __init__(self, inst, controls=True, xlim=None, ylim=None,
                 zlim=None, focus=None, style='dark_background',
                 **kwargs):
        plt.style.use(style)

        self.inst = inst
        self.controls = controls
        self._populate_offset_dict()

        self.fig = plt.figure(**kwargs)
        self.ax_zx = self.fig.add_subplot(2, 2, 1)
        self.ax_zy = self.fig.add_subplot(2, 2, 2)
        self.ax_xy = self.fig.add_subplot(2, 2, 3)
        self.ax_or = self.fig.add_subplot(2, 2, 4, projection='3d')
        self._draw_labels()

        self.offset = [0, 0, 0]

        if self.controls:
            plt.subplots_adjust(left=0.05, right=0.7)

            rectprops = dict(facecolor='blue', alpha=0.5)
            self.xlim_span_ax = plt.axes([0.72, 0.3, 0.25, 0.03])
            self.xlim_span_ax.set_yticks([])
            self.xlim_span = wid.SpanSelector(self.xlim_span_ax, self._inst_xlim_change, 'horizontal', rectprops=rectprops)
            xlim_span_label = self._create_text("xlim", [0.72, 0.25, 0.25, 0.03])

            self.ylim_span_ax = plt.axes([0.72, 0.2, 0.25, 0.03])
            self.ylim_span_ax.set_yticks([])
            self.ylim_span = wid.SpanSelector(self.ylim_span_ax, self._inst_ylim_change, 'horizontal', rectprops=rectprops)
            ylim_span_label = self._create_text("ylim", [0.72, 0.15, 0.25, 0.03])

            self.zlim_span_ax = plt.axes([0.72, 0.1, 0.25, 0.03])
            self.zlim_span_ax.set_yticks([])
            self.zlim_span = wid.SpanSelector(self.zlim_span_ax, self._inst_zlim_change, 'horizontal', rectprops=rectprops)
            zlim_span_label = self._create_text("zlim", [0.72, 0.05, 0.25, 0.03])

            self.comp_focus_textbox_ax = plt.axes([0.72, 0.35, 0.25, 0.63], facecolor='grey')
            self.comp_focus_buttons = wid.RadioButtons(self.comp_focus_textbox_ax, tuple(kr.comp_name for kr in self.inst.kernel_refs))
            self.comp_focus_buttons.on_clicked(self._inst_comp_focus)
            #self.comp_focus_textbox = wid.TextBox(self.comp_focus_textbox_ax, "", color='.1', hovercolor='.15')
            #self.comp_focus_textbox.on_text_change(self._inst_comp_focus)
            #comp_focus_label = self._create_text("Focussed component", [0.72, 0.45, 0.25, 0.03]

        if focus:
            self._inst_comp_focus(focus)

        if xlim:
            self._inst_xlim_change(xlim[0], xlim[1])

        if ylim:
            self._inst_ylim_change(ylim[0], ylim[1])
        
        if zlim:
            self._inst_zlim_change(zlim[0], zlim[1])

    def _populate_offset_dict(self):
        self.offset_dict = {}
        for kr in self.inst.kernel_refs:
            self.offset_dict[kr.comp_name] = kr.pos

    def _draw_labels(self):
        self.ax_zx.set_xlabel('z [m]')
        self.ax_zx.set_ylabel('x [m]')
        self.ax_zy.set_xlabel('z [m]')
        self.ax_zy.set_ylabel('y [m]')
        self.ax_xy.set_xlabel('x [m]')
        self.ax_xy.set_ylabel('y [m]')
        self.ax_or.set_xlabel('x [m]')
        self.ax_or.set_ylabel('z [m]')
        self.ax_or.set_zlabel('y [m]')

    def _update_limits(self):
        def_xlim = self.ax_zx.get_ylim()
        self.xlim_range = def_xlim[1] - def_xlim[0]
        def_zlim = self.ax_zx.get_xlim()
        self.zlim_range = def_zlim[1] - def_zlim[0]
        def_ylim = self.ax_xy.get_ylim()
        self.ylim_range = def_ylim[1] - def_ylim[0]

        if self.controls:
            self.xlim_span_ax.set_xlim([(def_xlim[0]-0.5*self.xlim_range), (def_xlim[1]+0.5*self.xlim_range)])
            self.ylim_span_ax.set_xlim([(def_ylim[0]-0.5*self.ylim_range), (def_ylim[1]+0.5*self.ylim_range)])
            self.zlim_span_ax.set_xlim([(def_zlim[0]-0.5*self.zlim_range), (def_zlim[1]+0.5*self.zlim_range)])

    def _update_offset(self, offset):
        self.offset = offset
        self.ax_or.clear()
        self.ax_xy.clear()
        self.ax_zx.clear()
        self.ax_zy.clear()
        self._plot_inst()
        self._update_limits()
        self._draw_labels()
        plt.draw()

    def _inst_comp_focus(self, comp_name):
        for kr in self.inst.kernel_refs:
            if kr.comp_name == comp_name:
                list_pos = [kr.comp["pos"]["s0"], kr.comp["pos"]["s1"], kr.comp["pos"]["s2"]]
                self._update_offset(np.multiply(list_pos, -1.0))

                self._inst_xlim_change(-self.xlim_range*0.1, self.xlim_range*0.1)
                self._inst_zlim_change(-self.zlim_range*0.1, self.zlim_range*0.1)

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

    def _update_plots(self):
        self.ax_zx.clear()

    def _plot_trace(self):
        trace_lines = []
        for block in self.inst.blocks:
            trace_lines += block.trace_lines

        for line in trace_lines:
            self.ax_or.plot(
                np.add(line[0], self.offset[0]),
                np.add(line[2], self.offset[2]),
                np.add(line[1], self.offset[1]),
                'r-',
                linewidth=1
            )
            self.ax_zx.plot(
                np.add(line[2], self.offset[2]),
                np.add(line[0], self.offset[0]),
                'r-',
                linewidth=1
            )
            self.ax_zy.plot(
                np.add(line[2], self.offset[2]),
                np.add(line[1], self.offset[1]),
                'r-',
                linewidth=1
            )
            self.ax_xy.plot(
                np.add(line[0], self.offset[0]),
                np.add(line[1], self.offset[1]),
                'r-',
                linewidth=1
            )

    def _plot_kernel(self, kernel):
        s_lines = self.inst.blocks[kernel.block].components[kernel.comp_name]["scat_kernel"].lines()
        g_lines = self.inst.blocks[kernel.block].components[kernel.comp_name]["geom_kernel"].lines()

        from .mcramp import frame_rotate
        g_lines = frame_rotate(g_lines, kernel.rot)
        s_lines = frame_rotate(s_lines, kernel.rot)

        g_lines_xz = [
            np.add(g_lines[2], kernel.pos[2]+self.offset[2]),
            np.add(g_lines[0], kernel.pos[0]+self.offset[0])
        ]
        self.ax_zx.plot(g_lines_xz[0], g_lines_xz[1])

        g_lines_zy = [
            np.add(g_lines[2], kernel.pos[2]+self.offset[2]),
            np.add(g_lines[1], kernel.pos[1]+self.offset[1])
        ]
        self.ax_zy.plot(g_lines_zy[0], g_lines_zy[1])

        g_lines_xy = [
            np.add(g_lines[0], kernel.pos[0]+self.offset[0]),
            np.add(g_lines[1], kernel.pos[1]+self.offset[1])
        ]
        self.ax_xy.plot(g_lines_xy[0], g_lines_xy[1])

        g_lines_xyz = [
            np.add(g_lines[0], kernel.pos[0]+self.offset[0]),
            np.add(g_lines[1], kernel.pos[1]+self.offset[1]),
            np.add(g_lines[2], kernel.pos[2]+self.offset[2])
        ]
        self.ax_or.plot(g_lines_xyz[0],
                   g_lines_xyz[2],
                   g_lines_xyz[1])

        s_lines_xz = [
            np.add(s_lines[2], kernel.pos[2]+self.offset[2]),
            np.add(s_lines[0], kernel.pos[0]+self.offset[0])
        ]
        self.ax_zx.plot(s_lines_xz[0], s_lines_xz[1])

        s_lines_zy = [
            np.add(s_lines[2], kernel.pos[2]+self.offset[2]),
            np.add(s_lines[1], kernel.pos[1]+self.offset[1])
        ]
        self.ax_zy.plot(s_lines_zy[0], s_lines_zy[1])

        s_lines_xy = [
            np.add(s_lines[0], kernel.pos[0]+self.offset[0]),
            np.add(s_lines[1], kernel.pos[1]+self.offset[1])
        ]
        self.ax_xy.plot(s_lines_xy[0], s_lines_xy[1])

        s_lines_xyz = [
            np.add(s_lines[0], kernel.pos[0]+self.offset[0]),
            np.add(s_lines[1], kernel.pos[1]+self.offset[1]),
            np.add(s_lines[2], kernel.pos[2]+self.offset[2])
        ]
        self.ax_or.plot(s_lines_xyz[0],
                   s_lines_xyz[2],
                   s_lines_xyz[1])

        self._update_limits()

    def _plot_inst(self):
        for d in self.inst.kernel_refs:
            if not d.vis:
                continue

            self._plot_kernel(d)

        self._plot_trace()

    def show(self):
        self._plot_inst()
            
        plt.show()
