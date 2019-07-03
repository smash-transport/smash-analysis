import matplotlib
from distutils.version import LooseVersion
modern_matplotlib = (LooseVersion(matplotlib.__version__) >= LooseVersion("2.0.0"))
import matplotlib.pyplot as plt

import ConfigParser

import argparse

if modern_matplotlib:
    from cycler import cycler
from itertools import cycle

def get_default_colors(config):
    """Get the matplotlib default colors from given config."""
    if modern_matplotlib:
        return config['axes.prop_cycle'].by_key()['color']
    else:
        return config['axes.color_cycle']


def set_default_colors(config, default_colors):
    """Set the default colors of the given matplotlib config."""
    if modern_matplotlib:
        config['axes.prop_cycle'] = cycler('color', default_colors)
    else:
        config['axes.color_cycle'] = default_colors


def errorcontour(x, y, yerr, **kwargs):
    """Plot a line with a filled are around it, representin the error.

    Will use a given `axis` if provided as a keyword arguement.
    """
    fill_kwargs = {}
    if 'color' in kwargs:
        fill_kwargs['color'] = kwargs['color']
    else:
        default_colors = get_default_colors(plt.rcParams)
        kwargs['color'] = default_colors[0]
        fill_kwargs['color'] = default_colors[0]
    axis = kwargs.pop('axis', plt)
    axis.fill_between(x, y - yerr, y + yerr, alpha=0.4, **fill_kwargs)
    axis.plot(x, y, **kwargs)


default_colors = ['r', 'b', 'c', 'm', 'y', 'k', 'g']

default_line_styles = [
    # dash: 20 points on, 20 points off
    [20, 10],
    # small-dash: 10 points on, 5 points off
    [10, 5],
    # dash-dotted line: 8 on, 4 off, 2 on, 4 off
    [8, 4, 2, 4],
    # dash-dot-dot line
    [8, 4, 2, 4, 2, 4],
    # dot-dash-dash line
    [2, 4, 8, 4, 8, 4],
    # dotted line: 2 points on, 2 points off
    [2, 2]
]


smash_style_default = {
        'backend' : 'pdf',

        'axes.labelsize': 60,
        'axes.titlesize': 50,
        'axes.linewidth' : 2.0,
        'axes.unicode_minus': True,

        'legend.frameon': False,
        'legend.fontsize': 39,
        'legend.numpoints': 1,
        'legend.loc': 'best',

        'lines.markeredgewidth'  : 5.0,
        'lines.linewidth' : 5.0,
        'lines.markersize'  : 10,
        'lines.markeredgewidth'  : 2.5,

        'mathtext.fontset': 'stix',
        #'text.usetex': True,  - crashes
        #'font.family':'serif',
        'font.serif':['Computer Modern Roman', 'CMU Serif'],
        'font.monospace': ['Computer Modern Typewriter', 'CMU Typewriter Text'],
        'font.sans-serif': ['DejaVu Sans'],   # on Mac Os, Arial is needed
        'font.size' : 30,
        'font.family':'sans-serif',
        'mathtext.default':'rm',

        'xtick.labelsize': 40,
        'xtick.major.size':10,
        'xtick.minor.size':5,
        'xtick.major.pad':15,
        'xtick.major.width' : 2.0,

        'ytick.labelsize': 40,
        'ytick.major.size':10,
        'ytick.minor.size':5,
        'ytick.major.pad':15,
        'ytick.major.width' : 2.0,

        'figure.figsize': (20., 15.),
        'figure.autolayout': True,
}

set_default_colors(smash_style_default, default_colors)

class SmashStyle(argparse.Namespace):
    """Smash style plot template for matplotlib

    The smash style provides proper figsize, fontsize, fonttype, linewidth,
    color cycles, ticksize and padding and legend style.

    Example usage:
    If one want to use the basic styles except linestyles, increased title
    and plot space, minorticks, one can import the template file in the
    beginning of python script:

          from common_plotting import smash_style

    or simply:

          import common_plotting

    Otherwise, call the set function before plt.legend() if there is one
    or before plt.show() or plt.savefig() if there is no plt.legend() like:

          smash_style.set()
    """

    def __init__(self):
        self.params = smash_style_default
        plt.rcParams.update(self.params)

    def get_all_figs(self):
        """Get all figs in the current plot."""
        fignums = plt.get_fignums()
        figs = [plt.figure(i) for i in fignums]

        return figs

    def get_all_subplots(self, fig):
        """Get all subplots in one fig."""
        return fig.get_axes()

    def get_all_lines(self, axes):
        """Get all lines in one subplots."""
        children = axes.get_children()
        lines = [l for l in children if type(l)==matplotlib.lines.Line2D]
        return lines

    def update_legends(self, ax):
        """Get all legends in one ax."""
        children = ax.get_children()
        legends = [l for l in children if type(l)==matplotlib.legend.Legend]

        for legend in legends:
            title = legend.get_title().get_text()
            handles, labels = ax.get_legend_handles_labels()
            if title == 'None':
                ax.legend(handles, labels)
            else:
                ax.legend(handles, labels, title=title)

    def set(self, line_styles=True, title_padding=1.03, minorticks_on=True,
            update_legends=False):
        """Apply the smash style to the current script.

        Args:
             line_styles (bool, optional): True to use smash line style cycles.
                 Defaults to True.
             title_padding (float, optional): To increase padding between title
                 and plot. Defaults to 1.03, bigger title_padding to move title
                 upper.
             minorticks_on: True to turn on minor ticks. Defaults to True.
                 False to switch off minorticks
        """
        figs = self.get_all_figs()
        for fig in figs:
            axes = self.get_all_subplots(fig)
            for ax in axes:
                if minorticks_on:
                    ax.minorticks_on()

                if title_padding:
                    #increase the title padding
                    title_text = ax.get_title()
                    ax.set_title(title_text, y=title_padding)

                if line_styles:
                    lines = self.get_all_lines(ax)

                    line_style_cycle = cycle(default_line_styles)
                    for i, line in enumerate(lines):
                        #skip if line.markers=None
                        if line.get_linestyle() == 'None':
                            continue
                        #for i==0, use solid line
                        if i > 0:
                            plt.setp(line, dashes=next(line_style_cycle))

                    if update_legends:
                        self.update_legends(ax)

smash_style = SmashStyle()
