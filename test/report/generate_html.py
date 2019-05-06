#!/usr/bin/env python2

import os
from os.path import sep
import sys
import fnmatch
import collections
from pprint import pprint
import json
import subprocess
from shutil import copyfile
from itertools import izip
import time
import re
from multiprocessing import Pool
from collections import defaultdict

import argparse

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from dicts_for_html_page import *

epilog = """
The script will recursively look for PDF files in the directories you
specified.  For each PDF file it will copy the PDF to the output directory,
generate a PNG file with the same name (different extension, of course), and
add text to the index.html file in the output directory.  For every directory
and subdirectory the script will add headings to the HTML file.
"""

class OrderedDefaultDict(collections.OrderedDict):
    def __init__(self, *args, **kwargs):
        if not args:
            self.default_factory = None
        else:
            if not (args[0] is None or callable(args[0])):
                raise TypeError('first argument must be callable or None')
            self.default_factory = args[0]
            args = args[1:]
        super(OrderedDefaultDict, self).__init__(*args, **kwargs)

    def __missing__ (self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = default = self.default_factory()
        return default

    def __reduce__(self):  # optional, for pickle support
        args = (self.default_factory,) if self.default_factory else ()
        return self.__class__, args, None, None, self.iteritems()

def gen_heading(title, depth, id=''):
    """Generate HTML for a heading."""
    if id:
        id = 'id="{}"'.format(id)
    return '<h{depth} {id}>{title}</h{depth}>'.format(title=title, depth=depth, id=id)

def gen_text(text, depth):
    """Generate HTML for a comment regarding known deviations."""
    return '<h{depth} style="font-weight:normal;">{text}</h{depth}>'.format(text=text, depth=depth)

def gen_image(src, width='40%', alt='', link=None):
    """Generate HTML for an image."""
    img = '<img alt="{}" src="{}" style="width:{};" />' \
        .format(alt, src, width)
    if link:
        img = '<a href="{}">{}</a>'.format(link, img)
    return img

def gen_image_group(sources):
    """Generate HTML for a group of images, optionally with links."""
    if len(sources[0]) == 2:
        images = (gen_image(src, link=link) for src, link in sources)
    elif len(sources[0]) == 1:
        images = (gen_image(src) for src in sources)
    else:
        raise ValueError('expected list of image names or list of (image, link) tuples')
    return '<p style="padding-top: 0.5cm">{}</p>'.format(' '.join(images))

def gen_config_links(configs):
    """Generate HTML for links to the SMASH configs."""
    links = ('<a href="{}">Config {}</a>'.format(c, i)
             for i, c in enumerate(configs))
    return '\n'.join(links)

def gen_document(body):
    """Generate HTML document given a body."""
    return '''<!DOCTYPE html>
    <html>
    <head>
    <style>

    #toc {{
        list-style: none;
        position: fixed;
        right: 0;
        top: 0;
        height:auto !important;
        background-color: #FFFFFF;
        box-shadow: 0 0 1em #E6E6E6;
        padding: 0.3em;
        max-height: 100%;
        overflow-y: auto;
    }}
    #toc ul {{
        margin-right: 0em;
        margin-left: 0em;
    }}
    #toc li {{
        overflow: hidden;
    }}
    .list {{
        display: none;
    }}
    .hide {{
        display: none;
    }}
    .show:target + .hide {{
        display: inline;
    }}
    .show:target {{
        display: none;
    }}
    .show:target ~ .list {{
        display: inline;
    }}
    .flexbox-style {{
      display: flex;
      margin-left: auto;
      margin-right: 5%;
    }}
    .black {{
      color: black;
    }}

    @media print {{
        .hide, .show {{ display: none; }}
    }}

    </style>
    </head>
    <meta charset="utf-8"/>
    <body link="black", width = 100%>
    <font size="4" face="Helvetica" color="black">
    {}
    <p style="font-size: small">
    <a href="https://www.gsi.de/en/bottommenu/impressum.htm?no_cache=1&page=69"> Impressum</a>
    <a href="https://www.gsi.de/en/bottommenu/data_privacy_protection.htm"> Data privacy protection</a>
    </p>
    </font>
    </body>
    </html>
    '''.format(body)

def find_files(pattern, top):
    """Iterator over all files matching `pattern` in directory `top`."""
    for path, dirlist, filelist in os.walk(top):
        for name in fnmatch.filter(filelist, pattern):
            yield os.path.join(path, name)

def concatenate(sources):
    """Concatenate iterators."""
    for s in sources:
        for item in s:
            yield item

def print_all(iterator):
    """Print all elements in an iterator."""
    for i in iterator:
        print i

def Tree():
    """Create a tree from nested dictionaries."""
    return OrderedDefaultDict(Tree)

def make_id(heading, counter):
    """Make a unique ID given a heading.

    `counter` is a collections.Counter and used to make sure the headings are
    unique.
    """
    assert isinstance(counter, collections.Counter)
    counter[heading] += 1
    return heading + '-' + str(counter[heading])

def make_html(lines, tree, counter, heading_depth=1, target = 'none'):
    """Recursively generate HTML from tree of plots."""
    for k, v in tree.iteritems():
        if heading_depth == 2:
            lines.append('<font color="maroon">')
            lines.append(gen_heading(target_dict[k], heading_depth, make_id(k, counter)))
            lines.append('</font>')
            lines.append(gen_text(descriptions_targetpage[k], heading_depth+2))
            lines.append('<hr/>')
            target = k
            if target == 'FOPI_pions':
                lines.append(gen_text('Comments:    '
                             + comments[target], heading_depth+2))
        elif heading_depth ==1:
            lines.append(gen_heading(k, heading_depth, make_id(k, counter)))
        else:
            #lines.append('<font color="maroon">')
            lines.append(gen_heading(target_headers[target][k], heading_depth, make_id(k, counter)))
            #lines.append('</font>')
            lines.append(gen_text('Comments:    '
                         + comments[target][k], heading_depth+2))

        lines.append('')
        if type(v) is dict:
            lines.append(gen_image_group(v['plots']))
            lines.append(gen_config_links(v['configs']))
            lines.append('<hr/>')
            lines.append('<span style="display:block; height: 0.5cm;"></span>')
            lines.append('')
        else:
            make_html(lines, v, counter, heading_depth=(heading_depth + 1), target = k)

def make_toc(lines, tree, counter, target):
    """Recursively generate HTML table of contents from tree."""
    first = not bool(counter)
    if first:
        lines.append('')
        lines.append('''<div id="toc">
        <a href="#expand-toc" class="show" id="expand-toc">[expand TOC]</a>
        <a href="#collapse-toc" class="hide" id="collapse-toc">[collapse TOC]</a>
        <br>
        ''')
        classstring = ' class="list"'
    else:
        classstring = ''
    lines.append('<ul{}>'.format(classstring))
    for k, v in tree.iteritems():
        #lines.append('<li><a href="#{}">{}</a></li>'.format(make_id(k, counter), k))
        if target == 'FOPI_pions':
            lines.append('<li><a href="#' + k + '-1">' + target_dict[k] +
                         '</a></li>'.format(make_id(k, counter), k))
        else:
            if isinstance(v, OrderedDefaultDict):
                lines.append('<li><a href="#' + k + '-1">' + target_dict[k] +
                             '</a></li>'.format(make_id(k, counter), k))
                make_toc(lines, v, counter, k)
            else:
                lines.append('<li><a href="#' + k + '-1">' +
                             target_headers[target][k]+ '</a></li>'.format(
                             make_id(k, counter), k))
    lines.append('</ul>')
    if first:
        lines.append('</div>')

def _make_remove_string(remover, doc=None):
    """Factory constructing a string remover based on a regular expression."""

    def remove(s, to_be_removed):
        match = re.match(remover.format(re.escape(to_be_removed)), s)
        if match:
            return match.group(1)
        return s

    if doc is not None:
        remove.__doc__ = doc

    return remove

remove_string_left = _make_remove_string('^{}(.*)',
    doc='Remove the second string from the left of the first string.')
remove_string_right = _make_remove_string('(.*){}$',
    doc='Remove the second string from the right of the first string.')

def as_png(pdfpath):
    """Replace the `.pdf` ending in the given path with `.png`."""
    return remove_string_right(pdfpath, '.pdf') + '.png'

def strip_upper_dirs(path, source_dir):
    """Only use last directory in source path."""
    to_be_stripped = sep.join(remove_string_right(source_dir, sep)
                              .split(sep)[:-1]) + sep
    return remove_string_left(path, to_be_stripped)

def convert_to_png(plot, plotname):
    """Convert the given plot from a pdf to a png file."""
    if not isinstance(plotname, str):
        # then, we have a list which is only the case for the energy_scan,
        # where we want a directory substructure.
        plotdir= os.path.join(output_dir, plotname[0])
        plotpath = os.path.join(plotdir, plotname[1])
    else:
        plotpath = os.path.join(output_dir, plotname)
        plotdir = os.path.split(plotpath)[0]

    try:
        os.makedirs(plotdir)
    except OSError:
        # The directory was already created, so we don't have to do anything.
        # Note that using os.path.exists would result in a race condition.
        pass
    subprocess.check_call(['pdftoppm', '-png', '-scale-to', '1600',
                           plot, plotpath])
    os.rename(plotpath + '-1.png', plotpath + '.png')
    copyfile(plot, plotpath + '.pdf')
    return plotpath + '.png'

def get_basedir_filename(filename, source_dir, ending):
    """Extract the name of the base directory and generate the file name."""
    basedir, filename = os.path.split(remove_string_right(filename, ending))
    basedir = strip_upper_dirs(basedir, source_dir)
    filename = basedir.replace(sep, '-') + '-' + filename
    return basedir, filename

def remove_redundant_configs(configs):
    """Remove redundant configs from the same run.

    `configs` has to be sorted.
    """
    new_configs = []
    prev_config = None
    unify = lambda x: x.split('/data/')[0]
    for c in configs:
        if prev_config is not None:
            if unify(c) != unify(prev_config):
                new_configs.append(c)
        else:
            new_configs.append(c)
        prev_config = c
    return new_configs

def include_back_button(lines):
    """Add a back button on the HTML page."""
    lines.append(''' <a href=".."> <button >Back</button> </a> ''')
    return lines

def make_frontpage():
    """Generate the front page of the analysis suite."""
    title = 'Analysis Results'
    lines = []
    explanation_header = \
            "This webpage contains a collections of results that are produced \
             by means of the SMASH-analysis suite for each tagged version of \
             SMASH. It consists of basic tests, to ensure SMASH works as \
             expected, as well as physics results. The produced results are \
             grouped into different categories:"
    lines.append('<hr/>')
    lines.append('<div>')
    lines.append('<div style="width:50%; height:200; float:left;">')
    lines.append('<font color="maroon">')
    lines.append(gen_heading(title, 1))
    lines.append('</font>')
    lines.append(gen_heading('for version ' + smash_version, 2))
    lines.append('</div>')
    lines.append('<div style="width:50%; height:200; float:right;">')
    lines.append('<img src="smash_logo.png" alt="SMASH_Logo" height="200" class="flexbox-style">')
    lines.append('</div>')
    lines.append('</div>')
    lines.append('<p style="margin-bottom:1.0cm;"> created on {}'.format(time.strftime("%Y-%m-%d %H:%M")) + '</p>')
    lines.append(explanation_header)
    lines.append('<ul style="margin-top:0.5cm; margin-left: 2%;" >')

    # Sort folders alphabetically
    dictionary_list = os.listdir(args.source_dirs[0])
    dictionary_list.sort(key=lambda v: v.upper())
    for directory in dictionary_list:
        if directory[0] != ".":
            lines.append('''<b> <li style="padding-bottom:0.5cm;color:#800000">\
                <a style="text-decoration: none; color:#800000" href="'''
                         + str(directory) + '/">'
                         + target_dict[directory] + '</b>' +
                         '<br/> <span class="black">' + descriptions_frontpage[directory] + '</span>' + '</a>')
    path_to_logo = os.path.dirname(os.path.abspath(__file__)) + '/smash_logo.png'
    copyfile(path_to_logo, output_dir + '/smash_logo.png')
    lines.append('</ul>')
    lines.append('<hr/>')
    return lines

def make_energy_scan_page(tree):
    """Generate the energy scan subpage (similar to the front page)."""
    title = 'Analysis Suite for ' + smash_version
    lines = []
    include_back_button(lines)
    lines.append('<hr/>')
    lines.append(gen_heading(title, 1))
    lines.append('<font color="maroon">')
    lines.append(gen_heading('Energy Scan', 2))
    lines.append('</font>')
    lines.append(gen_text(descriptions_targetpage['energy_scan'], 4))
    lines.append('<hr/ style="margin-bottom:0.6cm;">')
    lines.append('<ul style="margin-top:0.8cm;">')

    for observable in energy_scan_sorted_observables:
            lines.append('''<b> <li style="padding-bottom:0.4cm;color:#800000">\
                <a style="text-decoration: none; color:#800000" href="'''
                         + str(observable) + '/">'
                         + E_scan_observable_dict[observable] + '</b> </a>')
    lines.append('</ul>')
    lines.append('<hr/>')
    lines.append(gen_config_links((tree['energy_scan']['configs'])))
    lines.append('<hr/>')
    lines.append('<span style="display:block; height: 0.5cm;"></span>')
    lines.append('')
    include_back_button(lines)

    return lines

def make_target_page(tree, target):
    """Generate the subpage of a specific target."""
    lines = []
    include_back_button(lines)
    lines.append('<hr/ >')
    counter = collections.Counter()
    # skip first heading
    subtree = tree.itervalues().next()
    make_toc(lines, subtree, counter, target)
    lines.append('')

    counter = collections.Counter()
    make_html(lines, tree, counter)

    include_back_button(lines)

    return lines

def walk_tree(tree, source_directory, with_configs, with_plots):
    """ Walk tree of directories to find plots and configs. """
    if with_plots:
        print '-> Looking for plots ...',
        sys.stdout.flush()
        plots = list(find_files('*.pdf', source_directory))
        plots.sort()
        print 'done.'

    if with_configs:
        print '-> Looking for configs...',
        configs = list(find_files('config.yaml', source_directory))
        configs.sort()
        configs = remove_redundant_configs(configs)
        print 'done.'

    if with_plots:
        # convert plots to png
        print '-> Converting pdf to png...',
        sys.stdout.flush()
        pool = Pool()
        results = []
        for plot in plots:
            basedir, plotname = get_basedir_filename(plot, source_dir, '.pdf')
            if 'energy_scan' in basedir:
                if len(plotname.split('-')) == 3:
                    plotname = plotname.split('-')[1] + '-' + \
                               plotname.split('-')[2]
                if len(plotname.split('-')) == 4:
                    #additional '-' because of negative charge
                    plotname = plotname.split('-')[1] + '-' + \
                               plotname.split('-')[2] + '-' + \
                               plotname.split('-')[3]
                # we want a folder substructure for the energy scan,
                # therefore return list with subdirectory and plotname as
                # plotname
                plotname = [basedir.split("/")[1], plotname]
            results.append(pool.apply_async(convert_to_png, (plot, plotname)))
        pool.close()
        # join the processes and raise exceptions
        for result in results:
            result.get()
        print 'done.'

    new_configs = []
    if with_configs:
        print '-> Copying configs ... ',
        # copy configs
        for c in configs:
            left, right = c.split('/data/')
            config = os.path.join(left, right)
            basedir, filename = get_basedir_filename(config, source_directory, '')
            new_configs.append(filename)
            configpath = os.path.join(output_dir, filename)
            try:
                os.makedirs(output_dir)
            except OSError:
                # The directory was already created, so we don't have to do anything.
                # Note that using os.path.exists would result in a race condition.
                pass
            copyfile(c, configpath)

            if not with_plots:
                dirs = basedir.split(sep)
                for d in dirs:
                    tree[d] = {
                        'plots': [], 'configs': new_configs,
                    }
        print 'done.'

    if with_plots:
        print '-> Build tree of plots ... ',
        for plot in plots:
            basedir, plotname = get_basedir_filename(plot, source_directory, '.pdf')
            # only use last directory of source for heading
            dirs = basedir.split(sep)
            prevtree = tree[title]
            for d in dirs[:-1]:
                prevtree = prevtree[d]
            lastdir = dirs[-1]

            # Filter only those configs that belong to the specific subsection.
            # For energy scan and FOPI pions, there is no filtering necessary
            # since multiple SMASH runs contribute to the plotted quantities.
            if 'energy_scan' in basedir or 'FOPI_pions' in basedir:
                relevant_configs = new_configs
            else:
                relevant_configs = []
                for config in new_configs:
                    if config.split("-")[-3] == basedir.split("/")[-1]:
                        relevant_configs.append(config)

            if lastdir not in prevtree:
                prevtree[lastdir] = {
                    'plots': [], 'configs': relevant_configs,
                }
            prevtree[lastdir]['plots'].append((plotname + '.png', plotname + '.pdf'))
        print 'done.'

    print
    if not tree:
        raise ValueError('-> No pdfs or configs found!')

    print '-> Generating HTML...',
    sys.stdout.flush()

    return tree

def reformat_energy_scan_tree(tree, observable):
    ''' Reformat tree of energy scan to group the multiplets. Additional layer
        required, one for each multiplet.
        For mtspectra and yspectra, we also want a sorting according to
        collision system and collision energy. '''
    new_tree = Tree()
    subtree = tree[title][observable]

    if observable not in ['mtspectra', 'yspectra']:
        for pdg in pdgs_sorted:
            for plot in subtree['plots']:
                if pdg in plot[1]:
                    new_tree[title][observable][pdgs_to_name[pdg]] = {
                        'plots': [plot], 'configs': [], }

    else:
        for pdg in pdgs_sorted:
            for spectra in energy_scan_sorted_spectra:
                for plot in subtree['plots']:
                    if pdg in plot[1]:
                        if spectra in plot[1]:
                            if len(new_tree[title][observable][pdgs_to_name[pdg]]) == 0:
                                new_tree[title][observable][pdgs_to_name[pdg]] = {
                                    'plots': [plot], 'configs': [], }
                            else:
                                new_tree[title][observable][pdgs_to_name[pdg]]['plots'].append(plot)

    return new_tree

def reorder_tree(tree, target, observable = 'none'):
    ''' Change the order in which the plots appear on the wiki page to
        a custom format.
        `energy_scan` is reordered separately by reformat_energy_scan_tree. '''

    if target in ['angular_distributions', 'detailed_balance', 'dileptons', \
                  'elastic_box']:
        # enforce custom ordering, by importance.
        new_tree = Tree()
        subtree = tree[title][target]

        for subtarget in sorted_subtargets[target]:
            for plot in subtree[subtarget]['plots']:
                if len(new_tree[title][target][subtarget]) == 0:
                    new_tree[title][target][subtarget] = {
                        'plots': [plot], 'configs': subtree[subtarget]['configs'], }
                else:
                    new_tree[title][target][subtarget]['plots'].append(plot)

    elif target in ['cross_sections']:
        # we want a more specific ordering, where the xs_process_type plots and
        # the xs_grouped plots are displayed first.
        new_tree = Tree()
        subtree = tree[title][target]

        for subtarget in sorted_subtargets[target]:
            if len(subtree[subtarget]['plots']) != 0:
                plot_list = []
                # We want a specific order for the most important plots, they
                # collected in the dictionary 'cross_sections_order', so we
                # loop over it's entries.
                for x_sec_plot_name in cross_sections_order:
                    for plot in subtree[subtarget]['plots']:
                        if x_sec_plot_name in plot[1]:
                            plot_list.append(plot)
                            # Remove those that are already part of the new tree
                            subtree[subtarget]['plots'].remove(plot)
                # Append the remaining plots, about whose order we do not care
                # anymore, to the new tree.
                for plot in subtree[subtarget]['plots']:
                    plot_list.append(plot)

                new_tree[title][target][subtarget] = {
                    'plots': plot_list, 'configs': subtree[subtarget]['configs'], }

    else:
        new_tree = tree

    return new_tree

def write_html(output_dir, lines):
    """ Write the generated html lines to the corresponsing index.html file. """
    with open(os.path.join(output_dir, 'index.html'), 'w') as f:
        s = gen_document('\n'.join(lines))
        f.write(s)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='generate HTML report of the anlaysis',
        epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('source_dirs', metavar='SOURCE', nargs='+',
        help='source directories that will be searched for plots')
    parser.add_argument('-o', '--output', nargs=1, default=['report_output'],
        help='output directory for report')
    parser.add_argument('--version', default = 'X.xx',
        help='SMASH version for which the wiki page is generated.')
    parser.add_argument('--overwrite', action='store_true',
        help='overwrite existing files')
    parser.add_argument('--frontpage', action = 'store_true')
    args = parser.parse_args()
    output_dir = os.path.abspath(args.output[0])

    smash_version = args.version

    if not args.frontpage:
        if (not args.overwrite) and os.path.exists(output_dir) and os.listdir(output_dir):
            print '-> "{}" is not empty, please clean it up or choose a different' \
                  '-> output directory!'.format(output_dir)
            sys.exit(1)

        # walk given directories and build tree of plots
        tree = Tree()
        title = 'Analysis Suite for ' + smash_version
        for source_dir in args.source_dirs:
            if 'energy_scan' in source_dir:
                # generate subpage for energy scan
                tree = walk_tree(tree, source_dir, with_configs = True, with_plots = False)
                lines = make_energy_scan_page(tree)
                write_html(output_dir, lines)

                # generate target pages for energy scan observables
                for observable in energy_scan_sorted_observables:
                    tree_e_scan = Tree()
                    source_dir_e_scan = os.path.join(source_dir, observable)
                    output_dir_e_scan = os.path.join(output_dir, observable)
                    tree_e_scan = walk_tree(tree_e_scan, source_dir_e_scan, with_configs = False, with_plots = True)
                    tree_e_scan = reformat_energy_scan_tree(tree_e_scan, observable)
                    lines = make_target_page(tree_e_scan, observable)
                    write_html(output_dir_e_scan, lines)
            else:
                tree = walk_tree(tree, source_dir, with_configs = True, with_plots = True)
                tree = reorder_tree(tree, source_dir.split('/')[-2])
                lines = make_target_page(tree, source_dir.split('/')[-2])
                write_html(output_dir, lines)

    else:
        lines = make_frontpage()
        write_html(output_dir, lines)

    print 'done.'
    print
