#ty prof

from __future__ import print_function
import numpy as np
import matplotlib
# matplotlib.use('PDF')
import matplotlib.pyplot as plt
import sys

from optparse import OptionParser

# parameters
defLW = 1.2  # default line width
defMS = 7  # default marker size
dashes = ['-',  # solid line
    '--',  # dashed line
    '-.',  # dash-dot line
    ':',  # dotted line
    '-',
    '--']

markers = ['+', 'x', 's', '^', 'o', 'd']
colors = ['tab:orange', 'tab:green', 'tab:blue', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

class CmdLineParser(object):
    def __init__(self):
        self.parser = OptionParser(usage='usage: python perfprof.py [options] cvsfile.csv outputfile.pdf')
        # default options
        # self.parser.add_option("-D", "--delimiter", dest="delimiter", default=None, help="delimiter for input files")
        # self.parser.add_option("-M", "--maxratio", dest="maxratio", default=4, type=int, help="maxratio for perf. profile")
        # self.parser.add_option("-S", "--shift", dest="shift", default=0, type=float, help="shift for data")
        # self.parser.add_option("-L", "--logplot", dest="logplot", action="store_true", default=False, help="log scale for x")
        # self.parser.add_option("-T", "--timelimit", dest="timelimit", default=1e99, type=float, help="time limit for runs")
        # self.parser.add_option("-P", "--plot-title", dest="plottitle", default=None, help="plot title")
        # self.parser.add_option("-X", "--x-label", dest="xlabel", default='Time Ratio', help="x axis label")
        # self.parser.add_option("-B", "--bw", dest="bw", action="store_true", default=False, help="plot B/W")
        self.parser.add_option("-f", "--fields", dest="fields", default=None, help="additional fields")
        self.parser.add_option("-q", "--quiet", action="store_true", dest="noplot", default=False, help="don't plot")

    def addOption(self, *args, **kwargs):
        self.parser.add_option(*args, **kwargs)

    def parseArgs(self):
        (options, args) = self.parser.parse_args()
        options.input = args[0]
        # options.output = args[1]
        return options


def readTable(fp, delimiter):
    firstline = fp.readline().strip().split(delimiter)
    ncols = len(firstline)

    cnames = {}
    for i, cname in enumerate(firstline):
        cnames[cname] = i

    rows = []
    for row in fp:
        row = row.strip().split(delimiter)
        rdata = np.empty(ncols)
        for j in range(ncols):
            rdata[j] = float(row[j])
        rows.append(rdata)
    data = np.array(rows)
    return (cnames, data)


def main():
    parser = CmdLineParser()
    opt = parser.parseArgs()
    cnames, data = readTable(open(opt.input, 'r'), ',')
    nrows, ncols = data.shape

    print("plotting", opt.input)

    print("avaiable fields:")
    for f in cnames.keys(): print("\t-", f)
    print()

    additional_fields = [] if opt.fields  is None else opt.fields.split(',')
    print("additional fields:")
    for f in additional_fields: print("\t-", f)
    print()

    for f in additional_fields:
        if f not in cnames.keys():
            print("error: field", f, "not avaiable")
            exit()

    print("some stats: ")
    print("\t-", nrows + 1, "iterations")
    print("\t-", sum(data[:,cnames['conditional']]), "infeasible direction out of", nrows + 1)

    fig = plt.figure()
    nplots = len(additional_fields)+1
    gs = fig.add_gridspec(nplots, hspace=0)

    axs = gs.subplots(sharex=True, sharey=False)
    title = opt.input + ": plotting "
    title += ','.join(["phi"] + additional_fields)
    title += "(phistar = " + str(max(data[:, cnames["dual"]])) + ")"
    fig.suptitle(title)

    x = range(nrows)
    phis = data[:, cnames["dual"]] 

    monotone_phis = [0] * len(phis)
    monotone_phis[0] = phis[0]
    for i in range(1, len(phis)):
        monotone_phis[i] = max(monotone_phis[i-1], phis[i])

    monotonicity_points = []
    for i in range(len(monotone_phis)):
        if i == 0:
            continue
        if monotone_phis[i] != monotone_phis[i-1]:
            monotonicity_points.append(i)

    phiaxes = None
    if nplots == 1:
        phiaxes = axs
    else:
        phiaxes = axs[0]

    color_idx = 0
    phiaxes.plot(x, monotone_phis, alpha=0.7, color=colors[color_idx], marker='o', markevery=monotonicity_points)
    color_idx += 1
    phiaxes.plot(x, phis, color=colors[color_idx])
    color_idx += 1

    for i, col in enumerate(additional_fields):
        axs[i + 1].plot(x, data[:, cnames[col]], color=colors[color_idx])
        color_idx += 1

    if nplots > 1:
        for ax in axs:
            ax.label_outer()

    # fig, ax1 = plt.subplots()
    #
    #
    # color = 'tab:orange'
    # ax1.set_xlabel('iteration')
    # ax1.set_ylabel("phi", color=color)
    # ax1.plot(monotone_phis, alpha=0.7, color='tab:green', marker='o', markevery=monotonicity_points)
    # ax1.plot(phis, color=color)

    # color = 'tab:blue'
    # col = 'infeasible_dir'
    # ax2 = ax1.twinx()
    # ax2.set_ylabel(col, color=color)
    # ax2.plot(data[:cut_table, cnames[col]], color=color)

    fig.tight_layout()
    if not opt.noplot: plt.show()

if __name__ == '__main__':
    main()
