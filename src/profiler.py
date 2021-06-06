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
colors = ['r', 'b', 'y', 'g', 'm', 'c']

class CmdLineParser(object):
	def __init__(self):
		self.parser = OptionParser(usage='usage: python perfprof.py [options] cvsfile.csv outputfile.pdf')
		# default options
		self.parser.add_option("-D", "--delimiter", dest="delimiter", default=None, help="delimiter for input files")
		self.parser.add_option("-M", "--maxratio", dest="maxratio", default=4, type=int, help="maxratio for perf. profile")
		self.parser.add_option("-S", "--shift", dest="shift", default=0, type=float, help="shift for data")
		self.parser.add_option("-L", "--logplot", dest="logplot", action="store_true", default=False, help="log scale for x")
		self.parser.add_option("-T", "--timelimit", dest="timelimit", default=1e99, type=float, help="time limit for runs")
		self.parser.add_option("-P", "--plot-title", dest="plottitle", default=None, help="plot title")
		self.parser.add_option("-X", "--x-label", dest="xlabel", default='Time Ratio', help="x axis label")
		self.parser.add_option("-B", "--bw", dest="bw", action="store_true", default=False, help="plot B/W")

	def addOption(self, *args, **kwargs):
		self.parser.add_option(*args, **kwargs)

	def parseArgs(self):
		(options, args) = self.parser.parse_args()
		options.input = args[0]
		# options.output = args[1]
		return options


def readTable(fp, delimiter):
    """
    read a CSV file with performance profile specification
    the format is as follows:
    ncols algo1 algo2 ...
    nome_istanza tempo(algo1) tempo(algo2) ...
    ...
    """
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

    fig, ax1 = plt.subplots()

    # cut_table = nrows
    cut_table = nrows

    print("cols: ", cnames)
    print(sum(data[:cut_table,cnames['infeasible_dir']]), "infeasible out of", nrows)

    color = 'tab:orange'
    col = 'dual'
    ax1.set_xlabel('iteration')
    ax1.set_ylabel(col, color=color)
    ax1.plot(data[:cut_table, cnames[col]], color=color)

    # color = 'tab:blue'
    # col = 'infeasible_dir'
    # ax2 = ax1.twinx()
    # ax2.set_ylabel(col, color=color)
    # ax2.plot(data[:cut_table, cnames[col]], color=color)

    fig.tight_layout()
    plt.show()

if __name__ == '__main__':
	main()
