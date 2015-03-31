#!/usr/bin/env python
#-*- coding: utf-8 -*-

## Import common moduli
import matplotlib, sys, os, time
import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy.constants import c, hbar, pi

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--paramname', type=str, help='parameter by which the lines are sorted')
parser.add_argument('--paramunit',  type=float, default=1., help='prescaling of the parameter')
parser.add_argument('--xunit',  type=float, default=1., help='prescaling of the x-axis')
parser.add_argument('--yunit',  type=float, default=1., help='prescaling of the y-axis')
parser.add_argument('--paramlabel', type=str, default='parameter = %g', help='line label (use %d, %s, %f or %g to format the parameter, use LaTeX for typesetting)')
parser.add_argument('--xcol', type=str, default='0', help='number or exact name of the x-axis column') ## TODO or -- if it is to be generated
parser.add_argument('--ycol', type=str, default='1', help='number or exact name of the y-axis column')
parser.add_argument('--xlabel', type=str, default='', help='label of the x-axis (use LaTeX)')
parser.add_argument('--ylabel', type=str, default='', help='label of the y-axis (use LaTeX)')
parser.add_argument('--output', type=str, default='output.png', help='output file (e.g. output.png or output.pdf)')
parser.add_argument('--colormap', type=str, default='hsv', help='matplotlib colormap, available are: hsv (default), jet, gist_earth, greys, dark2, brg...')
parser.add_argument('filenames', type=str, nargs='+', help='CSV files to be processed')
args = parser.parse_args()


## Options 
#cmap = matplotlib.cm.gist_earth
cmap = getattr(matplotlib.cm, args.colormap)


## Use LaTeX
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', size=12)
matplotlib.rc('text.latex', preamble = '\usepackage{amsmath}, \usepackage{yfonts}, \usepackage{txfonts}') #, \usepackage{palatino} \usepackage{lmodern}, 
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman, Times']})  ## select fonts

def get_param(filename):             ## Load header to the 'parameters' dictionary#{{{
    parameters = {}
    with open(filename) as datafile:
        for line in datafile:
            if (line[0:1] in '0123456789') or ('column' in line.lower()): break    # end of parameter list
            key, value = line.replace(',', ' ').split()[-2:]
            try: value = float(value) ## Try to convert to float, if possible
            except: pass                ## otherwise keep as string
            parameters[key] = value
    return parameters
#}}}

## Start figure + subplot 
fig = plt.figure(figsize=(10,10))
fig.subplots_adjust(left=.05, bottom=.05, right=.99, top=.99, wspace=.05, hspace=.05) ## (for interactive mode)

## Sort arguments by a _numerical_ value in their parameter, keep the color order
print args.filenames
filenames = args.filenames

params  = [get_param(n)[args.paramname] for n in filenames]
datasets = zip(params, filenames)                               ## sort the files by the parameter
datasets.sort()
colors = cmap(np.linspace(0.0,0.9,len(filenames)+1)[:-1]) ## add the colors to sorted files
datasets = zip(colors, *zip(*datasets))

ax = plt.subplot(111, axisbg='w')

def loadtxt_columns(filename): #{{{
    columns     = []
    with open(filename) as datafile:
        for line in datafile:
            if ('column' in line.lower()): columns.append(line.strip().split(' ', 1)[-1]) # (todo) this may need fixing to avoid collision
    return columns
#}}}
def get_col_index(col, fn):#{{{
    columnnames = loadtxt_columns(fn)
    try:
        return int(col), columnnames[int(col)]      ## column number given, find its name
    except ValueError:
        try:
            return loadtxt_columns(fn).index(col), col      ## column name given, find its number
        except:
            raise ValueError, "Could not find column %s for the x-axis in file %s" % (col, fn)
#}}}

for color, param, filename in datasets:
    xcol, xcolname = get_col_index(args.xcol, filename)
    ycol, ycolname = get_col_index(args.ycol, filename)
    (x, y) = np.loadtxt(filename, usecols=[xcol, ycol], unpack=True)
    plt.plot(x/args.xunit, y/args.yunit, color=color, label=args.paramlabel % (param/args.paramunit))
plt.xlabel(xcolname if args.xlabel == '' else args.xlabel) 
plt.ylabel(ycolname if args.ylabel == '' else args.ylabel)
plt.grid()
plt.legend(prop={'size':12}, loc='upper left').draw_frame(False)


## ==== Outputting ====
plt.savefig(args.output, bbox_inches='tight')
