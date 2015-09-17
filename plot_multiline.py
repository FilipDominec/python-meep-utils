#!/usr/bin/env python
#-*- coding: utf-8 -*-

## Import common moduli
import matplotlib, sys, os, time
import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy.constants import c, hbar, pi

parser = argparse.ArgumentParser(description='Plot a selected column from each of multiple files into a single comparison plot. ')
parser.add_argument('--paramname',  type=str,               help='parameter by which the lines are sorted')
parser.add_argument('--paramunit',  type=float, default=1., help='prescaling of the parameter (if it is a number)')
parser.add_argument('--title',      type=str, default='', help='plot title')
# todo: remove the --xunit --yunit
parser.add_argument('--yunit',      type=float, default=1., help='prescaling of the y-axis')
parser.add_argument('--paramlabel', type=str,   default='', help='line label (use standard "printf percent substitutes" to format the parameter, use LaTeX for typesetting)')
parser.add_argument('--xcol',       type=str,   default='0', help='number or exact name of the x-axis column') ## TODO or -- if it is to be generated
parser.add_argument('--ycol',       type=str,   default='1', help='number or exact name of the y-axis column')
parser.add_argument('--xeval',      type=str,   default='x', help='any python expression to preprocess the `x`-values, e.g. `1e6*c/x` to convert Hertz to micrometers') 
parser.add_argument('--yeval',      type=str,   default='y', help='any python expression to preprocess the `y`-values, e.g. `y/x` to normalize against newly computed x') 
parser.add_argument('--xlim1',       type=str,   default='', help='start for the x-axis range')
parser.add_argument('--xlim2',       type=str,   default='', help='end for the x-axis range')
parser.add_argument('--ylim1',       type=str,   default='', help='start for the y-axis range')
parser.add_argument('--ylim2',       type=str,   default='', help='end for the y-axis range')
parser.add_argument('--xlabel',     type=str,   default='', help='label of the x-axis (can use LaTeX)')
parser.add_argument('--ylabel',     type=str,   default='', help='label of the y-axis (use LaTeX)')
parser.add_argument('--output',     type=str,   default='output.png', help='output file (e.g. output.png or output.pdf)')
parser.add_argument('--colormap',   type=str,   default='hsv', help='matplotlib colormap, available are: hsv (default), jet, gist_earth, greys, dark2, brg...')
parser.add_argument('filenames',    type=str,   nargs='+', help='CSV files to be processed')
args = parser.parse_args()


## Options 
#cmap = matplotlib.cm.gist_earth
cmap = getattr(matplotlib.cm, args.colormap)


## Use LaTeX
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', size=12)
matplotlib.rc('text.latex', preamble = '\usepackage{amsmath}, \usepackage{yfonts}, \usepackage{txfonts}, \usepackage{upgreek}') #, \usepackage{palatino} \usepackage{lmodern}, 
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman, Times']})  ## select fonts

def get_param(filename):             ## Load header to the 'parameters' dictionary#{{{
    parameters = {}
    with open(filename) as datafile:
        for line in datafile:
            if (line[0:1] in '0123456789') or ('column' in line.lower()): break    # end of parameter list
            ## key-value separator is either ',' or '='; take the left word from it as the param name, and everything on the right as the param value
            left, value = line.replace(',', '=').strip().split('=', 1)
            key = left.split()[-1]
            try: value = float(value) ## Try to convert to float, if possible
            except: pass                ## otherwise keep as string
            parameters[key] = value
    return parameters
#}}}

## Start figure + subplot 
fig = plt.figure(figsize=(8,4))
fig.subplots_adjust(left=.05, bottom=.05, right=.99, top=.99, wspace=.05, hspace=.05) ## (for interactive mode)

def loadtxt_columns(filename): #{{{
    """ Retrieves the names of columns from a CSV file header """
    columns     = []
    with open(filename) as datafile:
        for line in datafile:
            if ('column' in line.lower()): columns.append(line.strip().split(' ', 1)[-1]) # (todo) this may need fixing to avoid collision
    return columns
#}}}
def get_col_index(col, fn):#{{{
    columnnames = loadtxt_columns(fn)
    if columnnames == []: return int(col), ("column %s" % str(col))
    try:
        return int(col), columnnames[int(col)]      ## column number given, find its name
    except ValueError:
        try:
            return loadtxt_columns(fn).index(col), col      ## column name given, find its number
        except:
            raise ValueError, "Could not find column %s for the x-axis in file %s" % (col, fn)
#}}}


## Sort arguments by a _numerical_ value in their parameter, keep the color order
filenames = args.filenames
params  = [get_param(n)[args.paramname] for n in filenames]
datasets = zip(params, filenames)                               ## sort the files by the parameter
datasets.sort()
colors = cmap(np.linspace(0.0,0.9,len(filenames)+1)[:-1]) ## add the colors to sorted files
datasets = zip(colors, *zip(*datasets))

ax = plt.subplot(111, axisbg='w')

## Cycle through all files
for color, param, filename in datasets:
    # identify the x,y columns by its number or by its name, load them and optionally process them with an expression
    xcol, xcolname = get_col_index(args.xcol, filename) 
    ycol, ycolname = get_col_index(args.ycol, filename)
    (x, y) = np.loadtxt(filename, usecols=[xcol, ycol], unpack=True)
    x = eval(args.xeval)
    y = eval(args.yeval)

    # if the legend format is not supplied by user, generate it from the parameter name 
    if type(param) in (float, int):
        label = (args.paramlabel % (param/args.paramunit)) if args.paramlabel else ("%s = %.3g" % (args.paramname, (param/args.paramunit)))
    else:
        label = (args.paramlabel % (param)) if args.paramlabel else ("%s = %s" % (args.paramname, param))
    plt.plot(x, y, color=color, label=label)

if args.xlim1 != "": plt.xlim(left=float(args.xlim1))
if args.xlim2 != "": plt.xlim(right=float(args.xlim2))
if args.ylim1 != "": plt.ylim(left=float(args.ylim1))
if args.ylim2 != "": plt.ylim(right=float(args.ylim2))
if args.title: plt.title(args.title)
plt.xlabel(xcolname if args.xlabel == '' else args.xlabel) 
plt.ylabel(ycolname if args.ylabel == '' else args.ylabel)
plt.grid()
plt.legend(prop={'size':12}, loc='upper left').draw_frame(False)

## ==== Outputting ====
plt.savefig(args.output, bbox_inches='tight')
