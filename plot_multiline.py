#!/usr/bin/env python
#-*- coding: utf-8 -*-
# 
"""
plot_multiline.py  --  a script for ready-to-publish presentation of multiple data files in one plot
                                                                                                       
The data, stored as the `y' value depending on the `x' coordinate in each file, can be either plotted 
as multiple lines in the x-y plot, which is the default behaviour. Each curve is distinguished by the different
value of the `param', or simply by the file name if no --paramname option is specified.

Or, if the number of files is over 10 or 20, it may be preferable to plot their `y' values as a 2-D contour plot,
where the `y' value still represents the horizontal coordinate, but the `param' value now serves as the vertical one. 
The `--xeval' or `--parameval' expressions should return similar order of magnitude so that the data can be triangulated.

        1) --contours no              2) --contours yes
               (default)
           ^                             +---------------------+                                       
        y  |                          p  |                     |                                       
           |                          a  |                     |                                       
           |                          r  |          y          |                                       
           |                          a  |                     |                                       
           |                          m  |                     |                                       
           +------------------>          +---------------------+                                       
                           x                               x                                           

example use:
    ./plot_multiline.py *.dat    --paramname depth[nm] --contours yes  \\
                    --yeval 'np.abs(y)**2' --ylim1 0 \\
                    --xlabel 'Normalized grating period $\Lambda/\lambda$' --xeval x*1e-9/0.8e-6 --xlim1=0    \\
                    --paramlabel 'Grating depth $h$ ($\muup$m)' --parameval param*1e-9/1e-6  \\
                    --title "Reflectance $|r|^2$ of the grating at $\\\\lambda$=800 nm)" \\
                    --usetex yes --colormap 'Paired/2' \\
                    --output results_R_${PWD##*/}.png    

"""

## Import common moduli
import matplotlib, sys, os, time
import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy.constants import c, hbar, pi

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--paramname',  type=str, default='', help='parameter by which the lines are sorted (filename used if omitted)')
parser.add_argument('--title',      type=str, default='', help='plot title')
# todo: remove the --xunit --yunit
parser.add_argument('--yunit',      type=float, default=1., help='prescaling of the y-axis')
parser.add_argument('--paramlabel', type=str,   default='', help='line label; use standard printf percent substitutes to format the parameter, LaTeX for typesetting')
parser.add_argument('--xcol',       type=str,   default='0', help='number or exact name of the x-axis column') ## TODO or -- if it is to be generated
parser.add_argument('--ycol',       type=str,   default='1', help='number or exact name of the y-axis column')
parser.add_argument('--ycol2',      type=str,   default='', help='number or exact name of the auxiliary y-axis column; value can be accessed as `ycol2`; is plotted as a dashed line')
parser.add_argument('--xeval',      type=str,   default='x', help='Python expression to preprocess the `x`-values, e.g. `1e6*c/x` to convert Hertz to the wavelength in micrometers') 
parser.add_argument('--yeval',      type=str,   default='y', help='Python expression to preprocess the `y`-values, e.g. `y/x` to normalize against newly computed x') 
parser.add_argument('--y2eval',     type=str,   default='y2', help='Python expression to preprocess the auxiliary `y`-values (computed before `y` is processed)') 
parser.add_argument('--parameval',  type=str,   default='param', help='any python expression to preprocess the `param`-values, e.g. `param/1e-9` to convert it to nanometers') 
parser.add_argument('--xlim1',      type=str,   default='', help='start for the x-axis range')
parser.add_argument('--xlim2',      type=str,   default='', help='end for the x-axis range')
parser.add_argument('--ylim1',      type=str,   default='', help='start for the plotted value range')
parser.add_argument('--ylim2',      type=str,   default='', help='end for the plotteld value range')
parser.add_argument('--plim1',      type=str,   default='', help='start for the plotted parameter range')
parser.add_argument('--plim2',      type=str,   default='', help='end for the plotted parameter range')
parser.add_argument('--xlabel',     type=str,   default='', help='label of the x-axis (can use LaTeX)')
parser.add_argument('--ylabel',     type=str,   default='', help='label of the y-axis (use LaTeX)')
parser.add_argument('--output',     type=str,   default='output.png', help='output file (e.g. output.png or output.pdf)')
parser.add_argument('--colormap',   type=str,   default='default', help='matplotlib colormap, available are: hsv (default for lines), gist_earth (default for contours), jet, greys, dark2, brg...')
parser.add_argument('--overlayplot',type=str,   default='', help='one or more expressions, separated by comma, that are plotted to help guide the eye (e.g. 1/x)')
parser.add_argument('--numcontours',type=int,   default=50, help='number of levels in the contour plot (default 50)')
parser.add_argument('--contourresx',type=int,   default=200,help='row length of the internal interpolation matrix for contour plot (default 200)')
parser.add_argument('--contourresp',type=int,   default=200,help='column height of the internal interpolation matrix for contour plot (default 200)')

parser.add_argument('--figsizex',   type=float, default=8, help='figure width (inches), 8 is default')
parser.add_argument('--figsizey',   type=float, default=4, help='figure height (inches), 4 is default')
parser.add_argument('--contours',   type=str,   default='no', help='make a 2-D contour plot instead of multiple curves')
parser.add_argument('--usetex',     type=str,   default='yes', help='by default, LaTeX is used for nicer typesetting')
parser.add_argument('--verbose',    type=str,   default='', help='explicitly print out what happens')
parser.add_argument('filenames',    type=str,   nargs='+', help='CSV files to be processed')
args = parser.parse_args()

## Plotting style
if args.usetex.lower() in ('yes', 'true'): 
    if args.verbose: print "Selecting to use Latex, disable it by setting   --usetex no"
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('text.latex', preamble = '\usepackage{amsmath}, \usepackage{txfonts}, \usepackage{upgreek}') #, \usepackage{palatino} \usepackage{lmodern}, 
matplotlib.rc('font', size=12)
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman, Times']})  ## select fonts

if args.colormap == 'default': 
    cmap = matplotlib.cm.gist_earth if (args.contours=='yes') else matplotlib.cm.hsv
else:
    if (args.colormap[-2:] == '/2'):        ## allow palette halving,  (c) unutbu, stackoverflow:18926031
        cmap = getattr(matplotlib.cm, args.colormap[:-2]) 
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=.0, b=.47), cmap(np.linspace(.0, .47, 100)))
    else: 
        cmap = getattr(matplotlib.cm, args.colormap)  

def get_param(filename):             ## Load header to the 'parameters' dictionary#{{{
    parameters = {}
    with open(filename) as datafile:
        for line in datafile:
            try:
                if (line[0:1] in '0123456789.') or ('column' in line.lower()): break    # end of parameter list
                ## key-value separator is either ',' or '='; take the left word from it as the param name, 
                ## and everything on the right as the param value
                left, value = line.replace(',', '=').strip().split('=', 1)
                key = left.split()[-1]
                try: value = float(value) ## Try to convert to float, if possible
                except: pass                ## otherwise keep as string
                parameters[key] = value
            except:
                pass
    return parameters
#}}}
def reasonable_ticks(lim1, lim2, density=1, extend_to_lims=False): #{{{
    """ 
    Aims to improve the numbering of the axis or colormap in matplotlib plots.

    By default, about 7 to 14 human-friendly intermediate values are generated between lim1 and lim2. 
    Their mantissa increases by 1, 2, or 5 as usual in hand-made plots. Their total number can be adjusted by 
    the `density' parameter.

    If the value of lim1 or lim2 is nice enough, it is included in the range, too. 
    This behaviour can be forced by setting `extend_to_lims' is set to True, in which case the first and last 
    ticks are just set to the limits.

    >>> reasonable_ticks(2.4242, 4.4242)
    array([ 2.6,  2.8,  3. ,  3.2,  3.4,  3.6,  3.8,  4. ,  4.2,  4.4,  4.6])

    >>> reasonable_ticks(2.4242, 4.4242, density=.6)
    array([ 2.5,  3. ,  3.5,  4. ])

    >>> reasonable_ticks(2.4242, 4.4242, density=.6, extend_to_lims=True)
    array([ 2.4242,  3.    ,  3.5   ,  4.4242])
    """
    diff = float(lim2-lim1)
    decimal = 10**np.floor(np.log10(diff/(density*4)))                                 # get the order of magnitude 
    step = (decimal, 2*decimal, 5*decimal)[np.int(3*diff/(density*4)/decimal/10)]      # select the correct step
    newrange = np.arange(np.ceil(lim1/step-1e-6)*step, (np.ceil(lim2/step+1e-6))*step, step)         # generate the tick positions, including end points if rounded
    if extend_to_lims: newrange = np.hstack([lim1, newrange[1:-1], lim2])                          # optionally, force including the end points
    return newrange
#}}}

## Start figure + subplot 
fig = plt.figure(figsize=(args.figsizex,args.figsizey))
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
            return columnnames.index(col), col      ## column name given, find its number
        except:
            raise ValueError, "Could not find column '%s' for the x-axis in file %s;\n\tIndex it by number 0-%d or by name: '%s'" % (col, fn, len(columnnames), "', '".join(columnnames))
#}}}


## Sort arguments by a the value of the specified parameter, keep the color order
filenames = args.filenames
if args.paramname == '':    
    if args.usetex:
        ## plot underscores correctly in the file names; otherwise LaTeX complains
        params, paramname = [fn.replace('_', '\_') for fn in filenames],       'file name'
        print params
    else:
        params, paramname = filenames,                                         'file name'
else:                       
        params, paramname = [get_param(n)[args.paramname] for n in filenames], args.paramname
datasets = zip(params, filenames)                               ## sort the files by the parameter
datasets.sort()
colors = cmap(np.linspace(0.0,0.9,len(filenames)+1)[:-1]) ## add the colors to sorted files
datasets = zip(colors, *zip(*datasets))

ax = plt.subplot(111, axisbg='w')

## Cycle through all files
if args.contours == 'yes': 
    xs, ys, params  = [np.array([]) for _ in range(3)] ## three empty arrays
for color, param, filename in datasets:
    # identify the x,y columns by its number or by its name, load them and optionally process them with an expression
    xcol,  xcolname  = get_col_index(args.xcol,  filename) 
    ycol,  ycolname  = get_col_index(args.ycol,  filename)
    ycol2, ycolname2 = get_col_index(args.ycol2, filename) if args.ycol2 else (None,None)
    (x, y) = np.loadtxt(filename, usecols=[xcol, ycol], unpack=True)
    if args.verbose:
        print "Loading columns '%s' (%d), '%s' (%d) from the file %s" % (xcolname, xcol, ycolname, ycol, filename)
    x  = eval(args.xeval)
    if args.ycol2: 
        y2 = np.loadtxt(filename, usecols=[ycol2], unpack=True)
        y2 = eval(args.y2eval)
    y  = eval(args.yeval)

    # if the legend format is not supplied by user, generate it from the parameter name 
    if type(param) in (float, int):
        param = eval(args.parameval)
    else:
        if args.contours == 'yes': 
            raise ValueError("Parameter must be a number for contour plot, since it is used at the vertical axis")

    if not args.contours == 'yes':
        ## Plot a curve with a nice label, generated from the parameter
        if args.paramlabel.lower() == 'none':
            label = ''
        elif ("'" in args.paramlabel) or ('"' in args.paramlabel):
            label = args.paramlabel.strip('"').strip("'")
        elif args.paramlabel:
            if '%' in args.paramlabel:
                label = (args.paramlabel % param)           # manually formatted label
            else:
                label = (args.paramlabel+(" = %s" % param))           # manual label with parameter value appended
        else:
            label = ("%s = %s" % (args.paramname, param))   # automatic formatted label

        if args.ycol2: plt.plot(x, y2, color=color, label='', marker='s', markersize=(3 if len(x)<50 else 0), ls='--')
        plt.plot(x, y, color=color, label=label, marker='o', markersize=(3 if len(x)<50 else 0))
    else:
        ## Store the points for later interpolation and contour plot
        xs      = np.append(xs, x)
        ys      = np.append(ys, y)
        params  = np.append(params, np.ones_like(x)*param) # (fixme: fails if param is string)

if args.contours == 'yes':
    # Grid the data, produce interpolated quantities:
    from matplotlib.mlab import griddata
    xi      = np.linspace(min(xs),       max(xs),       args.contourresx)
    paramsi = np.linspace(min(params),   max(params),   args.contourresp)
    interp_anisotropy = 1       # value lower than 1. interpolates rather vertically; optimize if plot disintegrates
    yi      = griddata(xs, params*interp_anisotropy, ys, xi, paramsi*interp_anisotropy, interp='linear')

    # Standard contour plot
    cmaprange1 = float(args.ylim1) if (args.ylim1 != "") else np.min(yi) 
    cmaprange2 = float(args.ylim2) if (args.ylim2 != "") else np.max(yi) 
    levels = np.linspace(cmaprange1, cmaprange2, args.numcontours) 
    contours = plt.contourf(xi, paramsi, yi, cmap=cmap, levels=levels, extend='both')  
    for contour in contours.collections: contour.set_antialiased(False) ## fix aliasing for old Matplotlib
    cb=plt.colorbar().set_ticks(reasonable_ticks(cmaprange1, cmaprange2, density=.8)) 
    if args.plim1 != "": plt.ylim(ymin=float(args.plim1))
    if args.plim2 != "": plt.ylim(ymax=float(args.plim2))

## ==== Plot tuning and labeling ====
if args.xlim1 != "": plt.xlim(left=float(args.xlim1))
if args.xlim2 != "": plt.xlim(right=float(args.xlim2))
plt.xlabel(xcolname if args.xlabel == '' else args.xlabel) 

if args.contours == 'yes':
    plt.ylabel(args.paramname if args.paramlabel == '' else args.paramlabel) 
    plt.title(args.title if args.title else (ycolname if args.ylabel == '' else args.ylabel)) 
else:
    if args.ylim1 != "": plt.ylim(ymin=float(args.ylim1)) 
    if args.ylim2 != "": plt.ylim(ymax=float(args.ylim2))
    plt.ylabel(args.ylabel if args.ylabel != '' else (ycolname+" (solid), "+ycolname2+" (dashed)") if ycol2 else ycolname) 
    if args.title: plt.title(args.title)
plt.grid()
if args.overlayplot:
    for overlayfunc in args.overlayplot.split(','):
        plt.plot(x, eval(overlayfunc), color='#808080', lw=.5, ls='-.', scaley=False)

try: 
    if not args.contours == 'yes': 
        #plt.legend(prop={'size':12}, loc='upper left') #.draw_frame(False)
        plt.legend(prop={'size':12}, loc='best', fancybox=True, framealpha=0.8)
except:
    pass


## ==== Outputting ====
plt.savefig(args.output, bbox_inches='tight')
# (todo) optional: plt.savefig('%s_%s.png' % (quantity, os.path.split(os.path.dirname(os.getcwd()))[1]), bbox_inches='tight')
