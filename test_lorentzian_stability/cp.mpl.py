#!/usr/bin/env python
#-*- coding: utf-8 -*-
import matplotlib 
import scipy,numpy
import matplotlib.pyplot as plt

## Initialize 
#matplotlib.rc('text', usetex=True)
#matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman, Times']})
#matplotlib.rc('text.latex',unicode=True)
#colors = ("#BB3300", "#8800DD", "#2200FF", "#0099DD", "#00AA00", "#AA8800",
					   #"#661100", "#440077", "#000088", "#003366", "#004400", "#554400")


plt.figure(figsize=(8,5))

x=numpy.arange(0, 900e9, 10e9)

zdata = []
ylist = []
for y in numpy.arange(-50e-6, 50e-6, 5e-6):
    print y, y/50e-6*90e-9
    znew=numpy.sin(x/(100e9 + y/50e-6*90e9))
    zdata.append(znew)
    ylist.append(y)


CS = plt.contour(x, ylist, zdata, linewidths=0.5,colors='k')
CS = plt.contourf(x, ylist, zdata, numpy.arange(0., 1.03, 0.01), cmap=plt.cm.jet)
#for c in CS.collections: c.set_antialiased(False)  # prevents thin white lines between contours


## Finish the graph 
plt.xlabel(u"") 
plt.ylabel(u"")
#plt.xlim((, ))
plt.legend()
plt.grid()
plt.savefig("cp.png")

