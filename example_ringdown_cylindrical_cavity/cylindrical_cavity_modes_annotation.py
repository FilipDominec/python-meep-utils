#!/usr/bin/env python
#-*- coding: utf-8 -*-

## Import common moduli
import matplotlib, sys, os, time
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, hbar, pi

## User settings
radius=33.774e-3
height=122.36e-3

analytic_modes = {}
from scipy.special import jnyn_zeros
max_p, max_n, max_m = 4,4,4
with open('./annotate.txt', 'w') as f:
    for p in range(max_p+1):
        for n in range(max_n):
            S = " "*p           ## Spaces can be used to denote the vertical position of the label
            for m,B in enumerate(jnyn_zeros(n, max_m+1)[0]):
                freq = c/(2*pi) * np.sqrt((B/radius)**2 + (p*pi/height)**2)
                f.write("%s$TM_{%d%d%d}$%s\t%.6e\n" % (S,n,m+1,p,S, freq))
            for m,B in enumerate(jnyn_zeros(n, max_m+1)[1]):
                if p>0: ## TExx0 modes do not exist [Pozar]
                    freq = c/(2*pi) * np.sqrt((B/radius)**2 + (p*pi/height)**2)
                    f.write("   %s$TE_{%d%d%d}$%s   \t%.6e\n" % (S,n,m+1,p,S, freq))

