#!/usr/bin/env python
#-*- coding: utf-8 -*-

## Import common moduli
import matplotlib, sys, os, time
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, hbar, pi

## Load settings

import matplotlib, sys, os, argparse
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--radius',      type=float, default=33.774e-3, help='cylinder radius in meters')
parser.add_argument('--height',      type=float, default=122.36e-3, help='cylinder height in meters')
args = parser.parse_args()

## Generate and save frequencies
analytic_modes = {}
from scipy.special import jnyn_zeros
max_p, max_n, max_m = 4,4,4
with open('./annotate.txt', 'w') as f:
    for p in range(max_p+1):
        for n in range(max_n):
            S = " "*p           ## Spaces can be used to denote the vertical position of the label
            for m,B in enumerate(jnyn_zeros(n, max_m+1)[0]):
                freq = c/(2*pi) * np.sqrt((B/args.radius)**2 + (p*pi/args.height)**2)
                f.write("%s$TM_{%d%d%d}$%s\t%.6e\n" % (S,n,m+1,p,S, freq))
            for m,B in enumerate(jnyn_zeros(n, max_m+1)[1]):
                if p>0: ## TExx0 modes do not exist [Pozar]
                    freq = c/(2*pi) * np.sqrt((B/args.radius)**2 + (p*pi/args.height)**2)
                    f.write("   %s$TE_{%d%d%d}$%s   \t%.6e\n" % (S,n,m+1,p,S, freq))

