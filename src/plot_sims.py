#!/usr/bin/env python

from __future__ import print_function
import pandas as pd
import sys
import re, math, os
import numpy as np
import matplotlib.pyplot as plt


from matplotlib import rcParams
rcParams['font.family'] = 'Myriad Pro'

filename = sys.argv[1]

print(filename)
dat = pd.read_csv(filename, nrows=30, sep=";")

colnames = dat.columns.values

# sort out cl
pf_names = [ x for x in colnames if re.match("p_\d+$",x) != None ]

pf_names_i = [ int(re.sub("p_(\d+)","\\1",x)) for x in pf_names ]

sub = dat[pf_names]


def abssqrt(x):
    return(x*x)

# x-root transform the data
sub2 = sub.applymap(lambda x: x**.1)
sub2 = np.transpose(np.array(sub2))

num_rows = 4 

plt.figure(figsize=(8,16))
plt.subplot(num_rows,1,1)
plt.plot(dat["generation"], dat["mean_p"],"r")
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'allele frequency $p_{i}$')
plt.xlabel(r'generation, $t$')
plt.ylim((-0.1,1.1))
plt.xlim((min(dat["generation"]),max(dat["generation"])))

plt.subplot(num_rows,1,2)
plt.ylim((max(pf_names_i),min(pf_names_i)))
imgrplot = plt.imshow(sub2,cmap=plt.cm.jet, aspect='auto',interpolation='none',extent=[min(dat["generation"]),max(dat["generation"]),min(pf_names_i),max(pf_names_i)])
plt.ylabel(r'copy number $c_{i}$')
plt.xlabel(r'generation, $t$')

# plot first row with frequencies
rowdat = dat.loc[0,pf_names]


plt.subplot(num_rows,1,3)
plt.plot(range(0,len(rowdat)), rowdat,"b")
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'freq. cells')
plt.xlabel(r'number parent-1 mitochondria per cell, $i$')
plt.text(5,0.35,"within-spore distribution after 30 generations")
plt.ylim((-0.1,0.5))
plt.xlim((0,len(rowdat)))

rowdat = dat.loc[29,pf_names]
plt.subplot(num_rows,1,4)
plt.plot(range(0,len(rowdat)), rowdat,"m")
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'freq. cells')
plt.xlabel(r'number parent-1 mitochondria per cell, $i$')
plt.text(5,0.35,"within-cell distribution after 30 generations")
plt.ylim((-0.1,0.5))
plt.xlim((0,len(rowdat)))

graphname = "graph_" + os.path.basename(filename) + ".pdf"
plt.subplots_adjust(hspace=.3)
plt.savefig(graphname,format="pdf")

