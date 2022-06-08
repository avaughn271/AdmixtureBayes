'''
Author: Peter Wilton
'''


#!/usr/bin/env python

import sys
import glob
import os
import numpy as np
import numpy.random as npr
from functools import partial
from multiprocessing.dummy import Pool
from subprocess import call

cmd = ['treemix']

usage = '''
treemixrunner -p numproc -n numopt [...]

Run treemix numopt times with numproc processors, keeping only the results
files from the run with the best log-likelihood.

    -p numproc             number of processes used [1]
    -n numopt              number of optimizations [1]
    -o                     output prefix for best-likelihood run [Treemix]
    -h                     print this help message
    -seed                  ignored (each run uses a unique seed)
    -print                 print log-likelihoods for each rep
    [...]                  arguments passed to treemix

Requires a version of treemix that has a -seed option.
'''
if len(sys.argv) == 1:
    print(usage)
    sys.exit(0)

i = 1
num_opt = 1
num_processes = 1
out_prefix = 'Treemix'
print_lls = False
while i < len(sys.argv):
    tok = sys.argv[i]
    if tok == '-h':
        print(usage)
        sys.exit(0)
    if tok == '-n':
        i += 1
        num_opt = int(sys.argv[i])
    elif tok == '-print':
        print_lls = True
    elif tok == '-p':
        i += 1
        num_processes = int(sys.argv[i])
    elif tok == '-o':
        i += 1
        out_prefix = sys.argv[i]
    elif tok == '-seed':
        i += 1
    else:
        cmd.append(tok)
    i += 1

tmp_prefix = '_tmp_treemixrunner'

pids = npr.choice(1000000, size = num_opt, replace = False)
prefixes = [tmp_prefix + str(pid) for pid in pids]

cmds = []
for pid, prefix in zip(pids, prefixes):
    cmd_p = cmd[:] + ['-seed', str(pid), '-o', prefix]
    cmds.append(cmd_p)
print(cmds)

pool = Pool(num_processes) # two concurrent commands at a time
fnull = open(os.devnull, 'w')
#removedprin "tmp_prefix, best_prefix, out_prefix", tmp_prefix, best_prefix, out_prefix
try:
    for i, returncode in enumerate(pool.imap(partial(call, shell=False), ['pwd']*10)):
        if returncode != 0:
           print("command %d failed with code: %d" % (i, returncode))
    for i, returncode in enumerate(pool.imap(partial(call, shell=False, stdout = fnull, stderr = fnull), cmds)):
        if returncode != 0:
           print("command %d failed with code: %d" % (i, returncode))
    
    lliks = []
    for repidx, prefix in enumerate(prefixes):
        like_filename = prefix + '.llik'
        with open(like_filename, 'r') as fin:
            for line in fin:
                pass
            llik = float(line.split(':')[1].strip())
        lliks.append(llik)
        if print_lls:
            print('rep', repidx, llik)
    best_prefix_idx = np.array(lliks).argmax()
    best_prefix = prefixes[best_prefix_idx]
    for fn in glob.glob(best_prefix+'*'):
        new_fn = fn.replace(best_prefix, out_prefix)
        print("fn, new_fn",fn, new_fn)
        os.rename(fn, new_fn)

finally:
    for fn in glob.glob(tmp_prefix+'*'):
        os.remove(fn)
    fnull.close()