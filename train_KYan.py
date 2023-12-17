import numpy as np
import dihmm_ext
import bed_writer
import os
from os import listdir,mkdir
from os.path import isfile,join,exists

n_bin_states = 30
n_domain_states = 30
domain_size = 20
tolerance = 1e-6
max_iter = 500


output='/data/sysdir/home/yk890/diHMM/diHMM/chr17_3cells/model'
path='/data/sysdir/home/yk890/diHMM/diHMM/chr17_3cells/data'

tpath=path+'/'
trainfile=[tpath+f for f in listdir(tpath)]
x=dihmm_ext.run_dihmm(n_bin_states, n_domain_states, domain_size, max_iter, tolerance, trainfile)
toutput=output+''
if not os.path.exists(toutput):
    os.mkdir(toutput)
dihmm_ext.save_model(x,toutput)


