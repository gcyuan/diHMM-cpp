#!/usr/bin/env python
#--coding:utf-8 --


__author__ = "Xuan Cao"
__date__ = "2023-12-19"
__email__ = "caoxuan"

import argparse
from datetime import datetime
import numpy as np
import dihmm_ext
import bed_writer
import os
from os import listdir,mkdir
from os.path import isfile,join,exists


def _parse_args():
    """
    Help informations!
    """
    description = "Train diHMM runner."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", dest="input_dir", required=True,type=str, help="The input binarized files dir. File name: X1_chr1_binary.txt.")
    parser.add_argument("--clusters", dest="clusters", required=True, type=str, help="Clusters/cell_types names used to train model. Example: X1,X2 .")
    parser.add_argument("--chroms", dest="chroms", required=True, type=str, help="chrs used to train model. Example: chr1,chr2 .")
    parser.add_argument("-o", dest="out_dir", required=True, type=str, help=" Output dir." )
    parser.add_argument("--n_bin_states", dest="n_bin_states", default=2, type=int, help="Number of bin states. Default=2." )
    parser.add_argument("--n_domain_states", dest="n_domain_states", default=4, type=int, help="Number of domain states. Default=4." )
    parser.add_argument("--domain_size", dest="domain_size", default=8, type=int, help="Number of domain size. Default=8." )
    parser.add_argument("--tolerance", dest="tolerance", default=1000000,  type=int, help="Number of bin states. Default=1e-6." )
    parser.add_argument("--max_iter", dest="max_iter", default=500, type=int, help="Max iter number. Default=500." )
    parser.add_argument("--bin_res", dest="bin_res", default=500, type=int, help="bin length used to generate binarized files. Default=500." )
    opt = parser.parse_args()
    return opt


def run_dihmm_main(input_dir, cell_types, chrs, output_dir, n_bin_states, n_domain_states, domain_size, tolerance, max_iter, bin_res):
    output_dir = os.getcwd() + '/' + output_dir
    input_dir = os.getcwd() + '/' + input_dir
 
    anno_dir = output_dir + '/anno/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(anno_dir):
        os.makedirs(anno_dir)
    
    os.chdir(output_dir)
    trainfile = [input_dir+'/' + ct + '_' + chrom + '_binary.txt' for ct in cell_types for chrom in chrs]
    x=dihmm_ext.run_dihmm(n_bin_states, n_domain_states, domain_size, max_iter, tolerance, trainfile)
    
    
    print('Start annotation')
    #### annotation
    for ct in cell_types:
        for chrom in chrs:
            binarized_file = input_dir + '/' + ct + '_' + chrom + '_binary.txt' 
            a = dihmm_ext.annotate(x,[binarized_file])
            b = bed_writer.BedWriter(a[0], x)
            print('start write bed')
            b.write_bed_files(anno_dir,ct,chrom, bin_res)
            print( 'process'+ anno_dir+'/'+ct+ '_'+ chrom)



def main():
    opt = _parse_args()
    cell_types = opt.clusters.split(',')
    chroms = opt.chroms.split(',')
    run_dihmm_main(opt.input_dir, cell_types, chroms, opt.out_dir, opt.n_bin_states, opt.n_domain_states, opt.domain_size, opt.tolerance, opt.max_iter, opt.bin_res)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    time_caused = datetime.now() - start_time
    print("The time caused is : ", time_caused, '\n')




    





