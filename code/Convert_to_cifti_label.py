#!/usr/bin/env python

"""
Created on Tue March 23 12:30:54 2019

@author: Guoyuan Yang
"""

# script for calculating average of gifti sulcal depth maps

import numpy as np
import sys
from sys import argv
# sys.path.append('/disk2/guoyuan/software/conde/install_path/anaconda2/lib/python2.7/site-packages/nibabel'
import nibabel as nib
import nibabel.gifti as nibgif
import os
import copy
import csv
from scipy.io import loadmat


def usage():
    print ("Usage: " + argv[0] + " <folder where data metrics are> <weightsFile> <hemisphere> <output filename>")
    sys.exit(1)

if len(argv) < 2:
    usage()

mat_file=argv[1]
orig_label_name_l=argv[2]
orig_label_name_r=argv[3]
save_label_name_l=argv[4]
save_label_name_r=argv[5]

read_mat_file=loadmat(mat_file)
lh_label_new=read_mat_file["lh_labels"]
rh_label_new=read_mat_file["rh_labels"]

# lh_label_old=nibgif.giftiio.read(orig_label_name_l)
# rh_label_old=nibgif.giftiio.read(orig_label_name_r)
lh_label_old=nib.load(orig_label_name_l)
rh_label_old=nib.load(orig_label_name_r)

lh_label_old.darrays[0].data=np.int32(lh_label_new)
rh_label_old.darrays[0].data=np.int32(rh_label_new)

nib.save(lh_label_old,save_label_name_l)
nib.save(rh_label_old,save_label_name_r)
