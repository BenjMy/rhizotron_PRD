#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 12:38:56 2022

@author: ben
"""

import glob
import pathlib
import shutil

path = pathlib.Path('../imaging_ERT_MALM/rawData')
files_with_returns = glob.glob(str(path)+'/*8returns.csv')

b_return = [1,4,8,17,20,24,33,36,40,49,52,64]

for fr in files_with_returns:
    for br in b_return:
        my_file = pathlib.Path(fr)
        to_file = pathlib.Path(fr.split('.csv')[0] + '_B'+ str(br) + '.csv')
        
        shutil.copy(my_file, to_file)  # For Python 3.8+.
