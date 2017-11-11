# -*- coding: utf-8 -*-
"""
Created on Tue Nov 07 14:00:27 2017

@author: irham

This program is made to combine run_spectral_processing.pro .txt output tables
into one main table. This also generates list of file names to plot by using
gnuplot > load 'plot_list.plt' < command.
"""

# import necessary modules
import time
start = time.time()

import sys
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO
import pandas as pd
import numpy as np
import glob, os

# read list of file names
file_list = np.array([])
plot_list = np.array([])

os.chdir('table/') # change directory

for file in glob.glob('*.txt'):
    file_list = np.append(file_list, file)
    
    # generating list of file names to plot with gnuplot
    plot_list = np.append(plot_list, 
                          'load "plot/'+file.replace('.txt', '.gp"'))

os.chdir('../')

np.savetxt('plot_list.plt', plot_list, fmt='%s')

##file_list = np.loadtxt('result_table/file_name.txt', dtype=str)

# begin the loop
for i in range(len(file_list)):
    print 'Processing table number:', i
    
    # opening each file
    with open('table/'+file_list[i], 'r') as file :
        file_name = file.read()
        file_name = file_name.replace('"', '')
    
    # create new DataFrame or concatenate the exisiting one
    if i == 0:
        data = pd.read_csv(StringIO(file_name), delimiter='|')
        # ordering column names:
        columns = np.append(data.columns.values, 'file_name')
    else:
        temp = pd.read_csv(StringIO(file_name), delimiter='|')
        data = pd.concat((data, temp))
        
# add file names to table
data['file_name'] = file_list

# make a list of columns to keep
with open('result/columns_to_keep.txt', 'r') as file :
    columns_to_keep = file.read().replace('\n', '')\
        .replace(' ', '').replace("'", "").split(',')
        
# saving DataFrame to csv
data.to_csv('result/data_measured.csv', index=False, 
            sep=',', columns=columns_to_keep)

end = time.time()
print 'Finished with elapsed time:', end - start

'''
Further notes:
1. It is important to select only quality flag equals to 0
2. All of luminosities in units of 10**42 erg/s, 
3. while FWHMs and velocity offsets in km/s
'''