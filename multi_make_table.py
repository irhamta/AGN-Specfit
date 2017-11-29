# -*- coding: utf-8 -*-
"""
Created on Tue Nov 07 14:00:27 2017

@author: Irham
This the parallelized version of make_table.py to speed up the computation.
This program is made to combine process_spectra.pro .txt output tables
into one main table. This also generates list of file names to plot by using
gnuplot > load 'plot_list.plt' < command.
"""

# import necessary modules
import multiprocessing
import time
start = time.time() # to calculate elapsed time

import sys
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO
import pandas as pd
import numpy as np
import glob, os

import tqdm


# create process function to concat pandas tables
def process(name):        
    # opening each file
    with open('table/'+name, 'r') as file :
        file_name = file.read()
        file_name = file_name.replace('"', '')

    # load DataFrame and pass it to result
    temp = pd.read_csv(StringIO(file_name), delimiter='|')
    return temp


if __name__ == '__main__':

    
#==============================================================================
#     Part 1. Make list of file names and generate plots
#==============================================================================
    
    print 'Make list of file names and generate plots.....'
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

    # limit file_list for experiment
    # file_list = file_list[:1000]
    
    
#==============================================================================
#   Part 2. Prepare and execute multiprocess pool to concatenate pandas tables
#==============================================================================
    
    print 'Execute multiprocess pool to concatenate pandas tables.....'
    n_cpu = 3
    pool = multiprocessing.Pool(n_cpu)

    temp_data = []
    for _ in tqdm.tqdm(pool.imap(process, file_list), total=len(file_list)):
        temp_data.append(_)
        pass
    print 'Concatenating.....'
    data = pd.concat(temp_data)
    
    # close the pool
    pool.close()
    pool.join()
    
    print data

    
#==============================================================================
#   Part 3. Insert file name, clean data, and save to .csv file
#==============================================================================

    print 'Insert file name, clean data, and save to .csv file'
    # remove .txt extension in file_name
    for i in range(len(file_list)):
        file_list[i] = file_list[i][:-4]
    
    # make a list of columns to keep
    with open('result/columns_to_keep.txt', 'r') as file :
        columns_to_keep = file.read().replace('\n', '')\
            .replace(' ', '').replace("'", "").split(',')
    
    # add file names to table
    data['file_name'] = file_list
            
    # saving DataFrame to csv
    data.to_csv('result/data_result.csv', index=False, 
                sep=',', columns=columns_to_keep)

        
#==============================================================================
#   Delete following lines if you are not Irham
#==============================================================================
    
    print 'Combining tables.....'
    
    # left join master table with new result table
    data_1 = pd.read_excel('../../Data/QSO_Sample.xlsx')
    data_2 = pd.read_csv('result/data_result.csv', 
                         skipinitialspace=True) # to remove whitespace
    data_master = pd.merge(data_1, data_2, how='left', on='file_name')
    
    data_master.to_csv('result/QSO_Data_v0.csv', index=False, sep=',')
    
    end = time.time()
    print 'Finished with elapsed time:', end - start
    
#==============================================================================
#   Further notes:
#       1. It is important to select only quality flag equals to 0
#       2. All of luminosities in units of 10**42 erg/s, 
#       3. while FWHMs and velocity offsets in km/s
#==============================================================================



#==============================================================================
#   SDSS name generator
#==============================================================================
#    from astropy import units as u
#    from astropy.coordinates import SkyCoord
#    
#    def name(x, y):
#        c = SkyCoord(ra=x*u.degree, dec=y*u.degree, frame='icrs')
#        if y >= 0:
#            return str('SDSS J%02d%02d%05.2f+%02d%02d%04.1f' %(c.ra.hms[0], 
#            c.ra.hms[1], c.ra.hms[2], c.dec.dms[0], c.dec.dms[1], c.dec.dms[2]))
#        else:
#            return str('SDSS J%02d%02d%05.2f-%02d%02d%04.1f' %(c.ra.hms[0], 
#            c.ra.hms[1], c.ra.hms[2], abs(c.dec.dms[0]), abs(c.dec.dms[1]), 
#            abs(c.dec.dms[2])))
#    
#    object_name = []
#    for j in range(len(data)):
#        print j
#        object_name.append(name(data['ra'].loc[j], data['dec'].loc[j]))
#    
#    data['object_name'] = object_name
#    data['object_name'].to_csv('object_name.csv', index=False, sep=',')
#==============================================================================