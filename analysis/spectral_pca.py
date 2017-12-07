# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:47:53 2017

@author: irham
"""

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
from sklearn import preprocessing

columns_to_keep = [
'IRONUV__LUM_ERR',
'IRONUV__EW', 'IRONUV__EW_ERR', 'IRONUV__QUALITY',


'IRONOPT_BR__LUM', 'IRONOPT_BR__LUM_ERR',
'IRONOPT_BR__EW', 'IRONOPT_BR__EW_ERR', 'IRONOPT_BR__QUALITY',

'BR_MGII_2798__LUM', 'BR_MGII_2798__LUM_ERR',
'BR_MGII_2798__FWHM', 'BR_MGII_2798__FWHM_ERR',
'BR_MGII_2798__VOFF', 'BR_MGII_2798__VOFF_ERR',
'BR_MGII_2798__EW', 'BR_MGII_2798__EW_ERR',
'BR_MGII_2798__QUALITY',


'NA_OII_3727__LUM', 'NA_OII_3727__LUM_ERR',
'NA_OII_3727__FWHM', 'NA_OII_3727__FWHM_ERR',
'NA_OII_3727__VOFF', 'NA_OII_3727__VOFF_ERR',
'NA_OII_3727__EW', 'NA_OII_3727__EW_ERR',
'NA_OII_3727__QUALITY', 

'BR_HB__LUM', 'BR_HB__LUM_ERR',
'BR_HB__FWHM', 'BR_HB__FWHM_ERR',
'BR_HB__VOFF', 'BR_HB__VOFF_ERR',
'BR_HB__EW', 'BR_HB__EW_ERR',
'BR_HB__QUALITY',

'NA_HB__LUM', 'NA_HB__LUM_ERR',
'NA_HB__FWHM', 'NA_HB__FWHM_ERR',
'NA_HB__VOFF', 'NA_HB__VOFF_ERR',
'NA_HB__EW', 'NA_HB__EW_ERR',
'NA_HB__QUALITY',

'NA_OIII_5007__LUM', 'NA_OIII_5007__LUM_ERR',
'NA_OIII_5007__FWHM', 'NA_OIII_5007__FWHM_ERR',
'NA_OIII_5007__VOFF', 'NA_OIII_5007__VOFF_ERR',
'NA_OIII_5007__EW', 'NA_OIII_5007__EW_ERR',
'NA_OIII_5007__QUALITY',
        ]

data = pd.read_csv('../result/t1_agn.csv', usecols=columns_to_keep)

qcol = [s for s in data.columns if 'QUALITY' in s]
qcol_good = (data[qcol] == 0).all(axis=1)

ecol = [s for s in data.columns if 'ERR' in s]

#data['LOIII'] = data['NA_OIII_5007__LUM']
#data['LOIII_LOII'] = data['NA_OIII_5007__LUM'] - data['NA_OII_3727__LUM']
#data['LOIII_nHbeta'] = data['NA_OIII_5007__LUM'] - data['NA_HB__LUM']
#data['LFeII_LbHbeta'] = data['IRONOPT_BR__LUM'] - data['BR_HB__LUM']
#data['LNII_nHalpha'] = data['NA_NII_6583__LUM'] - data['NA_HA__LUM']
#data['LSII_nHalpha'] = data['NA_SII__LUM'] - data['NA_HA__LUM']
#data['FWHM_bHbeta'] = data['BR_HB__FWHM']
#data['FWHM_OIII'] = data['NA_OIII_5007__FWHM']

data.drop(np.concatenate((qcol, ecol)), axis=1, inplace=True)

#data.replace([-np.inf, np.inf], np.nan, inplace=True)



data.dropna(inplace=True)

print data.corr(method='spearman')

data_scaled = preprocessing.StandardScaler().fit_transform(data)


pca = PCA(n_components=5)
X = pca.fit_transform(data_scaled)

print('\n\nexplained variance ratio (first two components): \n%s'
      % str(pca.explained_variance_ratio_))


print '\n========================================\n\n'

# eigenvalue
print '\nVariance (Eigenvalue)'
print pca.explained_variance_

# Percentage of variance explained for each components
print '\nVariance Ratio (Proportion)'
print pca.explained_variance_ratio_

#print '\nMean'
#print pca2.mean_

# proportion of 
print '\nBR_HB__LUM', 'BR_HB__FWHM', 'IRONOPT_BR__LUM'

print '\nComponents'
print pca.components_

print '\nCumulative Variance (Cumulative)'
print pca.explained_variance_ratio_.cumsum()

import matplotlib.pyplot as plt

plt.figure()
plt.scatter(X[:, 0], X[:, 1])
plt.xlabel('Component 1')
plt.ylabel('Component 2')



o = open('../result/pca_result.csv', 'w')

o.write('---\tEV1\tEV2\tEV3\tEV4\tEV5\n')
o.write('Eigenvalue\t')
for i in pca.explained_variance_:
    o.write(str(np.round(i, 3)) + '\t')

o.write('\n')

o.write('Proportion\t')
for i in pca.explained_variance_ratio_:
    o.write(str(np.round(i, 3)) + '\t')

o.write('\n')

o.write('Cumulative\t')
for i in pca.explained_variance_ratio_.cumsum():
    o.write(str(np.round(i, 3)) + '\t')

o.write('\n')



for i in range(len(data.columns)):
    o.write(str(data.columns[i]) + '\t')
    for j in range(pca.n_components):
        o.write(str(np.round(pca.components_[j][i], 3)) + '\t')
    o.write('\n')

o.close()