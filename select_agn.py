# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 08:41:07 2017

@author: irham
"""
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
import matplotlib as mpl
import scipy.stats as st

start = time.time()

#==============================================================================
# Defining Eigenvector Diagram generator
#==============================================================================

def ev_diagram(a, b, c, mida, midb, number, mgrid, 
             labelx='x', labely='y', labelz='z'):
    xbin = []
    ybin = []
    zbin = []
    stdx = []
    stdy = []
    count = []
    
    for i in range(len(midb)):
        for j in range(len(mida)):
            x_member = [] # x location of member of each grid
            y_member = [] # y location of member of each grid
            z_member = []
                            
            
            for n in range(len(a)):
                
                # compute value if in the grid point
                if abs(a[n]-mida[j]) <= abs((mida[0]-mida[1])/2.)\
                and abs(b[n]-midb[i]) <= abs((midb[0]-midb[1])/2.):
                    x_member.append(a[n])
                    y_member.append(b[n])
                    z_member.append(c[n])
                                        
                                        
            if len(x_member) >= 10 and len(y_member) >= 10 :
                                    
                ybin.append((y_member))
                xbin.append((x_member))
                zbin.append([np.nanmedian(z_member)]*len(y_member))
                stdy.append(np.nanstd(y_member)*1.)
                stdx.append(np.nanstd(x_member)*1.)            
                count.append(len(y_member))
     

    X = [] # 10**(LFeII-LbHbeta)
    Y = [] # fwhm_bHbeta
    Z = [] # LO3

    for i in range(len(xbin)):
        for j in range(len(xbin[i])):
            X.append(xbin[i][j])
            Y.append(ybin[i][j])
            Z.append(zbin[i][j])

    
    fig = plt.figure(number)
    ax = fig.add_subplot(111)

    xx, yy = mgrid # create grid
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([X, Y])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    ax.contour(xx, yy, f, colors='k')
        
    sc = ax.scatter(X, Y, s=1, c=Z, marker='o', cmap=cm.jet, 
                    vmin=np.median(Z) - 3*np.std(Z), 
                    vmax=np.median(Z) + 3*np.std(Z))
    
    col = plt.colorbar(sc)
    col.set_label(labelz)
        
    plt.xlabel(labelx)
    plt.ylabel(labely)

#    plt.legend(loc='best', fontsize='small')

    plt.savefig('figures/fig_%i' %number)
#==============================================================================




#==============================================================================
# For each object with measureable H alpha, we can use diagnostic diagram to
# determine whether it is classified as AGN, Composite, or SF.
# Remember that we have to check the narrow lines fitting quality first.
#==============================================================================

# read data
data = pd.read_csv('result/QSO_Data_v0.csv', low_memory=False)

# compute [S II] total luminosity
data['NA_SII__LUM'] = data['NA_SII_6716__LUM'] + data['NA_SII_6731__LUM']

# convert luminosity unit to log10 unit
lum_col =  [s for s in data.columns if 'LUM' in s and 'ERR' not in s]
data[lum_col] = np.log10(data[lum_col] * 10.**42)

# check narrow lines quality and only select with good quality
# unmeasured lines also have quality flag of 0
nlq_col = ['NA_HB__QUALITY', 'NA_OIII_5007__QUALITY',
           'NA_HA__QUALITY', 'NA_NII_6583__QUALITY',
           'NA_SII_6716__QUALITY', 'NA_SII_6731__QUALITY',
           'NA_OII_3727__QUALITY']

nlq_good = (data[nlq_col] == 0).all(axis=1)

# Type 1 AGN criteria from modified Stern & Laor 2012
t1_agn = (data['BR_HB__QUALITY'] == 0) & (data['BR_HB__FWHM'] > 1000.)\
            & (data['IRONOPT_BR__QUALITY'] == 0)

print 'Number of good NL:', len(data.loc[nlq_good])
print 'Number of good BL:', len(data.loc[t1_agn])

# T1 AGN selection
#data = data.loc[nlq_good & t1_agn]

# define classification lines
def Ka03(x):#x = log(NII/Halpha), y = log(OIII/Hbeta)
    return 0.61/(x-0.05) + 1.3

def Ke01_a(x):#x = log(NII/Halpha), y = log(OIII/Hbeta)
    return 0.61/(x-0.47) + 1.19

def Ke01_b(x):#x = log(SII/Halpha), y = log(OIII/Hbeta)
    return 0.72/(x-0.32) + 1.3

def Ke06(x):#x = log(SII/Halpha), y = log(OIII/Hbeta)
    return 1.89*x + 0.76

#==============================================================================
# Plotting Section  1
# We need to prepare some things like:
#   1. set figure size to 8:6 with 100 dpi resolution
#   2. increase axes font size
#   3. enable text rendering with Latex
#==============================================================================

mpl_param = {'figure.figsize'   : [8.0, 6.0],
             'savefig.dpi'      : 100,
             'axes.titlesize'   : 'xx-large',
             'axes.labelsize'   : 'xx-large',
             'text.usetex'      : True,
             'font.family'      : 'serif'}
mpl.rcParams.update(mpl_param)



# select AGN by using Kewley et al. 2001 criteria
agn_low_z = data['NA_OIII_5007__LUM']-data['NA_HB__LUM'] \
            > Ke01_a(data['NA_NII_6583__LUM']-data['NA_HA__LUM'])


plt.figure(4)
plt.plot((data['NA_NII_6583__LUM']-data['NA_HA__LUM']).loc[nlq_good],
         (data['NA_OIII_5007__LUM']-data['NA_HB__LUM']).loc[nlq_good], 
         'k.', alpha=0.5)
plt.plot(np.linspace(-3., 0.2, 100), Ke01_a(np.linspace(-3., 0.4, 100)), 'r-')
plt.plot(np.linspace(-3., -0.2, 100), Ka03(np.linspace(-3., 0., 100)), 'b-')
plt.xlabel(r'[N II]/H$\alpha$')
plt.ylabel(r'[O III]/H$\beta$')
plt.xlim(-1.5, 1.2)
plt.ylim(-0.6, 2.)
plt.title(r'[N II] Diagnostic Diagram')

plt.figure(5)
plt.plot((data['NA_SII__LUM']-data['NA_HA__LUM']).loc[nlq_good], 
         (data['NA_OIII_5007__LUM']-data['NA_HB__LUM']).loc[nlq_good], 
         'k.', alpha=0.5)
plt.plot(np.linspace(-4, 0.2, 100), Ke01_b(np.linspace(-4, 0.2, 100)), 'r-')
plt.plot(np.linspace(-0.315, 1.5, 100), Ke06(np.linspace(-0.315, 1.5, 100)), 'g-')
plt.xlabel(r'[S II]/H$\alpha$')
plt.ylabel(r'[O III]/H$\beta$')
plt.title(r'[S II] Diagnostic Diagram')
plt.xlim(-1.6, 1.)
plt.ylim(-0.6, 2.)
plt.close('all')

#==============================================================================
# Plotting Section  2
# We will plot some Eigenvector Diagram (EV) here
#
# Further notes:
#   1. It is important to select only quality flag equals to 0
#   2. All of luminosities in units of 10**42 erg/s, 
#   3. FWHMs and velocity offsets in km/s
#==============================================================================


o3_col = ['NA_OIII_5007__QUALITY']
o3_good = (data[o3_col] == 0).all(axis=1) & t1_agn

o32_col = ['NA_OII_3727__QUALITY', 'NA_OIII_5007__QUALITY']
o32_good = (data[o32_col] == 0).all(axis=1) & t1_agn

# set parameter for [O III] EV diagram.
param_o3 =  {'a'        : 10.**(data['IRONOPT_BR__LUM']
                            - data['BR_HB__LUM']).loc[o3_good].values,
    
             'b'        : data['BR_HB__FWHM'].loc[o3_good].values,
             
             'c'        : data['NA_OIII_5007__LUM'].loc[o3_good].values,
             
             'mida'     : np.arange(-0.1, 8.5, 0.2),
             'midb'     : np.arange(0., 15000., 500.),
             'mgrid'    : np.mgrid[-0.1:8.5:43j, 0:15000:30j],
             'number'   : 11278,
             'labelx'   : r'$L_{\rm Fe~II}/L_{\rm bH\beta}$',
             'labely'   : r'FWHM bH$\beta$ (km/s)',
             'labelz'   : r'$\log L_{\rm [O~III]}$ (erg/s)'}

# set parameter for [O III]/[O II] EV diagram
param_o32 = {'a'        : 10.**(data['IRONOPT_BR__LUM']
                            - data['BR_HB__LUM']).loc[o32_good].values,

             'b'        : data['BR_HB__FWHM'].loc[o32_good].values,
             
             'c'        : 10**(data['NA_OIII_5007__LUM']
                            - data['NA_OII_3727__LUM']).loc[o32_good].values,
             
             'mida'     : np.arange(-0.1, 8.5, 0.2),
             'midb'     : np.arange(0., 15000., 500.),
             'mgrid'    : np.mgrid[-0.1:8.5:43j, 0:15000:30j],
             'number'   : 11279,
             'labelx'   : r'$L_{\rm Fe~II}/L_{\rm bH\beta}$',
             'labely'   : r'FWHM bH$\beta$ (km/s)',
             'labelz'   : r'$L_{\rm [O~III]}/L_{\rm [O~II]}$'}


param_o3w =  {'a'        : (data['IRONOPT_BR__EW']
                            / data['BR_HB__EW']).loc[o3_good].values,
    
             'b'        : data['BR_HB__FWHM'].loc[o3_good].values,
             
             'c'        : np.log10(data['NA_OIII_5007__EW'])\
                             .loc[o3_good].values,
             
             'mida'     : np.arange(-0.1, 11., 0.2),
             'midb'     : np.arange(0., 15000., 500.),
             'mgrid'    : np.mgrid[-0.1:11.:50j, 0:15000:30j],
             'number'   : 11280,
             'labelx'   : r'$R_{\rm Fe~II}$',
             'labely'   : r'FWHM bH$\beta$ (km/s)',
             'labelz'   : r'log EW$_{\rm [O~III]}$'}


param_o32w = {'a'        : (data['IRONOPT_BR__EW']
                            / data['BR_HB__EW']).loc[o32_good].values,
    
             'b'        : data['BR_HB__FWHM'].loc[o32_good].values,
             
             'c'        : np.log10(data['NA_OIII_5007__EW']
                                     / data['NA_OII_3727__EW'])\
                                     .loc[o32_good].values,
             
             'mida'     : np.arange(-0.1, 11., 0.2),
             'midb'     : np.arange(0., 15000., 500.),
             'mgrid'    : np.mgrid[-0.1:11.:50j, 0:15000:30j],
             'number'   : 11281,
             'labelx'   : r'$R_{\rm Fe~II}$',
             'labely'   : r'FWHM bH$\beta$ (km/s)',
             'labelz'   : r'log(EW$_{\rm [O~III]}$/EW$_{\rm [O~II]}$)'}



# plot EV diagram
#ev_diagram(**param_o3)
#ev_diagram(**param_o32)
#ev_diagram(**param_o3w)
#ev_diagram(**param_o32w)
plt.close('all')

end = time.time()
print 'Finished with elapsed time:', end - start

#==============================================================================