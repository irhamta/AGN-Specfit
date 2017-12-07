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
from scipy.stats import linregress, spearmanr
from astropy.cosmology import FlatLambdaCDM

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
data = pd.read_csv('../result/QSO_Data_v2.csv', low_memory=False)

# compute [S II] total luminosity
data['NA_SII__LUM'] = data['NA_SII_6716__LUM'] + data['NA_SII_6731__LUM']

# convert luminosity unit to log10 unit
lum_col =  [s for s in data.columns if 'LUM' in s and 'ERR' not in s]
lum_err_col =  [s for s in data.columns if 'LUM' in s and 'ERR' in s]

data[lum_col] = np.log10(data[lum_col] * 10.**42.)
data[lum_err_col] = data[lum_err_col] * 10.**42.

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
# Plotting Section  1.
# We will build Diagnostic Diagram here.
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
plt.plot((data['NA_NII_6583__LUM']-data['NA_HA__LUM']).loc[nlq_good & t1_agn],
         (data['NA_OIII_5007__LUM']-data['NA_HB__LUM']).loc[nlq_good & t1_agn], 
         'k.', alpha=0.5)
plt.plot(np.linspace(-3., 0.2, 100), Ke01_a(np.linspace(-3., 0.4, 100)), 'r-')
plt.plot(np.linspace(-3., -0.2, 100), Ka03(np.linspace(-3., 0., 100)), 'b-')
plt.xlabel(r'[N II]/H$\alpha$')
plt.ylabel(r'[O III]/H$\beta$')
plt.xlim(-1.5, 1.2)
plt.ylim(-0.6, 2.)
plt.title(r'[N II] Diagnostic Diagram')
plt.savefig('figures/fig_4.png')


plt.figure(5)
plt.plot((data['NA_SII__LUM']-data['NA_HA__LUM']).loc[nlq_good & t1_agn], 
         (data['NA_OIII_5007__LUM']-data['NA_HB__LUM']).loc[nlq_good & t1_agn], 
         'k.', alpha=0.5)
plt.plot(np.linspace(-4, 0.2, 100), Ke01_b(np.linspace(-4, 0.2, 100)), 'r-')
plt.plot(np.linspace(-0.315, 1.5, 100), Ke06(np.linspace(-0.315, 1.5, 100)), 'g-')
plt.xlabel(r'[S II]/H$\alpha$')
plt.ylabel(r'[O III]/H$\beta$')
plt.title(r'[S II] Diagnostic Diagram')
plt.xlim(-1.6, 1.)
plt.ylim(-0.6, 2.)
plt.savefig('figures/fig_5.png')
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

# define good mask for narrow lines
o2_col = ['NA_OII_3727__QUALITY']
o2_good = (data[o2_col] == 0).all(axis=1) & t1_agn

o3_col = ['NA_OIII_5007__QUALITY']
o3_good = (data[o3_col] == 0).all(axis=1) & t1_agn

o32_col = ['NA_OII_3727__QUALITY', 'NA_OIII_5007__QUALITY']
o32_good = (data[o32_col] == 0).all(axis=1) & t1_agn

hb_col = ['NA_HB__QUALITY']
hb_good = (data[hb_col] == 0).all(axis=1) & t1_agn

ha_col = ['NA_HA__QUALITY']
ha_good = (data[ha_col] == 0).all(axis=1) & t1_agn

n2_col = ['NA_NII_6583__QUALITY']
n2_good = (data[n2_col] == 0).all(axis=1) & t1_agn

s2_col = ['NA_SII_6716__QUALITY', 'NA_SII_6731__QUALITY']
s2_good = (data[s2_col] == 0).all(axis=1) & t1_agn

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



## plot EV diagram
#ev_diagram(**param_o3)
#ev_diagram(**param_o32)
#ev_diagram(**param_o3w)
#ev_diagram(**param_o32w)
#plt.close('all')

'''
# experiment on MgII space
mg2_col = ['BR_MGII_2798__QUALITY', 'IRONUV__QUALITY']
mg2_good = (data[mg2_col] == 0).all(axis=1) # & t1_agn

mg2_good = mg2_good & o3_good

param_mg2w =  {'a'        : (data['IRONUV__EW']
                            / data['BR_MGII_2798__EW']).loc[mg2_good].values,
    
             'b'        : data['BR_MGII_2798__FWHM'].loc[mg2_good].values,
             
             'c'        : np.log10(data['NA_OIII_5007__EW'])\
                             .loc[mg2_good].values,
             
             'mida'     : np.arange(2., 10., 0.2),
             'midb'     : np.arange(0., 10000., 500.),
             'mgrid'    : np.mgrid[2:10:20j, 0:10000:20j],
             'number'   : 21,
             'labelx'   : r'UV $R_{\rm Fe~II}$',
             'labely'   : r'FWHM Mg~II (km/s)',
             'labelz'   : r'log EW$_{\rm [O~III]}$'}

ev_diagram(**param_mg2w)
'''


# calculating black hole mass and Eddington ratio based on Tammour et al. 2015
M_BH = 0.91 + 0.5*np.log10(10**data['CONT5__LUM']/1e+44)\
        + 2*np.log10(data['BR_HB__FWHM']) # in log10 unit
L_bol = 9 * 10**data['CONT5__LUM']
L_Edd = 1.5e+38 * 10**M_BH
Edd_ratio = L_bol/L_Edd


cont_good = (data['CONT5__QUALITY'] == 0) & t1_agn


param_cont11 =  {'a'        : (data['IRONOPT_BR__EW']
                                / data['BR_HB__EW']).loc[cont_good].values,
    
                 'b'        : data['BR_HB__FWHM'].loc[cont_good].values,
                 
                 'c'        : Edd_ratio.loc[cont_good].values,
                 
                 'mida'     : np.arange(-0.1, 11., 0.2),
                 'midb'     : np.arange(0., 15000., 500.),
                 'mgrid'    : np.mgrid[-0.1:11.:50j, 0:15000:30j],
                 'number'   : 11,
                 'labelx'   : r'$R_{\rm Fe~II}$',
                 'labely'   : r'FWHM bH$\beta$ (km/s)',
                 'labelz'   : r'$L/L_{\rm Edd}$'}

param_cont12 =  {'a'        : (data['IRONOPT_BR__EW']
                                / data['BR_HB__EW']).loc[cont_good].values,
    
                 'b'        : data['BR_HB__FWHM'].loc[cont_good].values,
                 
                 'c'        : M_BH.loc[cont_good].values,
                 
                 'mida'     : np.arange(-0.1, 11., 0.2),
                 'midb'     : np.arange(0., 15000., 500.),
                 'mgrid'    : np.mgrid[-0.1:11.:50j, 0:15000:30j],
                 'number'   : 12,
                 'labelx'   : r'$R_{\rm Fe~II}$',
                 'labely'   : r'FWHM bH$\beta$ (km/s)',
                 'labelz'   : r'$\log M_{\rm BH}$'}

param_cont13 =  {'a'        : (data['IRONOPT_BR__EW']
                                / data['BR_HB__EW']).loc[cont_good].values,
    
                 'b'        : data['BR_HB__FWHM'].loc[cont_good].values,
                 
                 'c'        : np.log10(L_bol).loc[cont_good].values,
                 
                 'mida'     : np.arange(-0.1, 11., 0.2),
                 'midb'     : np.arange(0., 15000., 500.),
                 'mgrid'    : np.mgrid[-0.1:11.:50j, 0:15000:30j],
                 'number'   : 13,
                 'labelx'   : r'$R_{\rm Fe~II}$',
                 'labely'   : r'FWHM bH$\beta$ (km/s)',
                 'labelz'   : r'$\log L_{\rm bol}$'}


#ev_diagram(**param_cont11)
#ev_diagram(**param_cont12)
#ev_diagram(**param_cont13)


param_cont21 =  {'a'        : 10.**(data['IRONOPT_BR__LUM']
                                - data['BR_HB__LUM']).loc[cont_good].values,
    
                 'b'        : data['BR_HB__FWHM'].loc[cont_good].values,
                 
                 'c'        : Edd_ratio.loc[cont_good].values,
                 
                 'mida'     : np.arange(-0.1, 8.5, 0.2),
                 'midb'     : np.arange(0., 15000., 500.),
                 'mgrid'    : np.mgrid[-0.1:8.5:43j, 0:15000:30j],
                 'number'   : 21,
                 'labelx'   : r'$L_{\rm Fe~II}/L_{\rm bH\beta}$',
                 'labely'   : r'FWHM bH$\beta$ (km/s)',
                 'labelz'   : r'$L/L_{\rm Edd}$'}

param_cont22 =  {'a'        : 10.**(data['IRONOPT_BR__LUM']
                                - data['BR_HB__LUM']).loc[cont_good].values,
    
                 'b'        : data['BR_HB__FWHM'].loc[cont_good].values,
                 
                 'c'        : M_BH.loc[cont_good].values,
                 
                 'mida'     : np.arange(-0.1, 8.5, 0.2),
                 'midb'     : np.arange(0., 15000., 500.),
                 'mgrid'    : np.mgrid[-0.1:8.5:43j, 0:15000:30j],
                 'number'   : 22,
                 'labelx'   : r'$L_{\rm Fe~II}/L_{\rm bH\beta}$',
                 'labely'   : r'FWHM bH$\beta$ (km/s)',
                 'labelz'   : r'$\log M_{\rm BH}$'}


param_cont23 =  {'a'        : 10.**(data['IRONOPT_BR__LUM']
                                - data['BR_HB__LUM']).loc[cont_good].values,
    
                 'b'        : data['BR_HB__FWHM'].loc[cont_good].values,
                 
                 'c'        : np.log10(L_bol).loc[cont_good].values,
                 
                 'mida'     : np.arange(-0.1, 8.5, 0.2),
                 'midb'     : np.arange(0., 15000., 500.),
                 'mgrid'    : np.mgrid[-0.1:8.5:43j, 0:15000:30j],
                 'number'   : 23,
                 'labelx'   : r'$L_{\rm Fe~II}/L_{\rm bH\beta}$',
                 'labely'   : r'FWHM bH$\beta$ (km/s)',
                 'labelz'   : r'$\log L_{\rm bol}$'}


#ev_diagram(**param_cont21)
#ev_diagram(**param_cont22)
#ev_diagram(**param_cont23)


#==============================================================================
# Plotting Section  3
# We will investigate the Baldwin effect here
#==============================================================================

o = open('../result/output_parameters.csv', 'w')
o.write('number\tm\tc\tr\tp\tstd\trho\tP\n')
o.close()


def linear_plot(x, m, c):
    a = []
    for i in range(len(x)):
        a.append(x[i]*m + c)
    return a

def plot_average(number, x, y, group, labelx, labely):
    
    binned = []    
    median = []
    stdx = []
    stdy = []
    count = []
    
    for i in range(len(group)):
        point = []
        axis = []
        for j in range(len(x)):
            if abs(x[j]-group[i]) <= abs((group[0]-group[1])/2.): 
                axis.append(x[j])
                point.append(y[j])
        
        if len(point) >= 25:
            median.append(np.nanmedian(point))
            stdy.append(np.nanstd(point)*1.)
            stdx.append(np.nanstd(axis)*1.)            
            count.append(len(point))
            binned.append(group[i])
    
    m, c, r, p, std = linregress(x, y)
    rho, P = spearmanr(x, y)
    o = open('../result/output_parameters.csv', 'a')    
    o.write(str(number) + '\t' + str(m) + '\t' + str(c) + '\t' 
            + str(r) + '\t' + str(p) + '\t' + str(std) + '\t' 
            + str(rho) + '\t' + str(P) + '\n')
    o.close()
    plt.figure(number)
    plt.plot(binned, linear_plot(binned, m, c), 'r--', 
             label = ('$\\rho_s = %.2f$\n$ P_s = %.2f$' %(rho, P)))
    plt.legend()
    if number in [1021, 1022]:
        plt.plot(x, y, 'k.', alpha=0.02)
    else:
        plt.plot(x, y, 'k.', alpha=0.1)        
    plt.errorbar(binned, median, yerr=stdy, xerr=0, fmt='bo')
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.xlim(min(binned)-0.05, max(binned)+0.05)
    
    if (2000 < number < 3000) \
    and number not in [2101, 2112, 2201, 2212, 2301, 2312]:
        plt.ylim(41, 48)

    elif (1000 < number < 2000):
        plt.ylim(-2.5, 2)
#        plt.xlim(40, 44)

    else:
        pass
    
    plt.savefig('figures/fig_%i' %number)


LO2 = data['NA_OII_3727__LUM'].loc[o2_good].dropna()
LbHbeta = data['BR_HB__LUM'].loc[LO2.index].values
LO2 = LO2.values
bin_bHbeta = np.arange(min(LbHbeta), max(LbHbeta), 0.25)
plot_average(1021, LbHbeta, LO2-LbHbeta, bin_bHbeta, 
             '$\\log \ L_{\\rm bH\\beta} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_{\\rm [O \ II]}/L_{\\rm bH\\beta}$')

LO3 = data['NA_OIII_5007__LUM'].loc[o3_good].dropna()
LbHbeta = data['BR_HB__LUM'].loc[LO3.index].values
LO3 = LO3.values
bin_bHbeta = np.arange(min(LbHbeta), max(LbHbeta), 0.25)
plot_average(1022, LbHbeta, LO3-LbHbeta, bin_bHbeta, 
             '$\\log \ L_{\\rm bH\\beta} \ \\rm (erg \ s^{-1})$',  
             '$\\log \ L_{\\rm [O \ III]}/L_{\\rm bH\\beta}$')

LnHbeta = data['NA_HB__LUM'].loc[hb_good].dropna()
LbHbeta = data['BR_HB__LUM'].loc[LnHbeta.index].values
LnHbeta = LnHbeta.values
bin_bHbeta = np.arange(min(LbHbeta), max(LbHbeta), 0.25)
plot_average(1023, LbHbeta, LnHbeta-LbHbeta, bin_bHbeta, 
             '$\\log \ L_{\\rm bH\\beta} \ \\rm (erg \ s^{-1})$',  
             '$\\log \ L_{\\rm nH\\beta}/L_{\\rm bH\\beta}$')

LnHalpha = data['NA_HA__LUM'].loc[ha_good].dropna()
LbHbeta = data['BR_HB__LUM'].loc[LnHalpha.index].values
LnHalpha = LnHalpha.values
bin_bHbeta = np.arange(min(LbHbeta), max(LbHbeta), 0.25)
plot_average(1024, LbHbeta, LnHalpha-LbHbeta, bin_bHbeta, 
             '$\\log \ L_{\\rm bH\\beta} \ \\rm (erg \ s^{-1})$',  
             '$\\log \ L_{\\rm nH\\alpha}/L_{\\rm bH\\beta}$')


LNII = data['NA_NII_6583__LUM'].loc[n2_good].dropna()
LbHbeta = data['BR_HB__LUM'].loc[LNII.index].values
LNII = LNII.values
bin_bHbeta = np.arange(min(LbHbeta), max(LbHbeta), 0.25)
plot_average(1025, LbHbeta, LNII-LbHbeta, bin_bHbeta, 
             '$\\log \ L_{\\rm bH\\beta} \ \\rm (erg \ s^{-1})$',  
             '$\\log \ L_{\\rm [N \ II]}/L_{\\rm bH\\beta}$')

LSII = data['NA_SII__LUM'].loc[s2_good].dropna()
LbHbeta = data['BR_HB__LUM'].loc[LSII.index].values
LSII = LSII.values
bin_bHbeta = np.arange(min(LbHbeta), max(LbHbeta), 0.25)
plot_average(1026, LbHbeta, LSII-LbHbeta, bin_bHbeta, 
             '$\\log \ L_{\\rm bH\\beta} \ \\rm (erg \ s^{-1})$',  
             '$\\log \ L_{\\rm [S \ II]}/L_{\\rm bH\\beta}$')

plt.close('all')

#==============================================================================
# Plotting Section  4
# Radio loud and radio quiet AGN
#==============================================================================

# calculate distance in cm
cosmo = FlatLambdaCDM(H0=70., Om0=0.3)
distance = cosmo.luminosity_distance(data['z']).value * 3.08567758*1e+24

# calculate luminosity at 2500 Angstrom,
# determined from AGN continuum luminosity at 5500 Angstrom
def power_law_shift_freq(L2, freq2, freq1, alpha, z):
    L1 = L2 * (freq1/freq2)**alpha / (1.+z)**(1.+alpha)
    return np.log10(L1)

def power_law_shift_wave(L2, wave2, wave1, alpha, z):
    L1 = L2 * (wave1/wave2)**alpha / (1.+z)**(1.+alpha)
    return np.log10(L1)

# the unit is in erg/s
data['CONT_2500A__LUM'] = \
power_law_shift_wave(10.**data['CONT5__LUM'], 
                5100., 2500., data['CONT5__SLOPE'], data['z'])

# continuum flux density at 2500 Angstrom stated in mJy
# 2500 Angstrom = 3e18/2500. Hz
data['CONT_2500A__FLUX'] = 10**data['CONT_2500A__LUM'] / (3e+18/2500) \
                            / (4.*np.pi*distance**2.) * 1e+26

# radio flux still in mJy, need to be converted to erg/s/Hz
data['CONT_20CM__LUM'] = \
np.log10(data['flux_20cm_int'].values * 4.*np.pi*distance**2. * 1e-26)

# calculate flux density at 6 cm
data['CONT_6CM__FLUX'] = 10**power_law_shift_freq(data['flux_20cm_int'], 
                            1.5, 5., -0.5, data['z'])

data['radio_loudness'] = data['CONT_6CM__FLUX']-data['CONT_2500A__FLUX']

end = time.time()
print 'Finished with elapsed time:', end - start

#==============================================================================
# Notes:
#    1. There is something wrong with Radio Loudness calculation
#    2. Immediate solution is by using calculated value from Shen et al. 2011
#    3. CONT5__LUM is good (compared with Shen's data), remember that the unit
#       is in 10^42 erg/s.
#
#
# Useful conversions:
#    1. 20 cm = 1.49896229 GHz
#    2. 6 cm = 4.996540966667 GHz
#    3. c = 3e+10 cm/s = 3e+18 Angstrom/s

# A constant to convert luminosity density from erg/s/Hz to erg/s/Angstrom:
#    > hertz_to_angstrom_const = 3e+18 / (6e+8)**2.
#==============================================================================

# compare with Shen et al. 2011 data
# use calculated RL from Shen et al. 2011
# or define RL as FIRST detected radio sources
data_shen = pd.read_csv('../../../Data/Shen11_data.csv', 
                        usecols = ['RL', 'file_name', 'logL5100', 
                                   'F6cm', 'logFnu'])

data_stern = pd.read_csv('../../../Data/Stern12_data.csv', 
                         usecols = ['file_name', 'logLX'])
    
data = pd.merge(data, data_shen, how='left', on='file_name')
data = pd.merge(data, data_stern, how='left', on='file_name')

# rename X-ray luminosity column
data.rename(columns={'logLX': 'L_x'}, inplace=True)

#plt.close('all')
#plt.loglog(data['radio_loudness'], data['RL'], 'k.')
#plt.loglog(np.arange(0, 1e3, 0.1), np.arange(0, 1e3, 0.1), 'r--')

#plt.plot(np.log10(data['CONT_2500A__FLUX']*1e-26), data['logFnu'], 'k.')
#plt.plot(np.arange(-25, -28, -0.1), np.arange(-25, -28, -0.1), 'r--')

#plt.plot(np.log10(data['flux_20cm_int']), np.log10(data['F6cm']), 'k.')
#plt.plot(np.arange(-0.5, 3.5, 0.1), np.arange(-0.5, 3.5, 0.1), 'r--')

#rl_agn = data['RL'] > 70.
#rq_agn = data['RL'] <= 10.

#==============================================================================

rl_agn = data['radio_loudness'] > 10.
rq_agn = data['radio_loudness'] <= 10.

plt.figure(11270)
plt.plot(10**(data['NA_OIII_5007__LUM'] - data['NA_OII_3727__LUM'])\
         .loc[rq_agn & o32_good], 
         data['CONT_20CM__LUM'].loc[rq_agn & o32_good], 
         'bo', label = 'Radio Quiet', markersize=3)

plt.plot(10**(data['NA_OIII_5007__LUM'] - data['NA_OII_3727__LUM'])\
         .loc[rl_agn & o32_good], 
         data['CONT_20CM__LUM'].loc[rl_agn & o32_good], 
         'r^', label = 'Radio Loud', markersize=5)
plt.xlabel('$L_{\\rm [O \ III]}/L_{\\rm [O \ II]}$', fontsize='x-large')
plt.ylabel('$\\log \ L_{\\rm 1.5 \ GHz} \ \\rm  (erg \ s^{-1} \ \\rm Hz^{-1})$')
plt.legend(loc='best', fontsize='large')
plt.savefig('figures/fig_11270')


plt.figure(11272)
plt.plot(data['NA_OIII_5007__LUM'].loc[rq_agn & o3_good], 
         data['CONT_20CM__LUM'].loc[rq_agn & o3_good], 'bo', 
         label = 'Radio Quiet', markersize=3)

plt.plot(data['NA_OIII_5007__LUM'].loc[rl_agn & o3_good], 
         data['CONT_20CM__LUM'].loc[rl_agn & o3_good], 'r^', 
         label = 'Radio Loud', markersize=5)



###
index_11272_rq = data['NA_OIII_5007__LUM'].loc[rq_agn & o3_good].dropna().index\
    .intersection(data['CONT_20CM__LUM'].loc[rq_agn & o3_good].dropna().index)

print spearmanr(data['NA_OIII_5007__LUM'].loc[index_11272_rq],
                 data['CONT_20CM__LUM'].loc[index_11272_rq])


index_11272_rl = data['NA_OIII_5007__LUM'].loc[rl_agn & o3_good].dropna().index\
    .intersection(data['CONT_20CM__LUM'].loc[rl_agn & o3_good].dropna().index)

print spearmanr(data['NA_OIII_5007__LUM'].loc[index_11272_rl],
                 data['CONT_20CM__LUM'].loc[index_11272_rl])
###


#print linregress((LO3)[fig_11258[0]], Lradio[fig_11258[0]])
plt.xlabel('$\\log \ L_{\\rm [O \ III]} \ \\rm  (erg \ s^{-1})$')
plt.ylabel('$\\log \ L_{\\rm 1.5 \ GHz} \ \\rm  (erg \ s^{-1} \ \\rm Hz^{-1})$')
plt.legend(loc='best', fontsize='large')
plt.savefig('figures/fig_11272')


plt.figure(11274)
plt.plot(data['NA_OII_3727__LUM'].loc[rq_agn & o2_good], 
         data['CONT_20CM__LUM'].loc[rq_agn & o2_good], 'bo', 
         label = 'Radio Quiet', markersize=3)
plt.plot(data['NA_OII_3727__LUM'].loc[rl_agn & o2_good], 
         data['CONT_20CM__LUM'].loc[rl_agn & o2_good], 'r^', 
         label = 'Radio Loud', markersize=5)

###
index_11274_rq = data['NA_OII_3727__LUM'].loc[rq_agn & o2_good].dropna().index\
    .intersection(data['CONT_20CM__LUM'].loc[rq_agn & o2_good].dropna().index)

print spearmanr(data['NA_OII_3727__LUM'].loc[index_11274_rq],
                 data['CONT_20CM__LUM'].loc[index_11274_rq])


index_11274_rl = data['NA_OII_3727__LUM'].loc[rl_agn & o2_good].dropna().index\
    .intersection(data['CONT_20CM__LUM'].loc[rl_agn & o2_good].dropna().index)

print spearmanr(data['NA_OII_3727__LUM'].loc[index_11274_rl],
                 data['CONT_20CM__LUM'].loc[index_11274_rl])
###


plt.xlabel('$\\log \ L_{\\rm [O \ II]} \ \\rm  (erg \ s^{-1})$')
plt.ylabel('$\\log \ L_{\\rm 1.5 \ GHz} \ \\rm  (erg \ s^{-1} \ \\rm Hz^{-1})$')
plt.legend(loc='best', fontsize='large')
plt.savefig('figures/fig_11274')
plt.close('all')


plt.figure(11258)
plt.plot(10.**(data['IRONOPT_BR__LUM'] - 
               data['BR_HB__LUM']).loc[rq_agn & t1_agn], 
         data['BR_HB__FWHM'].loc[rq_agn & t1_agn], 'k.', 
         label = 'Radio Quiet', markersize=3, alpha=0.3)
plt.plot(10.**(data['IRONOPT_BR__LUM'] - 
               data['BR_HB__LUM']).loc[rl_agn & t1_agn], 
         data['BR_HB__FWHM'].loc[rl_agn & t1_agn], 'r^', 
         label = 'Radio Loud', markersize=5, alpha=1.)
plt.xlim(-0.1, 10)
plt.ylim(0, 18000)
plt.xlabel('$\\log \ L_{\\rm [O \ II]} \ \\rm  (erg \ s^{-1})$')
plt.ylabel('$\\log \ L_{\\rm 1.5 \ GHz} \ \\rm  (erg \ s^{-1} \ \\rm Hz^{-1})$')
plt.legend(loc='best', fontsize='large')
plt.savefig('figures/fig_11258')
plt.close('all')


#==============================================================================
# Plotting Section  5
# [O III] / [O II] ratio with AGN Luminosity
#==============================================================================

# lambda * F_lambda = F (integrated)

# SDSS asinh Softening Parameters in maggies
b_u = 1.4 * 1e-10
b_g = 0.9 * 1e-10
b_r = 1.2 * 1e-10
b_i = 1.8 * 1e-10
b_z = 7.4 * 1e-10

# Bandpass efective wavelengths in Angstrom
lambda_fuv = 1516.
lambda_nuv = 2267.

lambda_u = 3543.
lambda_g = 4770.
lambda_r = 6231.
lambda_i = 7625.
lambda_z = 9134.

lambda_j = 1.235e+4
lambda_h = 1.662e+4
lambda_k = 2.159e+4


# convert GALEX magnitude to flux in erg/s/cm^2/Angstrom
flux_fuv = 10**((data['mag_fuv'] - 18.82 - 8.24*data['ebv']) 
            / -2.5 ) * 1.40e-15 * lambda_fuv
flux_nuv = 10**((data['mag_nuv'] - 20.08 - 8.20*data['ebv']) 
            / -2.5 ) * 2.06e-16 * lambda_nuv

# convert 2MASS AB-magnitude to flux in erg/s/cm^2/Angstrom
flux_j = 10**((data['mag_j'] - 0.723*data['ebv'])/-2.5) \
                * 3.129e-13 * 1e-4*1e7 * lambda_j
flux_h = 10**((data['mag_h'] - 0.460*data['ebv'])/-2.5) \
                * 1.133e-13 * 1e-4*1e7 * lambda_h
flux_k = 10**((data['mag_k'] - 0.310*data['ebv'])/-2.5) \
                * 4.283e-14 * 1e-4*1e7 * lambda_k

# convert magnitude to flux
def flux_to_lum(flux, dist):
    return np.log10(flux * 4. * np.pi * dist**2.)

# convert SDSS magnitude to flux in erg/s/cm^2/Angstrom
def sdss_mag_to_lum (mag, b, dist, wave):

    # calculate flux density in erg/s/cm^2/Hz
    flux = (3631 * 1e-23 * 2*b 
            * np.sinh((mag/-2.5) * np.log(10) - np.log(b)))

    # convert flux density to erg/s/cm^2/Angstrom
    flux = flux * 3e+18 / wave**2. * wave

    # calculate luminosity density in unit erg/s/Angstrom
    lum = np.log10(flux * 4. * np.pi * dist**2.)

    return lum

data['L_fuv'] = flux_to_lum(flux_fuv, distance)
data['L_nuv'] = flux_to_lum(flux_nuv, distance)
data['L_u'] = sdss_mag_to_lum(data['mag_u'], b_u, distance, lambda_u)
data['L_g'] = sdss_mag_to_lum(data['mag_g'], b_g, distance, lambda_g)
data['L_r'] = sdss_mag_to_lum(data['mag_r'], b_r, distance, lambda_r)
data['L_i'] = sdss_mag_to_lum(data['mag_i'], b_i, distance, lambda_i)
data['L_z'] = sdss_mag_to_lum(data['mag_z'], b_z, distance, lambda_z)
data['L_j'] = flux_to_lum(flux_j, distance)
data['L_h'] = flux_to_lum(flux_h, distance)
data['L_k'] = flux_to_lum(flux_k, distance)

x_o3    = data[['NA_OIII_5007__LUM', 'L_x']].loc[o3_good].dropna()
fuv_o3  = data[['NA_OIII_5007__LUM', 'L_fuv']].loc[o3_good].dropna()
nuv_o3  = data[['NA_OIII_5007__LUM', 'L_nuv']].loc[o3_good].dropna()
u_o3    = data[['NA_OIII_5007__LUM', 'L_u']].loc[o3_good].dropna()
g_o3    = data[['NA_OIII_5007__LUM', 'L_g']].loc[o3_good].dropna()
r_o3    = data[['NA_OIII_5007__LUM', 'L_r']].loc[o3_good].dropna()
i_o3    = data[['NA_OIII_5007__LUM', 'L_i']].loc[o3_good].dropna()
z_o3    = data[['NA_OIII_5007__LUM', 'L_z']].loc[o3_good].dropna()
j_o3    = data[['NA_OIII_5007__LUM', 'L_j']].loc[o3_good].dropna()
h_o3    = data[['NA_OIII_5007__LUM', 'L_h']].loc[o3_good].dropna()
k_o3    = data[['NA_OIII_5007__LUM', 'L_k']].loc[o3_good].dropna()
rad_o3  = data[['NA_OIII_5007__LUM', 'CONT_20CM__LUM']].loc[o3_good].dropna()

x_o2    = data[['NA_OII_3727__LUM', 'L_x']].loc[o2_good].dropna()
fuv_o2  = data[['NA_OII_3727__LUM', 'L_fuv']].loc[o2_good].dropna()
nuv_o2  = data[['NA_OII_3727__LUM', 'L_nuv']].loc[o2_good].dropna()
u_o2    = data[['NA_OII_3727__LUM', 'L_u']].loc[o2_good].dropna()
g_o2    = data[['NA_OII_3727__LUM', 'L_g']].loc[o2_good].dropna()
r_o2    = data[['NA_OII_3727__LUM', 'L_r']].loc[o2_good].dropna()
i_o2    = data[['NA_OII_3727__LUM', 'L_i']].loc[o2_good].dropna()
z_o2    = data[['NA_OII_3727__LUM', 'L_z']].loc[o2_good].dropna()
j_o2    = data[['NA_OII_3727__LUM', 'L_j']].loc[o2_good].dropna()
h_o2    = data[['NA_OII_3727__LUM', 'L_h']].loc[o2_good].dropna()
k_o2    = data[['NA_OII_3727__LUM', 'L_k']].loc[o2_good].dropna()
rad_o2  = data[['NA_OII_3727__LUM', 'CONT_20CM__LUM']].loc[o2_good].dropna()


x_o32    = data[['NA_OII_3727__LUM', 'NA_OIII_5007__LUM', 'L_x']].loc[o32_good].dropna()
fuv_o32  = data[['NA_OII_3727__LUM', 'NA_OIII_5007__LUM', 'L_fuv']].loc[o32_good].dropna()
nuv_o32  = data[['NA_OII_3727__LUM', 'NA_OIII_5007__LUM', 'L_nuv']].loc[o32_good].dropna()
u_o32    = data[['NA_OII_3727__LUM', 'NA_OIII_5007__LUM', 'L_u']].loc[o32_good].dropna()
g_o32    = data[['NA_OII_3727__LUM', 'NA_OIII_5007__LUM', 'L_g']].loc[o32_good].dropna()
r_o32    = data[['NA_OII_3727__LUM', 'NA_OIII_5007__LUM', 'L_r']].loc[o32_good].dropna()
i_o32    = data[['NA_OII_3727__LUM', 'NA_OIII_5007__LUM', 'L_i']].loc[o32_good].dropna()
z_o32    = data[['NA_OII_3727__LUM', 'NA_OIII_5007__LUM', 'L_z']].loc[o32_good].dropna()
j_o32    = data[['NA_OII_3727__LUM', 'NA_OIII_5007__LUM', 'L_j']].loc[o32_good].dropna()
h_o32    = data[['NA_OII_3727__LUM', 'NA_OIII_5007__LUM', 'L_h']].loc[o32_good].dropna()
k_o32    = data[['NA_OII_3727__LUM', 'NA_OIII_5007__LUM', 'L_k']].loc[o32_good].dropna()
rad_o32  = data[['NA_OII_3727__LUM', 'NA_OIII_5007__LUM', 'CONT_20CM__LUM']].loc[o32_good].dropna()


LO32 = (data['NA_OIII_5007__LUM'] - 
        data['NA_OII_3727__LUM']).loc[o32_good].dropna()

bin_O32 = np.arange(min(LO32), max(LO32), 0.25/2.)
bin_O3 = np.arange(min(LO3), max(LO3), 0.5/2.)
bin_O2 = np.arange(min(LO2), max(LO2), 0.5/2.)

'''
plot_average(2101, x_o3['NA_OIII_5007__LUM'].values, 
             x_o3['L_x'].values, bin_O3, 
             '$\\log \ L_{\\rm [O \ III]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_{\\rm X} \ \\rm (erg \ s^{-1})$')
plot_average(2102, fuv_o3['NA_OIII_5007__LUM'].values, 
             fuv_o3['L_fuv'].values, bin_O3, 
             '$\\log \ L_{\\rm [O \ III]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_{\\rm FUV} \ \\rm (erg \ s^{-1})$')
plot_average(2103, nuv_o3['NA_OIII_5007__LUM'].values, 
             nuv_o3['L_nuv'].values, bin_O3, 
             '$\\log \ L_{\\rm [O \ III]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_{\\rm NUV} \ \\rm (erg \ s^{-1})$')
plot_average(2104, u_o3['NA_OIII_5007__LUM'].values, 
             u_o3['L_u'].values, bin_O3, 
             '$\\log \ L_{\\rm [O \ III]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_u \ \\rm (erg \ s^{-1})$')
plot_average(2105, g_o3['NA_OIII_5007__LUM'].values, 
             g_o3['L_g'].values, bin_O3, 
             '$\\log \ L_{\\rm [O \ III]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_g \ \\rm (erg \ s^{-1})$')
plot_average(2106, r_o3['NA_OIII_5007__LUM'].values, 
             r_o3['L_r'].values, bin_O3, 
             '$\\log \ L_{\\rm [O \ III]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_r \ \\rm (erg \ s^{-1})$')
plot_average(2107, i_o3['NA_OIII_5007__LUM'].values, 
             i_o3['L_i'].values, bin_O3, 
             '$\\log \ L_{\\rm [O \ III]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_i \ \\rm (erg \ s^{-1})$')
plot_average(2108, z_o3['NA_OIII_5007__LUM'].values, 
             z_o3['L_z'].values, bin_O3, 
             '$\\log \ L_{\\rm [O \ III]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_z \ \\rm (erg \ s^{-1})$')
plot_average(2109, j_o3['NA_OIII_5007__LUM'].values, 
             j_o3['L_j'].values, bin_O3, 
             '$\\log \ L_{\\rm [O \ III]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_J \ \\rm (erg \ s^{-1})$')
plot_average(2110, h_o3['NA_OIII_5007__LUM'].values, 
             h_o3['L_h'].values, bin_O3, 
             '$\\log \ L_{\\rm [O \ III]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_H \ \\rm (erg \ s^{-1})$')
plot_average(2111, k_o3['NA_OIII_5007__LUM'].values, 
             k_o3['L_k'].values, bin_O3, 
             '$\\log \ L_{\\rm [O \ III]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_{Ks} \ \\rm (erg \ s^{-1})$')
plot_average(2112, rad_o3['NA_OIII_5007__LUM'].values, 
             rad_o3['CONT_20CM__LUM'].values + np.log10(1.5e+9), bin_O3, 
             '$\\log \ L_{\\rm [O \ III]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_{\\rm radio} \ \\rm (erg \ s^{-1})$')
plt.close('all')


plot_average(2201, x_o2['NA_OII_3727__LUM'].values, 
             x_o2['L_x'].values, bin_O2, 
             '$\\log \ L_{\\rm [O \ II]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_{\\rm X} \ \\rm (erg \ s^{-1})$')
plot_average(2202, fuv_o2['NA_OII_3727__LUM'].values, 
             fuv_o2['L_fuv'].values, bin_O2, 
             '$\\log \ L_{\\rm [O \ II]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_{\\rm FUV} \ \\rm (erg \ s^{-1})$')
plot_average(2203, nuv_o2['NA_OII_3727__LUM'].values, 
             nuv_o2['L_nuv'].values, bin_O2, 
             '$\\log \ L_{\\rm [O \ II]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_{\\rm NUV} \ \\rm (erg \ s^{-1})$')
plot_average(2204, u_o2['NA_OII_3727__LUM'].values, 
             u_o2['L_u'].values, bin_O2, 
             '$\\log \ L_{\\rm [O \ II]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_u \ \\rm (erg \ s^{-1})$')
plot_average(2205, g_o2['NA_OII_3727__LUM'].values, 
             g_o2['L_g'].values, bin_O2, 
             '$\\log \ L_{\\rm [O \ II]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_g \ \\rm (erg \ s^{-1})$')
plot_average(2206, r_o2['NA_OII_3727__LUM'].values, 
             r_o2['L_r'].values, bin_O2, 
             '$\\log \ L_{\\rm [O \ II]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_r \ \\rm (erg \ s^{-1})$')
plot_average(2207, i_o2['NA_OII_3727__LUM'].values, 
             i_o2['L_i'].values, bin_O2, 
             '$\\log \ L_{\\rm [O \ II]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_i \ \\rm (erg \ s^{-1})$')
plot_average(2208, z_o2['NA_OII_3727__LUM'].values, 
             z_o2['L_z'].values, bin_O2, 
             '$\\log \ L_{\\rm [O \ II]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_z \ \\rm (erg \ s^{-1})$')
plot_average(2209, j_o2['NA_OII_3727__LUM'].values, 
             j_o2['L_j'].values, bin_O2, 
             '$\\log \ L_{\\rm [O \ II]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_J \ \\rm (erg \ s^{-1})$')
plot_average(2210, h_o2['NA_OII_3727__LUM'].values, 
             h_o2['L_h'].values, bin_O2, 
             '$\\log \ L_{\\rm [O \ II]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_H \ \\rm (erg \ s^{-1})$')
plot_average(2211, k_o2['NA_OII_3727__LUM'].values, 
             k_o2['L_k'].values, bin_O2, 
             '$\\log \ L_{\\rm [O \ II]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_{Ks} \ \\rm (erg \ s^{-1})$')
plot_average(2212, rad_o2['NA_OII_3727__LUM'].values, 
             rad_o2['CONT_20CM__LUM'].values + np.log10(1.5e+9), bin_O2, 
             '$\\log \ L_{\\rm [O \ II]} \ \\rm (erg \ s^{-1})$', 
             '$\\log \ L_{\\rm radio} \ \\rm (erg \ s^{-1})$')

plt.close('all')


plot_average(2301, x_o32['NA_OIII_5007__LUM'].values-x_o32['NA_OII_3727__LUM'].values, 
             x_o32['L_x'].values, bin_O32, 
             '$\\log \ L_{\\rm [O \ III]}/L_{\\rm [O \ II]}$',
             '$\\log \ L_{\\rm X} \ \\rm (erg \ s^{-1})$')
plot_average(2302, fuv_o32['NA_OIII_5007__LUM'].values-fuv_o32['NA_OII_3727__LUM'].values, 
             fuv_o32['L_fuv'].values, bin_O32, 
             '$\\log \ L_{\\rm [O \ III]}/L_{\\rm [O \ II]}$',
             '$\\log \ L_{\\rm FUV} \ \\rm (erg \ s^{-1})$')
plot_average(2303, nuv_o32['NA_OIII_5007__LUM'].values-nuv_o32['NA_OII_3727__LUM'].values, 
             nuv_o32['L_nuv'].values, bin_O32, 
             '$\\log \ L_{\\rm [O \ III]}/L_{\\rm [O \ II]}$',
             '$\\log \ L_{\\rm NUV} \ \\rm (erg \ s^{-1})$')
plot_average(2304, u_o32['NA_OIII_5007__LUM'].values-u_o32['NA_OII_3727__LUM'].values, 
             u_o32['L_u'].values, bin_O32, 
             '$\\log \ L_{\\rm [O \ III]}/L_{\\rm [O \ II]}$',
             '$\\log \ L_u \ \\rm (erg \ s^{-1})$')
plot_average(2305, g_o32['NA_OIII_5007__LUM'].values-g_o32['NA_OII_3727__LUM'].values, 
             g_o32['L_g'].values, bin_O32, 
             '$\\log \ L_{\\rm [O \ III]}/L_{\\rm [O \ II]}$',
             '$\\log \ L_g \ \\rm (erg \ s^{-1})$')
plot_average(2306, r_o32['NA_OIII_5007__LUM'].values-r_o32['NA_OII_3727__LUM'].values, 
             r_o32['L_r'].values, bin_O32, 
             '$\\log \ L_{\\rm [O \ III]}/L_{\\rm [O \ II]}$',
             '$\\log \ L_r \ \\rm (erg \ s^{-1})$')
plot_average(2307, i_o32['NA_OIII_5007__LUM'].values-i_o32['NA_OII_3727__LUM'].values, 
             i_o32['L_i'].values, bin_O32, 
             '$\\log \ L_{\\rm [O \ III]}/L_{\\rm [O \ II]}$',
             '$\\log \ L_i \ \\rm (erg \ s^{-1})$')
plot_average(2308, z_o32['NA_OIII_5007__LUM'].values-z_o32['NA_OII_3727__LUM'].values, 
             z_o32['L_z'].values, bin_O32, 
             '$\\log \ L_{\\rm [O \ III]}/L_{\\rm [O \ II]}$',
             '$\\log \ L_z \ \\rm (erg \ s^{-1})$')
plot_average(2309, j_o32['NA_OIII_5007__LUM'].values-j_o32['NA_OII_3727__LUM'].values, 
             j_o32['L_j'].values, bin_O32, 
             '$\\log \ L_{\\rm [O \ III]}/L_{\\rm [O \ II]}$',
             '$\\log \ L_J \ \\rm (erg \ s^{-1})$')
plot_average(2310, h_o32['NA_OIII_5007__LUM'].values-h_o32['NA_OII_3727__LUM'].values, 
             h_o32['L_h'].values, bin_O32, 
             '$\\log \ L_{\\rm [O \ III]}/L_{\\rm [O \ II]}$',
             '$\\log \ L_H \ \\rm (erg \ s^{-1})$')
plot_average(2311, k_o32['NA_OIII_5007__LUM'].values-k_o32['NA_OII_3727__LUM'].values, 
             k_o32['L_k'].values, bin_O32, 
             '$\\log \ L_{\\rm [O \ III]}/L_{\\rm [O \ II]}$', 
             '$\\log \ L_{Ks} \ \\rm (erg \ s^{-1})$')
plot_average(2312, rad_o32['NA_OIII_5007__LUM'].values-rad_o32['NA_OII_3727__LUM'].values, 
             rad_o32['CONT_20CM__LUM'].values + np.log10(1.5e+9), bin_O32, 
             '$\\log \ L_{\\rm [O \ III]}/L_{\\rm [O \ II]}$', 
             '$\\log \ L_{\\rm radio} \ \\rm (erg \ s^{-1})$')

plt.close('all')
'''

#==============================================================================
# Saving Final Table
#==============================================================================

#columns_to_keep = [
#        'j', 'object_name', 'plate', 'mjd', 'fiberID', 'z',
#        'L_x', 'L_fuv', 'L_nuv', 'L_u', 'L_g', 'L_r', 'L_i', 'L_z',
#        'L_j', 'L_h', 'L_k', 'CONT_20CM__LUM',
#        'BR_HB__LUM', 'BR_HB__FWHM', 'IRONOPT_BR__LUM', 'CONT5__LUM',
#        'NA_OII_3727__LUM', 'NA_OIII_5007__LUM', 'NA_HB__LUM',
#        'NA_HA__LUM', 'NA_NII_6583__LUM', 'NA_SII__LUM', 'NA_OIII_5007__FWHM'
#        ]


data.loc[t1_agn].to_csv('../result/t1_agn.csv', index=False, 
        sep=',', columns=data.columns)


#==============================================================================

#==============================================================================
# Error propagation calculator
#==============================================================================
#x = 10.**data['NA_OIII_5007__LUM']
#dx = data['NA_OIII_5007__LUM_ERR']
#
#y = 10.**data['BR_HB__LUM']
#dy = data['BR_HB__LUM_ERR']
#
## calculate error propagation for x/y in unit of log10
#def calc_log_err(x, dx):
#    return dx / (x * np.log(10.))
#
#def calc_div_err(x, dx, y, dy):
#    r = x/y # not in log10
#    dr = abs(r) * np.sqrt( (dx/x)**2. + (dy/y)**2. )
#
#    return calc_log_err(np.log10(r), dr)
#
#
#calc_div_err(x.loc[23436], dx.loc[23436], y.loc[23436], dy.loc[23436])
##-0.46999468134153233
##-0.044967431743671793