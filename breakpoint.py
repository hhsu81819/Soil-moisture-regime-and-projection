#!/usr/bin/env python
# coding: utf-8

# In[2]:


from scipy import optimize
import scipy.stats as st
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import sys
import scipy.io
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import os.path
from os import path
import netCDF4
from statsmodels.tsa.stattools import acf
import netCDF4 as nc
import warnings
warnings.filterwarnings("ignore")

# In[3]:


def piecewise_linear(x, x0, y0, k1, k2):
#   Compute optimization of fitting two connected linregs to data.
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

"""
        optimize.curve_fit does the heavy lifting. p & e are 4-element series of results and error/uncertainties for:
        [0] : X or soil moisture value for breakpoint
        [1] : Y (T, flux, etc.) value for breakpoint
        [2] : Slope on left side of breakpoint
        [3] : slope on right side of breakpoint
        piecewise_linear is defined above - the function to optimize over
        Next 2 arguments are X and Y series (sorted on X) of daily data
        p0 is an optional first guess for each of the 4 predicted parameters
        bounds sets limits on the acceptable ranges of each parameter
"""
def piecewise3sg_linear(x, x0, x1, y0, y1, k1, k2, k3):
    condlist = [x < x0, (x >= x0) & (x < x1), x >= x1]
    funclist = [lambda x:k1*(x-x0) + y0, lambda x:k2*(x-x0) + y0, lambda x:k2*(x1-x0) + y0 + k3*(x-x1)]
    return np.piecewise(x, condlist, funclist)
"""
    Hypothesis: Fitting 3-segment regression model to the data assume there are a wilting point and a critical point
    lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0]
    lambda x: k1*x + b, lambda x: k1*x + b + k2*(x-x0), lambda x: k1*x + b + k2*(x-x0) + k3*(x - x1)
    Note that this could generate results that x1<x0. This is not solved analytically but by removing the result of 
    this afterward.
"""



# In[4]:
MODELNAME="AWI-ESM-1-1-LR"
RUN="1pctCO2"
fileSM="/project/cmip5/hhsu/" + MODELNAME + "/mrsos_" + MODELNAME + "_" + RUN + "_r1i1p1f1_regrided2x2_ng.nc"
fileLE="/project/cmip5/hhsu/" + MODELNAME + "/hfls_" + MODELNAME + "_" + RUN + "_r1i1p1f1_regrided2x2_ng.nc"


mat1 = nc.Dataset(fileSM)
mat2 = nc.Dataset(fileLE)
latmax=np.size(mat1['mrsos'][1,:,1])+1
lonmax=np.size(mat1['mrsos'][1,1,:])+1
lenlat=latmax-1
lenlon=lonmax-1
miss_val = -9.99e08
Intercept_1Seg_=np.array([np.full((lenlat,lenlon),miss_val)])
Slope_1Seg_=np.array([np.full((lenlat,lenlon),miss_val)])
DOF_SM_=np.array([np.full((lenlat,lenlon),miss_val)])
RSS_=np.array([np.full((5,lenlat,lenlon),miss_val)])
BIC_=np.array([np.full((5,lenlat,lenlon),miss_val)])


miss_val = -9.99e08
BPx_2Seg_LHSflat_=np.array([np.full((lenlat,lenlon),miss_val)])
BPy_2Seg_LHSflat_=np.array([np.full((lenlat,lenlon),miss_val)])
RHSSlope_2Seg_LHSflat_=np.array([np.full((lenlat,lenlon),miss_val)])


miss_val = -9.99e08
BPx_2Seg_RHSflat_=np.array([np.full((lenlat,lenlon),miss_val)])
BPy_2Seg_RHSflat_=np.array([np.full((lenlat,lenlon),miss_val)])
LHSSlope_2Seg_RHSflat_=np.array([np.full((lenlat,lenlon),miss_val)])


miss_val = -9.99e08
BPx1_3Seg_=np.array([np.full((lenlat,lenlon),miss_val)])
BPy1_3Seg_=np.array([np.full((lenlat,lenlon),miss_val)])
BPx2_3Seg_=np.array([np.full((lenlat,lenlon),miss_val)])
BPy2_3Seg_=np.array([np.full((lenlat,lenlon),miss_val)])
MIDSlope_3Seg_=np.array([np.full((lenlat,lenlon),miss_val)])


for lat in range(14 ,latmax):
    for lon in range(1 ,lonmax):
        if np.array(mat1['mrsos'][1,lat-1,lon-1])<100:
            SM=np.array(mat1['mrsos'][11315:29584,lat-1,lon-1])
            LE=np.array(mat2['hfls'][11315:29584,lat-1,lon-1])
            try:        
                model = LinearRegression()
                x=SM.reshape((-1, 1))
                y=LE.reshape((-1, 1))
                model.fit(x,y)
                Intercept_1Seg_[0,lat-1,lon-1]=model.intercept_[0]
                Slope_1Seg_[0,lat-1,lon-1]=model.coef_[0]
                lacc=acf(SM)[1]
                tau = -1/np.log(np.sqrt(np.sum(lacc*np.abs(lacc))/40))
                dof=np.rint(float(len(SM)) / (tau + 1))
                DOF_SM_[0,lat-1,lon-1]=dof
                RSS_[0,0,lat-1,lon-1]=np.var(np.array(LE))*len(SM)
                Y=SM*float(Slope_1Seg_[0,lat-1,lon-1])+float(Intercept_1Seg_[0,lat-1,lon-1])
                RSS_[0,1,lat-1,lon-1]=np.sum(np.square(LE-Y))
                BIC_[0,0,lat-1,lon-1]=len(SM)*np.log(RSS_[0,0,lat-1,lon-1]/len(SM))+1*np.log(len(SM))
                BIC_[0,1,lat-1,lon-1]=len(SM)*np.log(RSS_[0,1,lat-1,lon-1]/len(SM))+2*np.log(len(SM))
            except:
                pass


# In[7]:

            try:        
                p , e= optimize.curve_fit(piecewise_linear, SM, LE,p0=[(np.max(SM)+np.min(SM))/2,np.median(LE),0,50], bounds=([np.min(SM),np.min(LE),-0.001,1], [np.max(SM),np.max(LE),0.001,1000.0]))
                BPx_2Seg_LHSflat_[0,lat-1,lon-1]=p[0]
                BPy_2Seg_LHSflat_[0,lat-1,lon-1]=p[1]
                RHSSlope_2Seg_LHSflat_[0,lat-1,lon-1]=p[3]
                p0=float(BPx_2Seg_LHSflat_[0,lat-1,lon-1])
                p1=float(BPy_2Seg_LHSflat_[0,lat-1,lon-1])
                p2=0
                p3=float(RHSSlope_2Seg_LHSflat_[0,lat-1,lon-1])
                for i in range(1 ,np.size(LE)+1):
                    if np.array(SM)[i-1]<p0:
                        Y[i-1]=(SM[i-1]-p0)*p2+p1
                    if np.array(SM)[i-1]>p0:
                        Y[i-1]=(SM[i-1]-p0)*p3+p1
                RSS_[0,2,lat-1,lon-1]=np.sum(np.square(LE-Y))
                dof=DOF_SM_[0,lat-1,lon-1]
                BIC_[0,2,lat-1,lon-1]=len(SM)*np.log(RSS_[0,2,lat-1,lon-1]/len(SM))+4*np.log(len(SM))
            except:
                pass


# In[16]:
            try:        
                p , e= optimize.curve_fit(piecewise_linear, SM, LE,p0=[(np.max(SM)+np.min(SM))/2,np.median(LE),1,0], bounds=([np.min(SM),np.min(LE),0,-0.001], [np.max(SM),np.max(LE),1000,0.001]))
                BPx_2Seg_RHSflat_[0,lat-1,lon-1]=p[0]
                BPy_2Seg_RHSflat_[0,lat-1,lon-1]=p[1]
                LHSSlope_2Seg_RHSflat_[0,lat-1,lon-1]=p[2]
                p0=float(BPx_2Seg_RHSflat_[0,lat-1,lon-1])
                p1=float(BPy_2Seg_RHSflat_[0,lat-1,lon-1])
                p2=float(LHSSlope_2Seg_RHSflat_[0,lat-1,lon-1])
                p3=0
                for i in range(1 ,np.size(LE)+1):
                    if np.array(SM)[i-1]<p0:
                        Y[i-1]=(SM[i-1]-p0)*p2+p1
                    if np.array(SM)[i-1]>p0:
                        Y[i-1]=(SM[i-1]-p0)*p3+p1
                RSS_[0,3,lat-1,lon-1]=np.sum(np.square(LE-Y))
                dof=DOF_SM_[0,lat-1,lon-1]
                BIC_[0,3,lat-1,lon-1]=len(SM)*np.log(RSS_[0,3,lat-1,lon-1]/len(SM))+4*np.log(len(SM))
            except:
                pass


# In[19]:
            try:        
                p , e= optimize.curve_fit(piecewise3sg_linear, np.array(SM), np.array(LE),p0=[(np.max(np.array(SM)+0.001)+np.min(np.array(SM)))/2,(np.max(np.array(SM)+0.001)+np.min(np.array(SM)))/2,np.median(np.array(LE)),np.median(np.array(LE)),0,50,0], bounds=([np.min(np.array(SM)),np.min(np.array(SM)),np.min(np.array(LE)),np.min(np.array(LE)),-0.001,1,-0.001], [np.max(np.array(SM)+0.001),np.max(np.array(SM)+0.001),np.max(np.array(LE))+0.001,np.max(np.array(LE))+0.001,0.001,1000.0,0.001]))
                BPx1_3Seg_[0,lat-1,lon-1]=p[0]
                BPy1_3Seg_[0,lat-1,lon-1]=p[2]
                BPx2_3Seg_[0,lat-1,lon-1]=p[1]
                BPy2_3Seg_[0,lat-1,lon-1]=p[2]+ p[5]*(p[1]-p[0])
                MIDSlope_3Seg_[0,lat-1,lon-1]=p[5]
                p0=float(BPx1_3Seg_[0,lat-1,lon-1])
                p1=float(BPy1_3Seg_[0,lat-1,lon-1])
                p2=float(BPx2_3Seg_[0,lat-1,lon-1])
                p3=float(BPy2_3Seg_[0,lat-1,lon-1])   
                p4=0
                p5=float(MIDSlope_3Seg_[0,lat-1,lon-1])
                p6=0
                for i in range(1 ,np.size(LE)+1):
                    if np.array(SM)[i-1]<p0:
                        Y[i-1]=SM[i-1]-p0+p1
                    if np.array(SM)[i-1]>p0 and np.array(SM)[i-1]<p2:
                        Y[i-1]=p3-p5*(p2-SM[i-1])
                    if np.array(SM)[i-1]>p2:
                        Y[i-1]=p3+SM[i-1]-p2
                RSS_[0,4,lat-1,lon-1]=np.sum(np.square(LE-Y))
                dof=DOF_SM_[0,lat-1,lon-1]
                BIC_[0,4,lat-1,lon-1]=len(SM)*np.log(RSS_[0,4,lat-1,lon-1]/len(SM))+6*np.log(len(SM))
            except:
                pass


# In[20]:

ncfilenm='/project/land/hhsu/03.CMIP6/BP_SMxLE_'+ MODELNAME + '_' + RUN + '.nc'
with netCDF4.Dataset(ncfilenm, mode="w", format='NETCDF4') as dds:
    # some file-level meta-data attributes:
    dds.Conventions = "CF-1.6" # if you comply with the convension -- which you 
    dds.title = 'Breakpoint of SM-LE'
    dds.history = ''
    dds.references = ''
    dds.comment = ''

    # defining the dimensions of your arrays:
    time = dds.createDimension('time', Intercept_1Seg_.shape[0])
    lat = dds.createDimension('lat', Intercept_1Seg_.shape[1])
    lon = dds.createDimension('lon', Intercept_1Seg_.shape[2])
    option = dds.createDimension('option', RSS_.shape[1])
    # variables for the columns -- you should use real names
    lat = dds.createVariable('lat',Intercept_1Seg_.dtype, ('lat'))
    lat[:]=0
    lon = dds.createVariable('lon',Intercept_1Seg_.dtype, ('lon'))
    lon[:]=0
    option = dds.createVariable('option',RSS_.dtype, ('option'))
    option[:]=0
    
    Intercept_1Seg = dds.createVariable('Intercept_1Seg',Intercept_1Seg_.dtype, ('time','lat','lon'))
    Intercept_1Seg[:] = Intercept_1Seg_[:]
    ## adds some attibutes
    Intercept_1Seg.units = ''
    Intercept_1Seg.long_name = 'Interception of the 1 segment linear regression'
    Intercept_1Seg.standard_name = 'Intercept_1Seg'  
    
    Slope_1Seg = dds.createVariable('Slope_1Seg',Slope_1Seg_.dtype, ('time','lat','lon'))
    Slope_1Seg[:] = Slope_1Seg_[:]
    ## adds some attibutes
    Slope_1Seg.units = ''
    Slope_1Seg.long_name = 'Slope of the 1 segment linear regression'
    Slope_1Seg.standard_name = 'Slope_1Seg'  

    BPx_2Seg_LHSflat = dds.createVariable('BPx_2Seg_LHSflat',BPx_2Seg_LHSflat_.dtype, ('time','lat','lon'))
    BPx_2Seg_LHSflat[:] = BPx_2Seg_LHSflat_[:]
    ## adds some attibutes
    BPx_2Seg_LHSflat.units = ''
    BPx_2Seg_LHSflat.long_name = 'SM Break Point of the 2 segment LHS-flat linear regression'
    BPx_2Seg_LHSflat.standard_name = 'BPx_2Seg_LHSflat'  
    
    BPy_2Seg_LHSflat = dds.createVariable('BPy_2Seg_LHSflat',BPy_2Seg_LHSflat_.dtype, ('time','lat','lon'))
    BPy_2Seg_LHSflat[:] = BPy_2Seg_LHSflat_[:]
    ## adds some attibutes
    BPy_2Seg_LHSflat.units = ''
    BPy_2Seg_LHSflat.long_name = 'y Break Point of the 2 segment LHS-flat linear regression'
    BPy_2Seg_LHSflat.standard_name = 'BPy_2Seg_LHSflat'
    
    RHSSlope_2Seg_LHSflat = dds.createVariable('RHSSlope_2Seg_LHSflat',RHSSlope_2Seg_LHSflat_.dtype, ('time','lat','lon'))
    RHSSlope_2Seg_LHSflat[:] = RHSSlope_2Seg_LHSflat_[:]
    ## adds some attibutes
    RHSSlope_2Seg_LHSflat.units = ''
    RHSSlope_2Seg_LHSflat.long_name = 'RHS Slope of the 2 segment LHS-flat segment linear regression'
    RHSSlope_2Seg_LHSflat.standard_name = 'RHSSlope_2Seg_LHSflat'  
    
    BPx_2Seg_RHSflat = dds.createVariable('BPx_2Seg_RHSflat',BPx_2Seg_RHSflat_.dtype, ('time','lat','lon'))
    BPx_2Seg_RHSflat[:] = BPx_2Seg_RHSflat_[:]
    ## adds some attibutes
    BPx_2Seg_RHSflat.units = ''
    BPx_2Seg_RHSflat.long_name = 'SM Break Point of the 2 segment RHS-flat linear regression'
    BPx_2Seg_RHSflat.standard_name = 'BPx_2Seg_RHSflat' 
    
    BPy_2Seg_RHSflat = dds.createVariable('BPy_2Seg_RHSflat',BPy_2Seg_RHSflat_.dtype, ('time','lat','lon'))
    BPy_2Seg_RHSflat[:] = BPy_2Seg_RHSflat_[:]
    ## adds some attibutes
    BPy_2Seg_RHSflat.units = ''
    BPy_2Seg_RHSflat.long_name = 'y Break Point of the 2 segment RHS-flat linear regression'
    BPy_2Seg_RHSflat.standard_name = 'BPy_2Seg_RHSflat'  
    
    LHSSlope_2Seg_RHSflat = dds.createVariable('LHSSlope_2Seg_RHSflat',LHSSlope_2Seg_RHSflat_.dtype, ('time','lat','lon'))
    LHSSlope_2Seg_RHSflat[:] = LHSSlope_2Seg_RHSflat_[:]
    ## adds some attibutes
    LHSSlope_2Seg_RHSflat.units = ''
    LHSSlope_2Seg_RHSflat.long_name = 'LHS Slope of the 2 segment RHS-flat linear regression'
    LHSSlope_2Seg_RHSflat.standard_name = 'LHSSlope_2Seg_RHSflat' 
    
    BPx1_3Seg = dds.createVariable('BPx1_3Seg',BPx1_3Seg_.dtype, ('time','lat','lon'))
    BPx1_3Seg[:] = BPx1_3Seg_[:]
    ## adds some attibutes
    BPx1_3Seg.units = ''
    BPx1_3Seg.long_name = 'First SM Break Point of the 3 segment linear regression'
    BPx1_3Seg.standard_name = 'BPx1_3Seg'  
    
    BPy1_3Seg = dds.createVariable('BPy1_3Seg',BPy1_3Seg_.dtype, ('time','lat','lon'))
    BPy1_3Seg[:] = BPy1_3Seg_[:]
    ## adds some attibutes
    BPy1_3Seg.units = ''
    BPy1_3Seg.long_name = 'First LE Break Point of the 3 segment linear regression'
    BPy1_3Seg.standard_name = 'BPy1_3Seg'  

    BPx2_3Seg = dds.createVariable('BPx2_3Seg',BPx2_3Seg_.dtype, ('time','lat','lon'))
    BPx2_3Seg[:] = BPx2_3Seg_[:]
    ## adds some attibutes
    BPx2_3Seg.units = ''
    BPx2_3Seg.long_name = 'Second SM Break Point of the 3 segment linear regression'
    BPx2_3Seg.standard_name = 'BPx2_3Seg'  

    BPy2_3Seg = dds.createVariable('BPy2_3Seg',BPy2_3Seg_.dtype, ('time','lat','lon'))
    BPy2_3Seg[:] = BPy2_3Seg_[:]
    ## adds some attibutes
    BPy2_3Seg.units = ''
    BPy2_3Seg.standard_name = 'BPy2_3Seg'  
    
    MIDSlope_3Seg = dds.createVariable('MIDSlope_3Seg',MIDSlope_3Seg_.dtype, ('time','lat','lon'))
    MIDSlope_3Seg[:] = MIDSlope_3Seg_[:]
    ## adds some attibutes
    MIDSlope_3Seg.units = ''
    MIDSlope_3Seg.standard_name = 'MIDSlope_3Seg'    

    DOF_SM = dds.createVariable('DOF_SM',DOF_SM_.dtype, ('time','lat','lon'))
    DOF_SM[:] = DOF_SM_[:]
    ## adds some attibutes
    DOF_SM.units = ''
    DOF_SM.standard_name = 'DOF_SM'
    
    RSS = dds.createVariable('RSS',RSS_.dtype, ('time','option','lat','lon'))
    RSS[:] = RSS_[:]
    ## adds some attibutes
    RSS.long_name = 'residual sum of square'
    RSS.standard_name = 'RSS'   
    
    BIC = dds.createVariable('BIC',BIC_.dtype, ('time','option','lat','lon'))
    BIC[:] = BIC_[:]
    ## adds some attibutes
    BIC.long_name = 'Bayesian Information Criteria in terms of RSS'
    BIC.standard_name = 'BIC'   
# In[ ]:
