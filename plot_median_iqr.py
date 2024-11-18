#code for plotting median normalized dispersal distance (figure 5 in text) and ratio in stochastic dispersal component (figure 6 in text)
import numpy as np
import pylab as p
import netCDF4 as nc
import scipy.stats as sps
from loadData import loadGrid, loadData #uses loadData.py
from sklearn.neighbors import KernelDensity 
import cartopy.crs as ccrs
import cartopy.feature as ftr
import seawater.extras as swe
from tqdm import tqdm

dDir = None #directory where dispersal vectors are stored

#load precomputed dispersal vectors:
dDict = np.load('drifterDistances_6_24_24.npz',allow_pickle=True) 
numerical = dDict['numerical'].item() #load numerical dispersal vectors
drifters = dDict['drifters'].item() #load GDP dispersal vectors
centroids = dDict['centroids'].item() #load centroid dispersal vector

pldArr = np.arange(1,61) #make array of drifter durations (again, PLD stands for pelagic larval duration, since we're thinking about larvae)

#lists to plot median normalized GDP dispersal
median_r = []
median_i = []

#lists to save ratios of iqrs (norm. GDP/norm. numerical)
iqr_ratio = []
iqr_ratio_i = []

#iterate over all drifter durations with progress bar
for pld in tqdm(pldArr):
    idArr, lonArr, latArr = loadData(pld) #load ids for numerical drifters

    dArr = [] #dArr is for GDP drifters (d for drifter)
    cArr = [] #cArr is for numerical centroids (c for centroid)
    mArr = [] #mArr is for numerical drifters (m for model)
    pArr = [] #pArr is for normalized GDP drifters (pArr for plot, because this is the quantity I was originally plotting with this code)
    #iterate over all ids:
    for trID in np.unique(idArr):
        key = (pld,trID) #key in distance dict
        if key in drifters.keys():
            dri = drifters[key][0].flatten() #get GDP drifter disp. vector & flatten
            cen = centroids[key][0].flatten() #get centroid disp. vector
            if cen==0: #ignore non-moving numerical drifters 
                continue
            num = numerical[key][0].flatten() #get numerical disp. vectors
            nMsk = ~np.isnan(num) #ignore nans in numerical disp vectors
            num = num[nMsk]
            if np.any(np.isnan(num/cen)): #catch any nans in normalization
                assert False #(this will throw an assertion error)

            dArr.append(dri.item()) #append GDP drifter disp. vector to dArr
            cArr.append(cen.item()) #append centroid disp vector to cArr
            mArr.extend(num/cen) #normalize numerical drifter dispersal vectors and append to mArr (normalization must happen in each loop for numerical drifters due to array size; probably a better way to handle)

    #make lists into arrays
    dArr = np.array(dArr)
    cArr = np.array(cArr)
    mArr = np.array(mArr)

    #check that the length of GDP and centroid array match
    assert len(cArr)==len(dArr), 'lengths do not match'
    
    pArr = dArr/cArr #normalize GDP drifter disp. vectors by centroid vector
    
    nMask = ~np.isnan(pArr) #mask to remove nans 
    pArr = pArr[nMask]

    iMask = ~np.isinf(pArr) #mask to remove infinite values
    pArr = pArr[iMask]

    realVal = np.real(pArr) #real component of pArr
    imagVal = np.imag(pArr) #imaginary component of pArr
    
    median_r.append(np.median(realVal)) #compute the median value of the real component and append to list for plotting
    median_i.append(np.median(imagVal)) #compute the median value of the imag. component and append to list for plotting
    rv = sps.iqr(realVal) #compute iqr of real component of normalized GDP drifters
    iv = sps.iqr(imagVal) #compute iqr of imag. component of normalized GDP drifters

    mStd_r = sps.iqr(np.real(mArr)) #compute iqr of real component of normalized numerical drifters
    mStd_i = sps.iqr(np.imag(mArr)) #compute iqr of imag. component of normalized numerical drifters
    
    rt = rv/mStd_r #compute ratio of GDP iqr over numerical iqr (real component)
    rti = iv/mStd_i #compute ratio of GDP iqr over numerical iqr (imag. component)
    iqr_ratio.append(rt) #append to list for plotting
    iqr_ratio_i.append(rti)

saveDir = None # directory to save figure

#plotting code for figure 5 (median normalized GDP dispersal) 
if True:
    p.style.use('ggplot')
    p.figure(figsize=(6,4))
    p.clf()
    p.title('Median  ' + r"$\vec {R_{d}^{obs}}}$",fontsize='x-large')
    p.plot(pldArr,median_r,label = 'real')
    p.plot(pldArr,median_i,label = 'imag')
    p.xlabel('drifter duration (days)',fontsize='large')
    p.ylabel(r"$\vec {R_{d}^{obs}}}$       ",rotation='horizontal',fontsize='large')
    p.legend(fontsize='large')
    p.tight_layout()
    p.show()
    p.savefig(saveDir+'medianVals.png')

#plotting code for figure 6 (ratio of iqrs of norm. GDP over norm. numerical dispersal)
if True:
    p.figure(figsize=(6,4))
    p.clf()
    p.title('Ratio of '+r"$L_{diff}^{obs}$"+' to '+r"$L_{diff}^{m}$",fontsize='x-large')
    p.plot(pldArr,std_ratio,label='real')
    p.plot(pldArr,std_ratio_i,label='imaginary')
    p.xlabel('drifter duration (days)',fontsize='large')
    p.ylabel(r"$\frac{L_{diff}^{obs}}{L_{diff}^{m}}$     ",rotation='horizontal',fontsize='xx-large')
    p.legend(fontsize='large')
    p.tight_layout()
    p.show()
    p.savefig(saveDir+'lDiff_ratios.png')
