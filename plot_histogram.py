#code for plotting histogram of dispersal vectors for 30-day drifter duration (figure 7)
import numpy as np
import pylab as p
import netCDF4 as nc
import scipy.stats as sps
from loadData import loadGrid, loadData
from sklearn.neighbors import KernelDensity 
import cartopy.crs as ccrs
import cartopy.feature as ftr
import seawater.extras as swe
from tqdm import tqdm

dDir = None #directory where dispersal vectors are stored

#load precomputed location dicts:
dDict = np.load('drifterDistances_6_24_24.npz',allow_pickle=True)
numerical = dDict['numerical'].item()
drifters = dDict['drifters'].item()
centroids = dDict['centroids'].item()

for pld in [30]: #code can be run for other drifter durations (1 to 60 days), currently just doing 30-day duration
    idArr, lonArr, latArr = loadData(pld) #load ids for numerical drifters

    dArr = [] #dArr is for GDP drifters (d for drifter)
    cArr = [] #cArr is for numerical centroids (c for centroid)
    mArr = [] #mArr is for numerical drifters (m for model)
    pArr = [] #pArr is for normalized GDP drifters (pArr for plot, because this is the quantity I was originally plotting with this code)

    #iterate over all ids:
    for trID in np.unique(idArr):
        key = (pld,trID) #key in distance dict
        if (key in drifters.keys()) and (key in centroids.keys()):
            dri = drifters[key][0] #get GDP drifter disp. vector
            cen = centroids[key][0] #get centroid disp. vector
            num = numerical[key][0] #get numerical disp. vectors
            if cen==0: #ignore non-moving numerical drifters 
                continue
            
            dArr.append(dri) #append GDP drifter disp. vector to dArr
            cArr.append(cen) #append centroid disp vector to cArr
            mArr.append(num/cen) #normalize numerical drifter dispersal vectors and append to mArr (normalization must happen in each loop for numerical drifters due to array size; probably a better way to handle)

    #make lists into arrays
    dArr = np.ravel(dArr)
    cArr = np.ravel(cArr)
    mArr = np.ravel(mArr)

    #check that the length of GDP and centroid array match
    assert len(cArr)==len(dArr), 'lengths do not match, check yer shit'

    pArr = dArr/cArr #normalize GDP drifter disp. vectors by centroid vector

    nMask = ~np.isnan(pArr) #mask to remove nans 
    pArr = pArr[nMask]

    iMask = ~np.isinf(pArr) #mask to remove infinite values
    pArr = pArr[iMask]

    rVals = np.real(pArr) #real component of pArr
    iVals = np.imag(pArr) #imaginary component of pArr

    nMask2 = ~np.isnan(mArr) #remove nans from numerical results
    mArr = mArr[nMask2]
    rvm = np.real(mArr) #real component of normalized numerical results
    ivm = np.imag(mArr) #imaginary component of normalized numerical results

    rnM = ~np.isnan(rVals) #removing nans from real and imaginary components
    inM = ~np.isnan(iVals)
    rVals = rVals[np.logical_and(rnM,inM)]
    iVals = iVals[np.logical_and(rnM,inM)]

    rnMv = ~np.isnan(rvm)
    inMv = ~np.isnan(ivm)
    rvm = rvm[np.logical_and(rnMv,inMv)]
    ivm = ivm[np.logical_and(rnMv,inMv)]

    bins = np.linspace(-10.,10.,100) #set up bins for histogram

    #plot histogram of normalized dispersal vectors (real components; figure 7)
    if True:
        p.figure(figsize=(8,7))
        p.clf()
        p.style.use('ggplot')
        p.title(r"$\vec {R_d^{obs}}}$" + "and " + r"$\vec {R_d^{m}}}$" + " histogram for %s-day drifter duration\n(real part only) "%(pld),fontsize='x-large')
        p.hist(rVals,bins=bins,density=True,label='GDP drifters')
        p.hist(rvm,bins=bins,alpha=.5,density=True,label='modeled drifters')
        p.axvline(x=1.,color='grey',linestyle='dashed')
        sX = np.linspace(-10,10,1000)
        p.plot(sX, sps.norm.pdf(sX, 1.0, 1.0),color='k', linestyle='dashed',label='normal distribution\ncentered at x=1')
        p.axis([-10,10,0,1.])
        p.legend(fontsize='large')
        p.xticks(list(p.xticks()[0]) + [1.0])
        p.ylabel('normalized density',fontsize='large')
        p.xlabel('normalized distance',fontsize='large')
        p.tight_layout()
        p.savefig('histogram_30day_realsOnly.png')
        p.show()
