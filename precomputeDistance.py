#Code to get drifter dispersal vectors for both numerical and GDP drifter trajectories and save as complex numbers
import numpy as np
import pylab as p
import netCDF4 as nc
from loadData import loadGrid, loadData
import seawater.extras as swe
from tqdm import tqdm

def getCentroid(lon,lat): #calculating centroid location on sphere
    lon = np.radians(lon)
    lat = np.radians(lat)
    
    x1 = np.cos(lat)*np.cos(lon)
    y1 = np.cos(lat)*np.sin(lon)
    z1 = np.sin(lat)

    xa = np.mean(x1)
    ya = np.mean(y1)
    za = np.mean(z1)

    cLon = np.arctan2(ya,xa)
    hyp = np.sqrt((xa**2) + (ya**2))
    cLat = np.arctan2(za,hyp)
    return(np.degrees(cLon),np.degrees(cLat))

dDir = './' #directory where simplified GDP trajectory data is stored

#load GDP drifter data
trajDict = np.load(dDir+'trajDict_04_04_24.npz',allow_pickle=True) #load simplified GDP drifter trajectory
dDict = trajDict['trajDict'].item() #extract dictionary from .npz file
dId = list(dDict.keys()) #get IDs (keys to dictionary)
dSt = [dDict[x]['startLoc'] for x in dId] #get starting locations
getDictLoc = dict(zip(dSt,dId))  #create dict for quick lookup of dict key given start location (for GDP drifters)

#array of drifter durations (PLD is pelagic larval duration, since we're interested in larvae)
pldArr = np.arange(1,61)

#load numerical drifter starting locations
tMask, tLon, tLat, st2id, id2st = loadGrid(return_starts=True)

#initialize dicts in which to save data
numerical = {}
centroids = {}
drifters = {}
starts = {}
#iterate over all drifter durations
for pld in pldArr:
    idArr, lonArr, latArr = loadData(pld) #load numerical trajectory endpoints

    print('pld is: %s days '%(pld)) #print pld
    #iterate over ids in idArr
    for trID in tqdm(np.unique(idArr)):
        idMask = idArr==trID #only look at matching IDs for num. drifters
        iLon = lonArr[idMask] #all lons for corresponding ID (num. drifters)
        iLat = latArr[idMask] #all lats for corresponding ID (num. drifters)

        cLon,cLat = getCentroid(iLon,iLat) #get centroid of numerical data

        stLon,stLat = id2st[trID] #get ID for GDP drifters using starting location

        mLen,mAng = np.array(list(zip(*[swe.dist([stLat,iLat[x]],[stLon,iLon[x]]) for x in np.arange(len(iLon))]))) #get magnitude and angle of dispersal of numerical drifter (m is for 'modeled')
        xMod = mLen*np.sin(np.radians(mAng)) #get x/real component of numerical dispersal vectors
        yMod = mLen*np.cos(np.radians(mAng)) #get y/imag component of numerical dispersal vectors
        
        wDid = getDictLoc[(stLon,stLat)] #get GDP drifter key for dict

        dAge = dDict[wDid]['age'] #get age of GDP drifters
        if np.amax(dAge)<pld*24*3600: 
            continue #ignore drifters that do not have data at the desired PLD (some GDP drifters didn't have all 60 days of data)
        else:
            dIx = np.ravel(np.argwhere(dAge==pld*24*3600)) #otherwise find data where the age matches (time is in sec, so convert PLD to seconds)

        if len(dIx)==0: #if PLD not present (missing data)
            continue

        dLon = dDict[wDid]['lon'][dIx] #GDP drifter longitude
        dLat = dDict[wDid]['lat'][dIx] #GDP drifter latitude
        if type(dLon)==np.ndarray: #handle array/float issue
            dLon = dLon.item()
            dLat = dLat.item()

        dLen,dAng = swe.dist([stLat,dLat],[stLon,dLon]) #get magnitude and angle of GDP drifter dispersal
        xDri = dLen*np.sin(np.radians(dAng)) #get x/real component of GDP dispersal vector
        yDri = dLen*np.cos(np.radians(dAng)) #get y/imag component of GDP dispersal vector
        
        cLen,cAng = swe.dist([stLat,cLat],[stLon,cLon]) #get magnitude and angle of numerical centroid
        xCen = cLen*np.sin(np.radians(cAng)) #get x/real component of centroid dispersal vector
        yCen = cLen*np.cos(np.radians(cAng)) #get y/imag component of centroid dispersal vector

        if np.any(np.isnan([xDri,yDri,xCen,yCen])): #ignore nans (divide by 0)
            continue

        #get numerical dispersal vectors as array of complex numbers:
        nArr = xMod+1j*yMod
        nArr = nArr.flatten()
        numerical[(pld,trID)] = nArr
        centroids[(pld,trID)] = xCen+1j*yCen #centroid as complex number
        drifters[(pld,trID)] = xDri+1j*yDri #GDP drifter as complex number

#save distance dict as .npz file (may revise in future for a better storage scheme)
if True:        
    np.savez('drifterDistances_6_24_24.npz',numerical=numerical,centroids=centroids,drifters=drifters)
