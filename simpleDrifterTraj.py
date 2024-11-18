import numpy as np
import pylab as p
import netCDF4 as nc
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as ftr
import xarray as xr
from loadData import loadGrid
from tqdm import tqdm
from functools import reduce

#load start positions (from /home/break/willlush/workfiles/drifter_validation/drifter_starts_and_traj
print('loading start positions and times (from drifter_starts_and_traj.py)')
saveDir = '/data/break/willlush/drifter_validation/particle_start_positions/'
saveName ='startPositions_10daySeparation_newBathy.npz' 
stLoad = np.load(saveDir+saveName, allow_pickle=True)

stLon = stLoad['lons']
stLat = stLoad['lats']
stID = stLoad['index_in_drifter_array'] #use as proper ids in runs...
stTime = stLoad['times']

#time back to seconds because this is a dumb way to store things...
unixDT = pd.Timestamp('1970-01-01')
stTimes = np.array([(x - unixDT)//pd.Timedelta('1s') for x in stTime])

dDir = '/data/break/willlush/drifter_validation/drifter_traj/'
dName = 'gdp_jul22_ragged_6h.nc'
trData = nc.Dataset(dDir+dName,'r')
idList_compact = trData['ID'][:].data
keepIds = trData['ID'][:][trData['DrogueCenterDepth'][:]==15.].data #15m drogue depth
rowSizeList = trData['rowsize'][:].data
drStat = trData['drogue_status'][:].data.astype(bool)

idList = np.repeat(idList_compact,rowSizeList)
stTraj = idList[stID]
started = np.unique(idList[stID])
assert np.all(np.isin(started,keepIds)) #double-check that start depth is ok

drLon = trData['lon'][:].data
drLat = trData['lat'][:].data
drTime = trData['time'][:].data
#readTime = pd.to_datetime(dr_time,unit='s')

maxAge = 60.*24.*60.*60. #60 days (in seconds)
maxLen = int(maxAge/6/3600)+1 #maximum length of record, assuming 6 hr intervals

lonArr = []
latArr = []
timeArr = []
ageArr = []
indIDArr = []
print('start main loop')
saveDict = {}
counter = 0
for ID in tqdm(started):
    stMsk = stTraj==ID
    lstLon = stLon[stMsk]
    lstLat = stLat[stMsk]
    lstTime = stTimes[stMsk]
    
    idMask = idList==ID
    idLon = drLon[idMask]
    idLat = drLat[idMask]
    idTime = drTime[idMask]
    idDr = drStat[idMask] #drifter attached?

    #get index of starting locations
    chArg1 = np.ravel(np.argwhere(np.isin(idLon,lstLon)))
    chArg2 = np.ravel(np.argwhere(np.isin(idLat,lstLat)))
    chArg3 = np.ravel(np.argwhere(np.isin(idTime,lstTime)))
    common = reduce(np.intersect1d, (chArg1, chArg2, chArg3))

    for ix in common: #make each trajectory
        tLon = idLon[ix:] #lons
        tLat = idLat[ix:] #lats
        tTime = idTime[ix:] #time (s)
        tDr = idDr[ix:] #drogue attached?

        tAge = tTime-tTime[0]
        ageMask = tAge<=maxAge

        cMsk = tDr & ageMask
        tLon = tLon[cMsk]
        tLat = tLat[cMsk]
        tTime = tTime[cMsk]
        tAge = tAge[cMsk]

        saveDict[counter]={'startLoc':(tLon[0],tLat[0]),
                           'lon':tLon,
                           'lat':tLat,
                           'time':tTime,
                           'age':tAge,
                           'traj':ID}
        counter+=1
        

    # #look for time jumps...
    # tDiff = np.diff(idTime)
    
    if False:
        p.figure()
        p.clf()
        cMap = p.axes(projection=ccrs.PlateCarree())
        cMap.add_feature(ftr.COASTLINE,linewidth=0.3,zorder=100)
        cMap.set_global()
        #cMap.plot(stLon,stLat,'r.')
        cMap.plot(idLon[idDr],idLat[idDr])
        cMap.plot(idLon[~idDr],idLat[~idDr])
        cMap.plot(lstLon,lstLat,'r*')
        cMap.plot(idLon[common],idLat[common],'bo')
        #cMap.plot(stLon[ID],stLat[ID],'r*')
        p.show(block=False)

if True:
    np.savez('trajDict_04_04_24.npz',trajDict=saveDict)
