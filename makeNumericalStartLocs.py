#code to find starting locations for numerical drifters - takes GDP drifter trajectories as input and locates points that are on the shelf and separated by other points on the trajectory by at least 10 days. Output is an array of starting lats and lons.
import numpy as np
import pylab as p
import pickle
import zarr
import netCDF4 as nc
import pandas as pd
from tqdm import tqdm
from sklearn.neighbors import BallTree,NearestNeighbors
import cartopy.crs as ccrs
import cartopy

#load GDP data
dDir = None #path to directory where GDP data is saved
dName = 'gdp_jul22_ragged_6h.nc' #GDP data name
trData = nc.Dataset(dDir+dName,'r')
drifterNoList = trData['ID'][:].data #drifter numbers
drifterIDs = np.arange(len(drifterNoList)) #local drifter IDs
keepIds = trData['ID'][:][trData['DrogueCenterDepth'][:]==15.] #select only drifters with 15m drogues
drogueOn = trData['drogue_status'][:].data.astype(bool) #select only drifters with attached drogues
rowSizeList = trData['rowsize'][:].data #rowsize, for data access from ragged array
unqIDs = np.arange(len(drogueOn)) #IDs for drifters with attached drogues

#load depth array, for selecting drifter locations on the shelf
depthDir = '/home/willlush/workfiles/drifter_validation/'
depth = np.load(depthDir+'drifter_depth_mask.npz')['depth']
shallower = depth>=-500. #shallower than 500m (~over shelf)

#load landmask, for preventing stuck drifters
lmDir = '/home/willlush/workfiles/drifter_validation/particle_tracking/'
lm = np.load(lmDir+'depth_avg_landMask.npz')
uM = np.ravel(lm['uMask']) #u mask
vM = np.ravel(lm['vMask']) #v mask

#make map of nearest neighbors for fast lookup of u and v locations (to avoid 'beached' particles)
uLocs = np.radians(list(zip(np.ravel(lm['uLat']),np.ravel(lm['uLon']))))
vLocs = np.radians(list(zip(np.ravel(lm['vLat']),np.ravel(lm['vLon']))))
uN = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(uLocs)
vN = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(vLocs)

#make drifter ID full length of trajectory array
drIDArr = []
for ix in np.arange(len(drifterNoList)):
    sh = rowSizeList[ix]
    drIDArr.extend(np.full(sh,drifterNoList[ix]))
drIDArr = np.array(drIDArr)
#make a mask to keep IDs where the drogue is on and is centered at 15m
drIDmask = np.isin(drIDArr,keepIds)

#extract lat, lon, and time from GDP drifter data
dr_lon = trData['lon'][:].data
dr_lat = trData['lat'][:].data
dr_time = trData['time'][:].data

#make mask for all times after 2007 (excliude old drifters)
readTime = pd.to_datetime(dr_time,unit='s')
timeMask = readTime.year>=2007

#combine masks to ensure drifters included are on the shelf, after 2007, have attached drogues, and drogues are 15m
depthMask = shallower
bigMask = timeMask & depthMask
bigMask = bigMask & drIDmask
bigMask = bigMask & drogueOn

#apply masks to arrays
drIDArr = drIDArr[bigMask]
dr_lon = dr_lon[bigMask]
dr_lat = dr_lat[bigMask]
dr_time = dr_time[bigMask]
readTime = readTime[bigMask]
uId = unqIDs[bigMask]

#make arrays of neighboring points for all lats and lons
dLocs = np.radians(list(zip(dr_lat,dr_lon)))
print('starting neighbor search')
uD, uInd = uN.kneighbors(dLocs)
vD, vInd = vN.kneighbors(dLocs)
print('done')

#ensure velocities exist on nearest grid cell walls - e.g. particle is not started on land or in water shallower than drogue length
uEx = np.any(uM[uInd],axis=1) #is there a u velocity on at least one side?
vEx = np.any(vM[vInd],axis=1) #is there a v velocity on at least one side?

#make mask from above
deep = uEx|vEx

#apply mask
drIDArr = drIDArr[deep]
dr_lon = dr_lon[deep]
dr_lat = dr_lat[deep]
dr_time = dr_time[deep]
readTime = readTime[deep]
uId = uId[deep]

#loop to get start positions
dayDiff = 10. #minimum separation of points on each trajectory
#empty lists to save data
stLon = []
stLat = []
stID = []
stTime = []
tAfter = []
whichDrifters = []
for dr in tqdm(np.unique(drIDArr)):
    td = pd.Timedelta(days=dayDiff) #convert dayDiff into Timedelta object
    idMask = drIDArr==dr #mask to look only at a single trajectory (by ID)
    #apply mask
    whichID = uId[idMask]
    idLon = dr_lon[idMask]
    idLat = dr_lat[idMask]
    idTime = dr_time[idMask]
    rTime = readTime[idMask]
    #find minimum and maximum date/time along trajectory
    minTime = np.amin(rTime)
    maxTime = np.amax(rTime)
    tc = minTime #time value used to loop over trajectory
    stInd = []
    #loop over trajectory to ensure starting locations are at least 10 days apart
    while tc<maxTime:
        #this loop finds the index of the smallest value of time along the trajectory, with a minimum 10-day separation between iterations (enforce by tc) - this ensures that there will always be at least a 10-day separation between starting points, even if the trajectory is discontinuous (e.g. exits and re-enters the shelf)
        tDiff = rTime-tc #time elapsed along trajectory
        gThan = tDiff>=pd.Timedelta(seconds=0) #look only at current and future times (along trajectory)
        minGT = np.amin(rTime[gThan]) #smallest time considered
        localInd = np.argwhere(rTime==minGT).item() #index of smallest time
        stInd.append(localInd) #append starting index to start Index 
        tc = rTime[localInd]+td #advance tc value by 10 days
    if len(stInd)>0: #use starting indices to create a list of starting lats, lons, IDs, and times (time is when the GDP drifter transited a given point)
        stInd = np.unique(stInd)
        stLon.extend(idLon[stInd])
        stLat.extend(idLat[stInd])
        stID.extend(whichID[stInd])
        stTime.extend(rTime[stInd])
        
if True: #save starting positions as .npz file
    print('saving')
    saveDir = None #where to save particle start positions
    saveName ='startPositions_10daySeparation_newBathy.npz' 
    np.savez(saveDir+saveName,lons= stLon,lats = stLat, times= stTime,index_in_drifter_array = stID)
    print('done')

#plot starting locations on the globe
if True:
    p.figure(1,figsize=(10,10))
    p.clf()
    ax = p.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND,zorder=50)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.plot(dr_lon,dr_lat,'.',transform=ccrs.PlateCarree())
    ax.plot(stLon,stLat,'.',transform=ccrs.PlateCarree())
    p.show()
