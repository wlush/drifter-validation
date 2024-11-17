#Module to load some necessary data - loads starts, ID dictionaries, and trajectory starts/endpoints (for a given drifter duration)
import numpy as np
import pylab as p
import netCDF4 as nc
import zarr, pickle

#funtion to load t-points and mask, as well as start location-ID mapping dicts
def loadGrid(return_starts=False):
    #load grid cell center (t-point) locations
    maskName = '/data/guppy2/willlush/Mercator/cGrid/MeshFiles/ext-PSY4V3R1_mask.nc'
    lMask = nc.Dataset(maskName,'r')
    tMask = np.ravel(lMask['tmask'][0,0,:,:].data.astype(bool)) #mask
    tLon = np.ravel(lMask['nav_lon'][:].data) #longitude of t-point
    tLat = np.ravel(lMask['nav_lat'][:].data) #latitude of t-point
    lMask.close() #close netcdf file
    
    if return_starts==True: #load starting location-ID dicts
        #starting location dicts
        dName = './' #directory where start/id dicts are located
        loadStarts = np.load(dName+'startPosition_idDicts.npz',allow_pickle=True)
        st2id = loadStarts['start2id'].item() #dict to go from start to ID
        id2st = loadStarts['id2start'].item() #dict to go from ID to start
        return(tMask,tLon,tLat,st2id,id2st)
    else:
        return(tMask,tLon,tLat)

#function to load trajectory endpoints after a given drifter duration
#pld in code is simply drifter duration, stands for pelagic larval duration
def loadData(pld):
    #load trajectory endpoints
    trajDict = None #where processed trajectories are stored...
    trajName = 'particlePositions_pld%02d.nc'%(pld) #processed trajectory name
    tName = trajDict+trajName
    lTraj = nc.Dataset(tName,'r')
    idArr = lTraj['idArr'][:].data #get IDs from .nc file
    lonArr = lTraj['lonArr'][:].data #get lons from .nc file
    latArr = lTraj['latArr'][:].data #get lats from .nc file
    lTraj.close() #close netcdf file
    return(idArr,lonArr,latArr)
