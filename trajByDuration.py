#Code to add IDs and organize by drifter duration, saves output files in netcdf format
import numpy as np
import pylab as p
import netCDF4 as nc
import xarray as xr
import zarr
import cartopy.crs as ccrs
import cartopy.feature as ftr

yrList = np.arange(2007,2021) #list of years to iterate over

#path to directory for saving output data
saveDict = None

#starting location dicts...
dName = './' #directory for dict of starting positions and corresponding IDs
loadStarts = np.load(dName+'startPosition_idDicts.npz',allow_pickle=True)
st2id = loadStarts['start2id'].item() #dict to go from start location to ID
id2st = loadStarts['id2start'].item() #dict to go from ID to starting location

#iterate over all drifter durations (PLD is for pelagic larval duration, since I'm thinking about larvae)
for PLD in np.arange(1,61,1):
    pldS = PLD*24.*3600. #get drifter duration in seconds

    #initialize lists for saving data
    idArr = []
    lonArr = []
    latArr = []
    
    #loop over years
    for yr in yrList:
        dataDir = None #path to zarr arrays
        if (yr != 2007) and (yr != 2008): #if statement to deal with typo in file naming
            dataName = dataDir+'drifterValidation_run_%s.zarr'%(yr)
        else:
            dataName = dataDir+'drifterValidation_run__%s.zarr'%(yr)
        print(yr) #print year

        #get trajectory endpoints from zarr output files
        dat = zarr.convenience.open(dataName, 'r')
        time = dat.time[:]
        nanMask = ~np.all(np.isnan(time),axis=1)
        time = time[nanMask]
        lat = dat.lat[:][nanMask]
        lon = dat.lon[:][nanMask]
        age = dat.age[:][nanMask]
        maxLen = pldS/(6.*3600.)
        nzLength = np.count_nonzero(~np.isnan(age),axis=1)
        nzMask = nzLength<maxLen

        lat = lat[~nzMask]
        lon = lon[~nzMask]
        age = age[~nzMask]

        ageMask = age==pldS
        doubles = np.count_nonzero(ageMask,axis=1)==2
        lat = lat[~doubles]
        lon = lon[~doubles]
        age = age[~doubles]
        ageMask = ageMask[~doubles]

        ageLon = lon[ageMask]
        ageLat = lat[ageMask]
        included = np.any(ageMask,axis=1)
        #get IDs from each starting position
        lon0 = lon[included][:,0]
        lat0 = lat[included][:,0]
        list2id = list(zip(lon0,lat0))

        #add ids to whichID
        whichID = []
        for pos in list2id:
            if pos in st2id.keys():
                whichID.append(st2id[pos])
            else:
                whichID.append(np.nan)
        #check to make sure lens match
        assert len(ageLon)==len(whichID), 'ID and lon/lat arrs do not correspond'

        #add lats, lons, ids to lists
        idArr.extend(whichID)
        lonArr.extend(ageLon)
        latArr.extend(ageLat)
        
    #convert to arrays, sort by ID,  and save as netcdf
    idArr = np.array(idArr)
    lonArr = np.array(lonArr)
    latArr = np.array(latArr)
    idSort = np.argsort(idArr)

    idArr = idArr[idSort]
    lonArr = lonArr[idSort]
    latArr = latArr[idSort]

    ds = xr.Dataset(data_vars=dict(idArr=(["obs"],idArr),
                                   lonArr=(["obs"], lonArr),
                                   latArr=(["obs"], latArr),),
                    attrs=dict(description="simulated particle positions after %s days"%(PLD)),)
    ds.to_netcdf(saveDict+'particlePositions_pld%02d.nc'%(PLD))


        
