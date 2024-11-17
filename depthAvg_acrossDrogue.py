#Code for averaging u and v mercator velocities across 6m drogue length, centered at 15m
#generates daily depth averaged netcdf files, then combines into a single file for each year.
import numpy as np
import pylab as p
import xarray as xr
import netCDF4 as nc
import datetime as dt
import matplotlib.dates as mdates
from dask.diagnostics import ProgressBar
import pandas as pd
import time
import sys
import os

from scipy import interpolate

if True:
    print('code will not run without a path to Mercator velocity files')

#load velocity data
#dataDir='/data/ripBig/pringle/mercator/'
dataDir=None #path to velocity data from Mercator physical model

#where to save data + data name
saveDir = None #where to save data
whatYear = 2010 #year over which to iterate depth-averaging code
sName = 'depthAvg_6mDrogue_15mCenter_%s_'%(whatYear) #name to save files

#test to see if file exists (to avoid overwriting completed files)
sName_full_u = saveDir+sName+'u.nc'
sName_full_v = saveDir+sName+'v.nc'
if (os.path.exists(sName_full_u) or os.path.exists(sName_full_v)):
    assert False, 'file exists for %s'%(whatYear)

#start and end date
start = dt.datetime(whatYear,1,1,0,0)
end = dt.datetime(whatYear,12,31,0,0)

#file names to combine in final step:
uFileNameList = []
vFileNameList = []

#depth average across drogue length for each day in whatYear:
for date in pd.date_range(start,end):
    #names for averaged u and v files
    uFname = 'temp_data/uAvg_15mCent_6mLen_%s_%s_%s.nc'%(date.year,date.month,date.day)
    vFname = 'temp_data/vAvg_15mCent_6mLen_%s_%s_%s.nc'%(date.year,date.month,date.day)

    #get name for daily u and v files from Mercator data (thisDay is an older naming convention)
    thisDay=mdates.date2num(date)- mdates.date2num(np.datetime64('0000-12-31')) 
    uDataName = dataDir+'umerc_phy_%0.0d.nc'%(thisDay,)
    vDataName = dataDir+'vmerc_phy_%0.0d.nc'%(thisDay,)

    #open u and v files; get velocities
    uField = xr.open_dataset(uDataName)
    vField = xr.open_dataset(vDataName)
    uF = uField['vozocrtx'] #u velocity (zo is for zonal)
    vF = vField['vomecrty'] #v velocity (me is for meridional)
    
    #get coords for dataArray
    xArr = uF['x'].values
    yArr = uF['y'].values
    time_ct = uF['time_counter'].values
    dpth = [15.]

    #latitude and longitude of u and v points
    vnLon = vF['nav_lon'].values
    vnLat = vF['nav_lat'].values
    unLon = uF['nav_lon'].values
    unLat = uF['nav_lat'].values
    
    #interpolate velocities across drogue length
    interp_dep = uField['deptht'][:].values
    iDep = interp_dep[8:12]
    u_used = uField['vozocrtx'][0,8:12,:,:].values
    v_used = vField['vomecrty'][0,8:12,:,:].values
    int_fun_u = interpolate.interp1d(iDep,u_used,axis=0)
    int_fun_v = interpolate.interp1d(iDep,v_used,axis=0)#faster than xarray's inbuilt interpolation
    
    #6m drogue length - based on Pacific Gyre manufactured drogues
    #interpolate at top, bottom, and middle of drogue
    uTop = int_fun_u(12)
    uMid = int_fun_u(15)
    uBot = int_fun_u(18)
    
    vTop = int_fun_v(12)
    vMid = int_fun_v(15)
    vBot = int_fun_v(18)
    
    #compute average velocities across drogue depth and add dims for DataArray
    uAvg = np.expand_dims((uTop+uMid+uBot)/3.0,axis=(0,1))
    vAvg = np.expand_dims((vTop+vMid+vBot)/3.0,axis=(0,1))

    #create DataArrays and append
    uDs = xr.DataArray(data=uAvg,name = 'uAvg',dims=['time','depth','y','x'],
                       coords=dict(u_lon=(['y','x'],unLon),
                                   u_lat=(['y','x'],unLat),
                                   time=time_ct,
                                   depth=dpth))
    vDs = xr.DataArray(data=vAvg,name = 'vAvg',dims=['time','depth','y','x'],
                       coords=dict(v_lon=(['y','x'],vnLon),
                                   v_lat=(['y','x'],vnLat),
                                   time=time_ct,
                                   depth=dpth))
    
    uFileNameList.append(saveDir+uFname)
    vFileNameList.append(saveDir+vFname)

    #save individual u and v files as netcdf
    uDs.to_netcdf(saveDir+uFname)
    vDs.to_netcdf(saveDir+vFname)
    print('%s done in %s'%(date, time.time()-t1))


#now combine netcdf files for each year and add attributes

#create attribute dicts:
at_institution = 'University of New Hampshire'
at_source = 'Original model data from the Mercator Ocean 1/12 degree phsyical model PSY4V3R1 on the native model grid'
at_comment = 'contact: wl1039@wildcats.unh.edu'
at_units = 'm s-1'

u_at_title = 'Depth-averaged u-velocity (6m drogue length centeredt at 15m)'
u_at_long_name = 'Zonal Velocity'
u_at_standard_name = 'sea_water_x_velocity'
u_at_short_name = 'uAvg'

v_at_title = 'Depth-averaged v-velocity (6m drogue length centeredt at 15m)'
v_at_long_name = 'Zonal Velocity'
v_at_standard_name = 'sea_water_x_velocity'
v_at_short_name = 'vAvg'

u_attrs = dict(title = u_at_title,
               institution = at_institution,
               source = at_source,
               comment = at_comment,
               units = at_units,
               short_name = u_at_short_name,
               long_name = u_at_long_name,
               standard_name = u_at_standard_name)
v_attrs = dict(title = v_at_title,
               institution = at_institution,
               source = at_source,
               comment = at_comment,
               units = at_units,
               short_name = v_at_short_name,
               long_name = v_at_long_name,
               standard_name = v_at_standard_name)

#open all daily .nc files as single xarray dataset and save as .nc file
print('making combined files...')
print('opening u files as single dataset for %s'%(whatYear))
u_bigFile = xr.open_mfdataset(uFileNameList,parallel=True)
print('done opening u dataset')
u_bigFile.attrs = u_attrs
print('saving u values to netcdf for %s'%(whatYear))
write_job = u_bigFile.to_netcdf(sName_full_u,compute=False) #save

with ProgressBar():
    print(f"Writing to %s"%(sName_full_u))
    write_job.compute()

u_bigFile.close()

print('opening v files as single dataset for %s'%(whatYear))
v_bigFile = xr.open_mfdataset(vFileNameList,parallel=True)
print('done opening v dataset')
v_bigFile.attrs = v_attrs
print('saving v values to netcdf for %s'%(whatYear))
write_job = v_bigFile.to_netcdf(sName_full_v) #save

with ProgressBar():
    print(f"Writing to %s"%(sName_full_v))
    write_job.compute()

v_bigFile.close()

print('Done with year %s'%(whatYear))
