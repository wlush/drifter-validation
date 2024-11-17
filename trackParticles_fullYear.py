#Code to track all particles starting during a given year (uses particleTracking_core.py)

#import necessary libraries
import numpy as np
import particleTrackingCore as ptc

whichYear = 2022 #manually input year

saveDir = '/home/willlush/particleTracking_temp' #temporary directory for saving particle tracking results
print('starting ', whichYear) #print year
loadDir = '/data/break/willlush/drifter_validation/particle_start_positions/' #directory where startlist lives
loadName = loadDir+'startLocsForRun_10daySep_newBathy%s.npz'%(whichYear) #startlist name
loadLists = np.load(loadName,allow_pickle=True) #load startList

lon = loadLists['lon'] #start longitudes
lat = loadLists['lat'] #start latitudes
time = loadLists['time'] #start times
#ID = loadLists['ID']

tSort = np.argsort(time) #indices sorted by time 

lon = lon[tSort] #sort longitudes by time
lat = lat[tSort] #sort latitudes by time
time = time[tSort] #sort time by time
#ID = ID[tSort]

data = {'lon':lon, 'lat':lat, 'time':time,}#'ID':ID} #organize data into dict (possibly unnecessary)

ptc.track_particles(data['lon'],data['lat'],data['time'],60,6.0,saveDir+'drifterValidation_run_%s'%(whichYear)) #run particle tracking using particle tracking core code

print('done with tracking for %s'%(whichYear)) #note that run is finished
