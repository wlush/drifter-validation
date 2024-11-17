#Code to organize data, input startlist, and define a function to use Parcels to track particles in depth-averaged velocity fields

import numpy as np
import pylab as p
import pandas as pd
import datetime as dtm
import xarray as xr
import time
from parcels import AdvectionRK4, FieldSet, JITParticle, Variable, ParticleFile, ParticleSet, ScipyParticle, ErrorCode

#Data location
dataDir = None #Path to depth-averaged data files

#grid file location
gridDir = None #path to horizontal grid file for NEMO grid (from Mercator)
gridFile=gridDir+'ext-PSY4V3R1_mesh_hgr.nc'

#velocity filename patterns
fnameU = dataDir+'depthAvg_6mDrogue_15mCenter_*_u.nc'
fnameV = dataDir+'depthAvg_6mDrogue_15mCenter_*_v.nc'

#set up data for fieldset for 2d data on NEMO grid
filenames = {'U': {'lon': gridFile,
                   'lat': gridFile,
                   'data': fnameU },
             'V': {'lon': gridFile,
                   'lat': gridFile,
                   'data': fnameV}}
variables = {'U': 'uAvg',
             'V': 'vAvg'}
dimensions = {'lon': 'glamf', 'lat': 'gphif','time': 'time'}

#core function to track particles using Parcels
def track_particles(release_lon, release_lat, release_time, duration, output_frequency, outFile):
    #get first release in startlist
    firstRelease = np.amin(release_time)
    print(firstRelease) #print the first release
    lastRelease = np.amax(release_time) #get final release time
    #make timedelta object between first and last releases, with added duration of tracking
    t_diff = pd.Timedelta(lastRelease-firstRelease).total_seconds()+pd.Timedelta(days=duration).total_seconds()
    
    #make fieldset from grid using Parcels
    print('making fieldset')
    fSet = FieldSet.from_nemo(filenames, variables, dimensions, chunksize = 'auto', allow_time_extrapolation=True)
    print('done with fieldset')
    
    #define a particle class where particles have a defined age
    class AgeParticle(JITParticle):
        age=Variable('age',initial=0.,dtype=np.float32)
        
    #add a kernel that kills the particle after durationDays
    fSet.add_constant('maxage',pd.Timedelta(days=duration).total_seconds())
    #define behavior to delete particles older than maxage
    def SampleAge(particle,fieldset,time):
        particle.age=particle.age+math.fabs(particle.dt)
        if particle.age>fieldset.maxage:
            particle.delete()
            
    #define particle deletion behavior
    def DeleteParticle(particle, fieldset, time):
        print('deleting...')
        particle.delete()

    #get release times in proper datetime64 format
    time_release = np.array([np.datetime64(x)for x in release_time])

    #create particleset
    print('making particleset')
    pSet=ParticleSet.from_list(fieldset=fSet,pclass=AgeParticle,lon=np.array(release_lon),
                               lat=np.array(release_lat),time=np.array(time_release))
    print('done with particleset')

    #create particle output file and chunking (for performance - machine dependent)
    outFile=pSet.ParticleFile(name=outFile,outputdt=dtm.timedelta(hours=output_frequency),chunks=(2200,40))

    #make SampleAge into a kernel object
    k_SampleAge=pSet.Kernel(SampleAge)

    #track particles using Parcels, if particle runs out of bounds, delete
    pSet.execute(AdvectionRK4+k_SampleAge,runtime=t_diff,dt=pd.Timedelta(minutes=60.0).total_seconds(),output_file=outFile,verbose_progress=True,recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
    print('done with execution in %ss'%(time.time()-t1))
    outFile.close()

if __name__=='__main__':
    #some sample test cases...
    if False:
        #sample input data...
        sLon = [-74.0,-74.0,-74.0,-74.0,-74.0]
        sLat = [39.0,39.0,39.0,39.0,39.0]
        rel_dates = [dtm.datetime(2015,12,25,5,30),dtm.datetime(2019,7,7,10,45),dtm.datetime(2019,9,8,10,45),dtm.datetime(2015,12,25,5,30),dtm.datetime(2008,12,25,5,30)]
        of_name = 'pTracking_test'
        out_freq = 6.0

        track_particles(sLon,sLat,rel_dates,15.0, out_freq, of_name)
        
    if False:
        num_samples = 1000
        sLon = np.full(num_samples,-74.0)
        sLat = np.full(num_samples,39.0)
        rel_dates = pd.date_range(end='2/17/2019',periods=num_samples)
        of_name = 'pTracking_test'
        out_freq = 12.0

        track_particles(sLon,sLat,rel_dates,15.0, out_freq, of_name)

    if False:
        num_samples = 100
        sLon = np.full(num_samples,-68.0)
        sLat = np.full(num_samples,43.0)
        sID = np.arange(num_samples)
        rel_dates = pd.date_range(end='10/17/2007',periods=num_samples)
        of_name = 'pTracking_test'
        out_freq = 12.0

        track_particles(sLon,sLat,rel_dates,15.0, out_freq, of_name)
