import numpy as np
import giapy as gp

# American Northeast map parameters
amne = {'area_tresh': 1000.0,
        'llcrnrlat': 35,
        'llcrnrlon': -85,
        'urcrnrlat': 50,
        'urcrnrlon': -60,
        'projection': 'merc',
        'resolution': 'i',
        'rsphere': (6378137.0, 6356752.3142)}


def compute_gia(giasimconfig=None, out_times=None):
    """Compute the glacial isostatic adjustment.
    
    Computes it by thousand year increments, with a point far in future (30ka).
    """
    giasimconfig = {'ice': 'AA2_Tail_nochange5_hightres_Pers_288_square',
                'earth': '75km0p04Asth_4e23Lith',
                'topo': 'sstopo288'}

    sim = gp.sim.configure_giasim(giasimconfig)
    result = sim.performConvolution(eliter=5,
                                out_times=np.union1d(np.arange(25), [-30, -0.1, 0.1]))

    return result

def comp_diff_emerge(result, loc1, loc2, time, obs='sstopo'):
    """Compute differential glacial isostatic adjustment between locs 1 and 2.

    Parameters
    ----------
    result : GiaSimOutput
    loc1, loc2 : tuples of (lon, lat)
    time : the time for comparison
    obs : str
        The type of observation, defaults to 'sstopo', which is emergence. See
        print(result) for all possible fields. Note: fields computed in
        spectral space must be transformed prior to comparison, use
        result.transformObservers().
    """
    obs = getattr(result, obs)
    field = obs.nearest_to(time).T
    interper = result.inputs.grid.create_interper(field)
    xx, yy = result.inputs.grid.basemap((loc1[0], loc2[0]), (loc1[1], loc2[1]))
    diffs = interper.ev(xx, yy)
    return diffs[1] - diffs[0]
    
# Locations
ithaca = (-76.5, 42.443333)
rome = (-75.5, 43.2)
lake_hall = ((-77, -11, -37), (42, 51, 25)), 15.5, 'Lake Hall'
lake_warren = ((-77, -40), (42, 32)), 14, 'Lake Warren'

if __name__ == '__main__':

    import os, fnmatch
    import cPickle as pickle
    from scipy.interpolate import RectBivariateSpline

    # Check if the results are existent or already loaded.
    # If neither, compute them and save.
    def find(pattern, path):
        result = []
        for root, dirs, files in os.walk(path):
            for name in files:
                if fnmatch.fnmatch(name, pattern):
                    result.append(os.path.join(root,name))
        return result
    files = find('*.giaresult', './data/')
    try:
        result.sstopo
    except NameError:
        if not files:
            result = compute_gia()
            pickle.dump(result,
                    open('./data/'+gp.timestamp()[:10]+'.giaresult', 'w'))
        else:
            result = pickle.load(open(files[0], 'r'))

    # Load a higher-resolution topography
    topo576 = np.loadtxt('./data/576Topo.txt', skiprows=1, delimiter=',')
    topo576 = topo576.T.reshape(-1, 576, 576)
    topo576interp = RectBivariateSpline(topo576[0,0,:], 
                                        topo576[1,:,0],
                                        topo576[2].T)
    
    # Transform uplift, etc., fields
    result.transformObservers()

    lakestr = '{} at {} ka was at elevation {} m, {} m relative to Ithaca'
    uplstr = 'Height of paleosurface above present day sea level: {}.'

    for loc in [lake_warren, lake_hall]:
        dd = gp.maps.dms2dd(loc[0][0]), gp.maps.dms2dd(loc[0][1])
        sstopo = result.sstopo.nearest_to(loc[1])
        elev = result.inputs.grid.interp(sstopo.T, dd[0], dd[1])
        diff = comp_diff_emerge(result, dd, ithaca, loc[1])
        print(lakestr.format(loc[2], loc[1], elev, diff))
        
        topo = topo576interp.ev(dd[0], dd[1])
        upl  = result.inputs.grid.interp(result.upl.nearest_to(loc[1]).T,
                                         dd[0], dd[1])
        print(uplstr.format(upl+topo))
                                        

