import numpy as np
import giapy as gp

# American Northeast map parameters
amne = {'area_tresh': 1000.0,
        'llcrnrlat': 35,
        'llcrnrlon': -85,
        'urcrnrlat': 50,
        'urcrnrlon': -60,
        'projection': 'merc',
        'resolution': 'l',
        'rsphere': (6378137.0, 6356752.3142)}

giasimconfig = {'ice': 'AA2_Tail_nochange5_hightres_Pers_288',
                'earth': '75km0p04Asth_4e23Lith',
                'topo': 'sstopo288'}

gp.sim.GiaSimGlobal()
sim.out_times = np.union1d(ice.times, [-30, -0.1, 0.1])
result = sim.performConvolution(topo=topo, eliter=5)
