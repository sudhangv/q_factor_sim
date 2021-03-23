import meep as mp
from meep import mpb
import numpy as np
import datetime
import math
import sys
import matplotlib.pyplot as plt
from IPython.display import Video
import warnings


from cavitysimulations.geometry.waveguide import *
from cavitysimulations.geometry.cavity import *
from cavitysimulations.visualization import *
from cavitysimulations.utilities.sweep_util import *


def main():
    
    wvg_height = 0.1
    resolution = 40
    sx, sy, sz =  20, 8, 8

    geom = a_poly_tapered_cavity(substrate_parameters= {'wvg_height': wvg_height, 'substrate_height': 5},
                                 waveguide_parameters = {'wvg_width': 0.7, 'wvg_height'  : wvg_height},
                                 cavity_parameters = {'hx' : 0.175, 'hy': 0.325, 'w' : 0.7, 'wvg_height': wvg_height,
                                                      'a_center': 0.429, 'a_mirror': 0.449, 
                                                      'Lx': 20, 'number_of_tapered_holes': 10,
                                                      'substrate': True, 'mode': 'yO' },
                                 sweep_parameters = {'hx_min_sweep': 0.05, 'hy_min_sweep': 0.05},
                                 filename = "bandstructure_data/perturb_sub_100_yO.hdf5")
    

    boundary_layers = get_boundary_layer(sim2d=False)
    
    fcen = 0.6412
    df = 0.005
    
    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), 
                         component=mp.Ey, 
                         center=mp.Vector3())]
    
    #symmetries = [mp.Mirror(mp.X,+1), mp.Mirror(mp.Y,-1), mp.Mirror(mp.Z,+1)]
    #symmetries = [mp.Mirror(mp.X, +1), mp.Mirror(mp.Y,-1)]
    symmetries = [mp.Mirror(mp.X, +1)]   
    
    sim = mp.Simulation(resolution = resolution, 
                        cell_size = mp.Vector3(sx, sy, sz), 
                        geometry = geom, 
                        boundary_layers = boundary_layers, 
                        sources = sources,
                        symmetries = symmetries,
                        progress_interval = 100, dimensions = 3)
    
    h = mp.Harminv(mp.Ey, mp.Vector3(0, 0, 0), fcen, df)
    time_after_source = 500

    #f = open("timed_harminv.txt", "w")
    
    sim.use_output_directory()

    sim.run(mp.at_beginning(mp.output_epsilon),
            mp.after_sources(h),
            mp.when_true(lambda x: sim.round_time()%50 == 0 and sim.round_time()> 200, lambda x: print("MODES: " + str([m.freq for m in h.modes])) ),
            until_after_sources=time_after_source)
    
    #f.close()

    sim.run(mp.at_every(1/fcen/3, mp.output_efield), until=1/fcen)

#    print("--------------------------------")
#    print([m.freq for m in h.modes])
#    print("---------------------------------")
    
    print("Modal Volume: {}".format(sim.modal_volume_in_box()))
    
    visualize_sim_cell(sim)
    
if __name__ == '__main__':
    main()
