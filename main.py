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
    
    resolution = 30                                      
    filename = "sub_190_yO.hdf5"       # enter the filename of the desired bandstructure data here
    
    geom = a_poly_tapered_cavity(substrate_parameters= {'waveguide_height'  : 0.19, 'substrate_height'  : 5},
                                 waveguide_parameters = {'wvg_width' : 0.7, 'wvg_height'  : 0.19},
                                 filename = filename) 
    
    
    
    boundary_layers = get_boundary_layer(sim2d=False)
    
    fcen = 1/1.54
    df = 0.2
    
    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), 
                         component=mp.Ey, 
                         center=mp.Vector3())]
    
    symmetries = [mp.Mirror(mp.X, +1)]   
    
    sim = mp.Simulation(resolution = resolution, 
                        cell_size=mp.Vector3(20, 8, 8), 
                        geometry=geom, 
                        boundary_layers=boundary_layers, 
                        sources=sources,
                        symmetries = symmetries,
                        progress_interval=100, dimensions = 3)
    
    h = mp.Harminv(mp.Ey, mp.Vector3(0, 0, 0), fcen, df)
    time_after_source = 500
    
    #mp.when_true(lambda x: sim.round_time()%50 == 0 and sim.round_time()> 200, lambda x: print(h.modes)
    
    sim.run(mp.at_beginning(mp.output_epsilon),
            mp.after_sources(h),
            until_after_sources=time_after_source)

    sim.run(mp.at_every(1/(fcen/5), mp.output_efield), until=1/fcen)  # if you want to generate the h5 files of the field values,
    
    visualize_sim_cell(sim)
    
    print("Modal Volume: {}".format(sim.modal_volume_in_box()))
 
    
if __name__ == '__main__':
    main()
