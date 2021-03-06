import meep as mp
import warnings

from ..utilities.sweep_util import *


def add_waveguide_1d(geom=None, wvg_width=.65, wvg_height=.25, center=mp.Vector3(0, 0, 0),
                     material=mp.Medium(index=3.45), d_tuning=0, material_tuning=mp.Medium(index=1.025)):
    """
    Creates a 1D waveguide.
    The LHe shell can be adopted such that N2 ice layers could be investigated
    as well.
    """

    if geom is None:
        geom = []
    if isinstance(center, mp.Vector3):
        _center = center
    elif len(center) == 3:
        _center = mp.Vector3(center[0], center[1], center[2])
    else:
        warnings.warn("Variable center not understood but passed")
        _center = center

    # LHe shell
    if d_tuning != 0:
        geom.append(mp.Block(material=index_to_material(material_tuning),
                             center=_center,
                             size=mp.Vector3(mp.inf, wvg_width + 2 * d_tuning, wvg_height + 2 * d_tuning)))

    # Si waveguide
    geom.append(mp.Block(material=index_to_material(material),
                         center=_center,
                         size=mp.Vector3(mp.inf, wvg_width, wvg_height)))

    return geom


def add_waveguide_1d_on_substrate(geom=None, wvg_height=.22, substrate_height=5, substrate_material=1.44, **kwargs):

    geom = add_waveguide_1d(geom=geom, wvg_height=wvg_height, **kwargs)

    geom = add_substrate(geom=geom, substrate_height=substrate_height, material=substrate_material)

    return geom


def add_substrate(geom=None, wvg_height=0.22, substrate_height=5, material=mp.Medium(index=1.44),
                  embed_in_substrate = True):
    """
    Creates a (by default SiO2) substrate.
    If the unlike case occurs that we need to change the center, we could again make everything relative to the center
    of the substrate.
    """
    if geom is None:
        geom = []

    _center = mp.Vector3(0, 0, -wvg_height/2 - substrate_height/2)

    geom.append(mp.Block(material=index_to_material(material),
                         center=_center,
                         size=mp.Vector3(mp.inf, mp.inf, substrate_height)))

    return geom

def add_filled_substrate(geom=None, wvg_height=0.22, material=mp.Medium(index=1.44), substrate_height=5,
                         embed_in_substrate = True, sim_shape=(20, 8, 8)):
    
    """
    Return geom after completely filling its surroundings with the substrate material. 
    If geom is empty, than just returns a cell sized (sx, sy, sz) filled with the background material.
    """
    sx, sy, sz = sim_shape
    geom_temp = 0
                         
    if geom is not None:
        geom.pop()
        geom_temp = geom
        
   
    substrate_geom = [mp.Block(material=index_to_material(material),
                               center=mp.Vector3(0, 0, 0),
                               size=mp.Vector3(sx, sy, sz))]
    substrate_geom.extend(geom_temp)
    
    geom = substrate_geom
    
    center_below = mp.Vector3(0, 0, -wvg_height/2 - substrate_height/2)
    geom.append(mp.Block(material=index_to_material(material),
                         center=center_below,
                         size=mp.Vector3(mp.inf, mp.inf, substrate_height)))
    
    center_up = mp.Vector3(0, 0, wvg_height/2 + substrate_height/2)
    geom.append(mp.Block(material=index_to_material(material),
                         center=center_up,
                         size=mp.Vector3(mp.inf, mp.inf, substrate_height)))
    
    return geom

