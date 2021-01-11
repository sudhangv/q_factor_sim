import h5py
import numpy as np

z_index = 124    # Position will be 120(z = 0) + (z_index-120)*1/res um
sim_id = 56
res = 30

filename = 'slice_{}.h5'.format(sim_id)
e_h5file = h5py.File('main-e-000700.00.h5','r')

ex_slice = e_h5file['ex'][:,:, z_index]
ey_slice = e_h5file['ey'][:,:, z_index]
ez_slice = e_h5file['ez'][:,:, z_index]

e_h5file.close()

with h5py.File(filename,'w') as f:
    f.create_dataset('ex', data=ex_slice)
    f.create_dataset('ey', data=ey_slice)
    f.create_dataset('ez', data=ez_slice)

    f.attrs['z_index'] = z_index
    f.attrs['sim_id'] = sim_id
    f.attrs['res'] = res 
