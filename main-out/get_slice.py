import h5py
import numpy as np

e_h5file = h5py.File('main-e-000600.00.h5','r')

ex_slice = e_h5file['ex'][:,:,123]
ey_slice = e_h5file['ey'][:, 123]
ez_slice = e_h5file['ez'][:,:,123]


e_h5file.close()

with h5py.File('slice.h5','w') as f:
    f.create_dataset('ex', data=ex_slice)
    f.create_dataset('ey', data=ey_slice)
    f.create_dataset('ez', data=ez_slice)
    
