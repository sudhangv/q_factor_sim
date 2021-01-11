import h5py
import sys

filename = sys.argv[1]
print("Filename: ",filename)

with h5py.File(filename, 'r') as f:
	print("Resolution: ", f.attrs['res'])
	print("Z Index: ", f.attrs['z_index'])
	print("Sim ID: ", f.attrs['sim_id'])
