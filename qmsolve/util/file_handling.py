# used for loading and storing eigenstates in large computations. Requires h5py package
import numpy as np
from ..eigenstates import Eigenstates
import h5py

def save_eigenstates(eigenstates, name:str):
    #create dataset
    with h5py.File(name, 'a') as f:
        dset = f.create_dataset('energies', (eigenstates.energies.shape),dtype='float64')
        f['energies'][:] = eigenstates.energies
        dset = f.create_dataset('array', (eigenstates.array.shape),dtype='float64')
        f['array'][:] = eigenstates.array
        f.attrs['number'] = eigenstates.number
        f.attrs['extent'] = eigenstates.extent
        f.attrs['N'] = eigenstates.N
        f.attrs['type'] = eigenstates.type

def load_eigenstates(name:str):
    with h5py.File(name, 'a') as f:
        energies = np.copy(f['energies'])
        array = np.copy(f['array'])
        extent = f.attrs['number']
        number = f.attrs['extent']
        N = f.attrs['N']
        type = f.attrs['type']
    return Eigenstates(energies, array, extent, N, type) 
