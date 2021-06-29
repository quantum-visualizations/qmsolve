# used for loading and storing eigenstates in large computations. Requires h5py package

def save(eigenstates, name:str):
    import h5py
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

def load(eigenstates, name:str):
    import h5py
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
