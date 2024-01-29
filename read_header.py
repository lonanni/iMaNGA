import h5py

def read_header():
    """ Read various attributes from the header group. """
    f       = h5py.File('/home/nannil/iMaNGA/RefL0012N0188/snapshot_028_z000p000/snap_028_z000p000.0.hdf5', 'r')
    a       = f['Header'].attrs.get('Time')         # Scale factor.
    h       = f['Header'].attrs.get('HubbleParam')  # h.
    boxsize = f['Header'].attrs.get('BoxSize')      # L [Mph/h].
    f.close()

    return a, h, boxsize
