import numpy as np
# try:
#     from hdf5storage import loadmat # use for loading v7.3 mat files (>2 GB)
# except ImportError:
#     from scipy.io import loadmat
from scipy.io import loadmat
from scipy.interpolate import RegularGridInterpolator
import os
from datetime import datetime
def loadMatlabData(data_file_path): # WARNING: do not include periods in file name
    # Load .mat data file from absolute file path name
    data_dir = os.path.dirname(data_file_path)
    if data_dir == "": # file in current dir
        data_dir = os.getcwd()
    dataset_name, ext = os.path.splitext(os.path.basename(data_file_path))
    if not ext: # No extension included in filename, assume .mat
        ext = '.mat'
    elif ext != '.mat': # try appending '.mat' to end of full file name
        dataset_name = dataset_name+ext
        ext = '.mat'
    print("Loading dataset from {}".format(data_dir))
    try:  # include for older datafiles, obsolete now
        start_time = datetime.now()
        data = loadmat(os.path.join(data_dir,dataset_name+ext))
        cell_data = data['dataCell']

        efield = np.stack(cell_data[:,0],axis=0).squeeze()
        pol = np.stack(cell_data[:,1],axis=0).squeeze() # num_samples x 1 x num_funcs 3D tensor -> num_samples x num_funcs matrix

        # remove samples with nans (coordinate was outside GM/WM)
        Efield, NeuralResponse = remove_nan_inds(efield,pol)

        print("Loaded dataset using scipy.loadmat: {} (Time elapsed = {})".format(dataset_name,
                                                              datetime.now()-start_time))
    except:
        print('Using h5py to load v7.3 .mat file')
        import h5py # load v7.3 matlab files (>2 GB)
        start_time = datetime.now()
        f = h5py.File(os.path.join(data_dir,dataset_name+ext),'r')
        # dataCell = f['dataCell'] # N x 2 cell array {Efield,Pol}
        Eref = f['Efield']
        if 'Pol' in f.keys():
            Pref = f['Pol'] # Polarization (subthreshold stim)
        elif 'Prob' in f.keys():
            Pref = f['Prob'] # Firing probability (suprathreshold stim)
        elif 'Thresh' in f.keys():
            Pref = f['Thresh'] # Threshold (prob normalized by |E| center)
        # nrows = len(dataCell[0])

        Efield = np.zeros(Eref.shape)
        NeuralResponse = np.zeros(Pref.shape)
        Eref.read_direct(Efield)
        Pref.read_direct(NeuralResponse)

        Efield = Efield.transpose()
        NeuralResponse = NeuralResponse.transpose() # num_samples x 1 x num_funcs 3D tensor -> num_samples x num_funcs matrix

        Efield, NeuralResponse = remove_nan_inds(Efield, NeuralResponse)
        print("Loaded dataset using h5py: {} (Time elapsed = {})".format(dataset_name,
                                                              datetime.now()-start_time))

    return Efield, NeuralResponse

def loadEfield(data_file_path):
    # Load .mat data file from absolute file path name - Efield only
    data_dir = os.path.dirname(data_file_path)
    if data_dir == "": # file in current dir
        data_dir = os.getcwd()
    dataset_name, ext = os.path.splitext(os.path.basename(data_file_path))
    if not ext: # No extension included in filename, assume .mat
        ext = '.mat'

    print('Using h5py to load v7.3 .mat file')
    import h5py # load v7.3 matlab files (>2 GB)
    start_time = datetime.now()
    f = h5py.File(os.path.join(data_dir,dataset_name+ext),'r')
    # dataCell = f['dataCell'] # N x 2 cell array {Efield,Pol}
    Eref = f['Efield']

    Efield = np.zeros(Eref.shape)
    Eref.read_direct(Efield)

    Efield = Efield.transpose()

    Efield = remove_nan_inds1(Efield)
    print("Loaded dataset using h5py: {} (Time elapsed = {})".format(dataset_name,
                                                          datetime.now()-start_time))

    return Efield
# remove elements with nans from numpy.ndarrays x and y
def remove_nan_inds(x,y):

    x_inds = tuple([i for i in range(1,x.ndim)])
    y_inds = tuple([i for i in range(1,y.ndim)])
    nan_inds = np.logical_or(np.any(np.isnan(x),axis=x_inds),np.any(np.isnan(y),axis=y_inds))
    x = x[np.logical_not(nan_inds),:]
    y = y[np.logical_not(nan_inds),:]

    return x, y

# remove elements with nans from single numpy.ndarray x
def remove_nan_inds1(x):

    x_inds = tuple([i for i in range(1,x.ndim)])
    nan_inds = np.any(np.isnan(x),axis=x_inds)
    x = x[np.logical_not(nan_inds),:]
    return x

def upsampleEfield(Efield,model_shape):
    # Upsample Efield sampled on 3D grid, assumes cubic dimensions
    N1 = Efield.shape[1] # N per dimensions of E-field data
    N2 = model_shape[1] # N per dimensions of 3D CNN

    if N2 > N1 and N2 % 2 == 1 and Efield.ndim == 5: # must be higher resolution, odd, and 3D grid (n samples x N x N x N x n_channels)        
        n_channels = Efield.shape[-1]        
        x1, y1, z1 = np.linspace(1,N1,N1),np.linspace(1,N1,N1),np.linspace(1,N1,N1) # original grid
        x2, y2, z2 = np.linspace(1,N1,N2),np.linspace(1,N1,N2),np.linspace(1,N1,N2) # upsampled grid
        X2,Y2,Z2 = np.meshgrid(x2,y2,z2,indexing='xy') # rows Y, columns X, matches matlab
        # N2^3 x 3 array of grid points (x,y,z) for sampling interpolator
        pts2 = np.concatenate((np.atleast_2d(X2.flatten('F')).T,np.atleast_2d(Y2.flatten('F')).T,np.atleast_2d(Z2.flatten('F')).T),axis=1)
        Efield2 = np.zeros((Efield.shape[0],N2,N2,N2,n_channels)) # upsampled Efield samples
        for c in range(n_channels):
            for i in range(Efield.shape[0]):
                Einti = RegularGridInterpolator((x1,y1,z1),Efield[i,:,:,:,c],method='linear')
                Efield2[i,:,:,:,c] = Einti(pts2).reshape(N2,N2,N2,order='F').transpose((1,0,2)) # transpase to same ordering as 3D arrays output from matlab                
        print('Upsampled E-field input from N={} to N={}, new dimensions: {}'.format(N1,N2,Efield2.shape[1:]))
        return Efield2    
    elif N2 == N1:
        return Efield # nothing to be done
    else:
        raise NotImplementedError('Using CNN resolution (N={}) lower than E-field input (N={}) not implemented yet'.format(
                                                  N2,N1))