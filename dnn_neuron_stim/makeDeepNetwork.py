import tensorflow
from tensorflow.keras.layers import Dense, Flatten, Conv3D, BatchNormalization,\
                                    Activation,Dropout
from tensorflow.keras.models import Sequential
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import load_model
import os
import json, pprint
def makeDNN(x_train_shape = 81, y_len = 3, # defined by data structure
            num_layers = 4, nodes_per_layer = 50, lr = 5e-08, # optional model architecture settings
            activation = 'relu', shrink_rate = .75, dropout_rate = 0.5,
            output_activation = 'linear', batch_norm_on = 1):

    if activation == 'leaky_relu':
        activation = tensorflow.nn.leaky_relu

    model = Sequential()
    #build the layers
    for i in range(num_layers): # Adds num_layers dense layers, followed by output layer of length y_len
        if i == 0: # first dense layer
            model.add(Dense(nodes_per_layer,input_shape=(x_train_shape,)))
        elif i <= num_layers - 1: # all intermediate layers and last layer
            model.add(Dense(max(int(nodes_per_layer*(shrink_rate**i)),y_len), kernel_initializer='normal')) #max taken to prevent ensure num_neurons > y_len
        # After all dense layers add BatchNormalization, Activation, and Dropout layers
        if batch_norm_on: # add batch norm before activation
            model.add(BatchNormalization(center=True, scale=True))
        model.add(Activation(activation))
        if batch_norm_on: # or, add batch norm after activation
            model.add(BatchNormalization(center=True, scale=True))
        model.add(Dropout(dropout_rate))

    model.add(Dense(y_len,activation='linear')) # Add layer for outputs

    model.compile(loss='mse',optimizer=Adam(lr=lr),metrics = ['mse'])
    return model

def make3DCNN(N_per_dim = 7, y_len = 3, num_channels = 3, # defined by data structure
              num_conv = 4,num_dense = 2,nodes_per_layer = 50, # optional model architecture settings
              filters = 80,ker_size = 2,lr = 5e-08, activation = 'relu',
              shrink_rate = .75, dropout_rate = 0.5, output_activation = 'linear',
              batch_norm_on = 1,kernel_initializer='glorot_uniform',kernel_regularizer=None,
              padding = 'valid'):

    if activation == 'leaky_relu':
        activation = tensorflow.nn.leaky_relu
    if output_activation == 'leaky_relu':
        output_activation = tensorflow.nn.leaky_relu

    model = Sequential()
    #build the layers
    for i in range(num_conv): # Adds num_conv convolutional layers, followed by flatten
        if i == 0: # first layer
            model.add(Conv3D(filters,kernel_size = ker_size, padding = padding,
                             input_shape=(N_per_dim,N_per_dim,N_per_dim,num_channels),
                             kernel_initializer = kernel_initializer,
                             kernel_regularizer = kernel_regularizer))
        elif i <= num_conv - 1: # all intermediate layers and last layer
            model.add(Conv3D(max(int(filters*(shrink_rate**i)),y_len),
                             kernel_size = ker_size, padding = padding,
                             kernel_initializer = kernel_initializer,
                             kernel_regularizer = kernel_regularizer))
        # After all convolutional layers, add BatchNormalization and Activation layers
        if batch_norm_on == 1: # add batch norm before activation
            model.add(BatchNormalization(center=True, scale=True))
        model.add(Activation(activation))
        if batch_norm_on == 2: # or, add batch norm after activation
            model.add(BatchNormalization(center=True, scale=True))

    model.add(Flatten()) # Add flatten to go from 3D conv layers to 1D dense layers

    for i in range(num_dense): # Adds num_dense dense layers, followed by output layer of length y_len
        if i == 0: # first dense layer
            model.add(Dense(nodes_per_layer, kernel_initializer = kernel_initializer,
                            kernel_regularizer = kernel_regularizer))
        elif i <= num_dense - 1: # all intermediate layers and last layer
            model.add(Dense(max(int(nodes_per_layer*(shrink_rate**i)),y_len),
                            kernel_initializer = kernel_initializer,
                            kernel_regularizer = kernel_regularizer))
        # After all dense layser add Batch Normalization, Activation, and Dropout layers
        if batch_norm_on == 1: # add batch norm before activation
            model.add(BatchNormalization(center=True, scale=True))
        model.add(Activation(activation))
        if batch_norm_on == 2: # or, add batch norm before activation
            model.add(BatchNormalization(center=True, scale=True))
        if dropout_rate > 0:
            model.add(Dropout(dropout_rate))

    model.add(Dense(y_len,activation=output_activation)) # Add layer for outputs

    model.compile(loss='mse',optimizer=Adam(lr=lr),metrics = ['mse'])
    return model

def loadNet(weights_file_path):
    # Load model from weights file (must contain model as well,
    # i.e. save_weights_only = False)
    weights_dir = os.path.dirname(weights_file_path)
    weights_name, ext = os.path.splitext(os.path.basename(weights_file_path))
    if not ext: # No extension included in filename, assume .hdf5
        ext = '.hdf5'
    elif ext != '.hdf5': # catch case when filename has period
        ext = ext + '.hdf5'
        
    model = load_model(os.path.join(weights_dir,weights_name+ext),
                       custom_objects={'leaky_relu': tensorflow.nn.leaky_relu})
    return model

def writeNet(model,model_json_file_path):
    with open(model_json_file_path,'w') as f:
        f.write(model.to_json())


# 3 methods for loading models:
    # 1) Load arguments to make3DCNN/makeDNN in net_params field of json (model_format = 'json')
    # 2) Load model saved by keras model_to_json (makeDeepNetwork.writeNet())
    #    with model_from_json (model_format = 'json', automatically detected) - can include pretrained weights
    # 3) Load model from hdf5 file output by keras model.save() (model_format = 'hdf5') - should be identical to 1
def loadModelFromFile(model_params_name,model_format,x_shape,y_shape,
                      E_mode='3D',params_dir='data/model_params/',weights_dir='data/weights/'):

    if model_format == 'json':
        model_params_file = os.path.join(params_dir,model_params_name + '.json')
        with open(model_params_file) as f:
            params = json.load(f)

        if 'net_params' in params:
            net_params = params['net_params']
            print('Loaded net_params from {} file {}'.format(model_format,model_params_file))
            pp = pprint.PrettyPrinter(indent=2)
            print('net_params:')
            pp.pprint(net_params)

            # Make Networks
            if '3D' in E_mode:
                # Check if input different grid dimension for CNN than E-field data (should be higher)
                if 'N_per_dim' in net_params:
                    N_per_dim = net_params['N_per_dim'] # use higher resolution
                    # remove from args and input below
                    del net_params['N_per_dim']
                else:
                    N_per_dim = x_shape[1] # use sample sampling grid as input E-field data
                print('Running make3DCNN...')
                model = make3DCNN(N_per_dim = N_per_dim,
                                  y_len = y_shape[1],
                                  num_channels = x_shape[4],
                                  **net_params)

            else:
                print('Running makeDNN...')
                model = makeDNN(x_train_shape=x_shape[1],
                                y_len=y_shape[1],
                                **net_params)
        else: # Assume saved as keras json file
            from tensorflow.keras.models import model_from_json
            print('Loaded model from keras formatted {} file {}'.format(model_format,model_params_file))
            model = model_from_json(model_params_file)

    elif model_format == 'hdf5': # Load model file output by keras model.save
        model_params_file = os.path.join(weights_dir,model_params_name + '.hdf5')
        model = loadNet(model_params_file)
        print('Loaded model from {} file {}'.format(model_format, model_params_file))

    return model

def saveFilters(weights_file):
    import scipy.io
    import numpy as np

    model = loadNet(weights_file)
    filters = []
    for layer in model.layers:
	    # check for convolutional layer
	    if 'conv' not in layer.name:
		    continue
	    # get filter weights
	    filt, biases = layer.get_weights()
	    filters.append(filt)
    filt_obj = np.empty((len(filters),), dtype=np.object)
    for i in range(len(filters)):
        filt_obj[i] = filters[i]
    name = weights_file.split('.hdf5')
    weights_name, ext = os.path.splitext(os.path.basename(weights_file))
    scipy.io.savemat('data/mat_data/'+ weights_name + '.mat', {'data': filt_obj})
