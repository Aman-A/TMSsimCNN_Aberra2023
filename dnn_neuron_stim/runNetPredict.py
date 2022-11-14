#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 21:02:42 2021

@author: amanaberra
"""
from makeDeepNetwork import loadNet 
from loadMatlabData import loadEfield
from datetime import datetime
import os, sys
import h5py # load v7.3 matlab files (>2 GB)
def main(weights_file_path,Efield_file_path):
    # Load model from weights file, run model.predict on Efield from file, and save to file
    model = loadNet(weights_file_path)
    # Load Efield from temp file    
    Efield = loadEfield(Efield_file_path)
    start_time = datetime.now()
    NeuralResponse = model.predict(Efield)
    print('Ran model prediction on {} samples in {}'.format(Efield.shape[0],
                                                            datetime.now()-start_time))
    out_dir = os.path.dirname(Efield_file_path)    
    # Get temp file name to save output
    out_name, ext = os.path.splitext(os.path.basename(Efield_file_path))    
    out_file = os.path.join(out_dir,out_name + '_response' + ext)
    with h5py.File(out_file,'w') as f:
        f.create_dataset('NeuralResponse',NeuralResponse.shape,data=NeuralResponse)


if __name__ == '__main__':
    weights_file_path = sys.argv[1]
    Efield_file_path = sys.argv[2]
    main(weights_file_path,Efield_file_path)