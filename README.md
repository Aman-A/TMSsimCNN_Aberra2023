# TMSsimCNN_Aberra2022
Code associated with Aberra AS, Lopez A, Grill WM, Peterchev AV. (2022). "Rapid estimation of cortical neuron activation thresholds by transcranial magnetic stimulation using convolutional neural networks".  bioRxiv

1) Run init.m to set up paths to dependencies
2) Run addPaths_dnn_neuron_stim to add matlab code (dnn_neuron_stim/) to path 
3) OPTIONAL: To interpolate E-field vectors at all sampling grids from test dataset, run interpEfieldSample_all.m. Note: Running interpolation for all 15 neuron models, 12 rotations, and ~5,000 positions serially on a single CPU would take several weeks, parallelization highly recommended. interpEfieldSample.m uses all available CPUs to parallelize within model across all positions. 

To generate Fig. 3, Supp. Fig 1, and Supp. Fig 3, see plot_Fig3_panels.m, plot_SuppFig1.m, and plot_SuppFig3.m, respectively in plot_manuscript_figs/

To run threshold estimation on single neuron model/position, see plot_manuscript_figs/run_estimation_single_neuron.m script

Dependencies:

Python (installed in new environment by init.m): 
    keras/tensorflow
    h5py
    scipy
    numpy (check if installed with keras/tensorflow)
    # maybe?    
    matplotlib (check if come with keras/tensorflow)

MATLAB toolboxes (user install):
    Parallel Computing Toolbox  
    version >r2022a, specifically for:
        exportgraphics.m 
        >2020b:
        max with 'omitnan' flag

FEM E-field simulation was conducted in SimNIBS v3.1. The necessary MATLAB library functions to load/save/process this file is included in this repository. 

NEURON model coordinates correspond to models adapted from Blue Brain Project as part of Aberra AS, Peterchev AV & Grill WM (2018). Biophysically Realistic Neuron Models for Simulation of Cortical Stimulation. Journal of Neural Engineering 15, 066023. Model code can be found at https://senselab.med.yale.edu/ModelDB/ShowModel?model=241165#tabs-1

  