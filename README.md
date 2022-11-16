# TMSsimCNN_Aberra2022
Code associated with Aberra AS, Lopez A, Grill WM, Peterchev AV. (2022). "Rapid estimation of cortical neuron activation thresholds by transcranial magnetic stimulation using convolutional neural networks".  bioRxiv

1) Install python dependencies (recommended to install into a virtual environment), and add path to python executable to python_exec.m
2) Run init.m to set up paths to dependencies and download test dataset 
3) OPTIONAL: To interpolate E-field vectors at all sampling grids from test dataset, run interpEfieldSample_all.m. Note: Running interpolation for all 15 neuron models, 12 rotations, and ~5,000 positions serially on a single CPU would take several weeks and >40 GB storage, parallelizing on high-performance computing cluster strongly recommended. interpEfieldSample.m uses all available CPUs to parallelize within model across all positions. 

To generate Fig. 3, Supp. Fig 1, and Supp. Fig 3, see plot_Fig3_panels.m, plot_SuppFig1.m, and plot_SuppFig3.m, respectively in plot_manuscript_figs/

To run threshold estimation on single neuron model/position, see plot_manuscript_figs/run_estimation_single_neuron.m script

Dependencies:

Python: 
    keras/tensorflow
    h5py
    SciPy
    Numpy

MATLAB toolboxes (user install):
    Parallel Computing Toolbox  
    version >r2022a, specifically for:
        exportgraphics.m 
        >2020b:
        for max/mean/min functions with 'omitnan' flag (vs. nanmax)

FEM E-field simulation was conducted in SimNIBS v3.1 (https://simnibs.github.io/simnibs/build/html/index.html). The necessary MATLAB library functions to load/save/process output files are included in this repository (written by Andre Antunes and Axel Thielscher). Original code can be found at the SimNIBS github repository: https://github.com/simnibs/simnibs/tree/master/simnibs/matlab 

NEURON model data derived from models adapted from Blue Brain Project as part of Aberra AS, Peterchev AV & Grill WM (2018). Biophysically Realistic Neuron Models for Simulation of Cortical Stimulation. Journal of Neural Engineering 15, 066023. Original model code can be found at https://senselab.med.yale.edu/ModelDB/ShowModel?model=241165#tabs-1

  