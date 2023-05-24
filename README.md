# TMSsimCNN_Aberra2023
Code associated with Aberra AS, Lopez A, Grill WM, Peterchev AV. (2023). "Rapid estimation of cortical neuron activation thresholds by transcranial magnetic stimulation using convolutional neural networks".  Neuroimage

1) Install python dependencies (recommended to install into a virtual environment)
2)  Add path to python executable to `python_exec.m`
3) Download test dataset from doi.org/10.5281/zenodo.7326394 (~5.6 GB)
4) Run init.m to set up paths to dependencies (set download_test_dataset = 0 if step 2 is skipped)
5) OPTIONAL: To interpolate E-field vectors at all sampling grids from test dataset, run interpEfieldSample_all.m. Note: Running interpolation for all 15 neuron models, 12 rotations, and ~5,000 positions serially on a single CPU would take several weeks and >40 GB storage, parallelizing on high-performance computing cluster strongly recommended. `interpEfieldSample.m` uses all available CPUs to parallelize within model across all positions. Interpolated Efields for all neuron models at a single azimuthal rotation included in the downloadable dataset (found in `dnn_neuron_stim/output_data/nrn_efields/layer_set_1/M1_PA_MagVenture_MC_B70_ernie/nrn_pop1_maxH`)

To generate Fig. 3, Supp. Fig 1, and Supp. Fig 3, see `plot_Fig3_panels.m`, `plot_SuppFig1.m`, and `plot_SuppFig3.m`, respectively in `plot_manuscript_figs/`

To run threshold estimation on single neuron model/position, see `plot_manuscript_figs/run_estimation_single_neuron.m` script

Dependencies:

Python: 
    keras/tensorflow
    h5py
    scipy
    numpy

MATLAB toolboxes:
    Parallel Computing Toolbox  
    version >r2022a, specifically for:
        `exportgraphics.m `
        >2020b:
        for `max`/`mean`/`min`/`median` functions with 'omitnan' flag (vs. nanmax)

FEM E-field simulation was conducted in SimNIBS v3.1 (https://simnibs.github.io/simnibs/build/html/index.html). The necessary MATLAB library functions to load/save/process output files are included in this repository (written by Andre Antunes and Axel Thielscher). Original code can be found at the SimNIBS github repository: https://github.com/simnibs/simnibs/tree/master/simnibs/matlab 

NEURON model data derived from models adapted from Blue Brain Project as part of Aberra AS, Peterchev AV & Grill WM (2018). Biophysically Realistic Neuron Models for Simulation of Cortical Stimulation. Journal of Neural Engineering 15, 066023. Original model code can be found at https://senselab.med.yale.edu/ModelDB/ShowModel?model=241165#tabs-1

  