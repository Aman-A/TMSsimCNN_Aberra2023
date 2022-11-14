# TMSsimCNN_Aberra2022
Code associated with Aberra AS, Lopez A, Grill WM, Peterchev AV. (2022). "Rapid estimation of cortical neuron activation thresholds by transcranial magnetic stimulation using convolutional neural networks".  bioRxiv

Startup:
1) Run init.m to set up paths to dependencies
2) Run cnn_neuron_stim/addPaths_dnn_neuron_stim
3) Create python environment with necessary packages from requirements.txt
OR install dependencies into existing environment and specify path in ...
4) Download data:
    -neuron population data for nrn_model_ver = 'maxH', nrn_pop1_maxH.mat - nrn_pop12_maxH.mat go in dnn_neuron_stim/output_data/layer_data/<mesh_name>/<mesh_roi_name>/<layer_set_name>/<nrn_model_ver>
    -Interpolated E-field vectors for test dataset (5 clones of L2/3 PCs, L4 LBCs, L5 PCs at 4,999-5,000 positions within respective layers) for nrn_pop1_maxH.mat - nrn_pop12_maxH.mat, go in dnn_neuron_stim/output_data/nrn_efields/layer_set_1/<Efield_name>/<nrn_pop_name>/
    -ernie mesh and solution data for test dataset (TMS of M1 handknob with PA directed induced current)

5) Optional: To interpolate E-field vectors at all sampling grids from test dataset, run interpEfieldSample_all.m. Note: Running interpolation for all 15 neuron models, 12 rotations, and ~5,000 positions serially on a single CPU would take several weeks, parallelization highly recommended. interpEfieldSample.m uses all available CPUs to parallelize within model across all positions. 
    
Dependencies:
Python: 
    keras/tensorflow
    h5py
    scipy
    numpy (check if installed with keras/tensorflow)
    # maybe?    
    matplotlib (check if come with keras/tensorflow)

MATLAB toolboxes:
    Parallel Computing Toolbox  
    version >r2022a, specifically for:
        exportgraphics.m 
        >2020b:
        max with 'omitnan' flag
SimNIBS MATLAB library:
    mesh_load_gmsh4
    mesh_get_tetrahedron_centers
    get_field_idx
    mesh_save_gmsh4