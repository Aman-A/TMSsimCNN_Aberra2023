% Run this function first to add functions to path
function [main_dir,varargout] = addPaths_dnn_neuron_stim()
    main_dir =  [fileparts(which('addPaths_dnn_neuron_stim.m')) filesep 'dnn_neuron_stim'];
    pathCell = regexp(path, pathsep, 'split');
    test_dir = fullfile(main_dir,'mat_util'); % check if this subfolder is on path
    if ispc  % Windows is not case-sensitive
        onPath = any(strcmpi(test_dir, pathCell));
    else
        onPath = any(strcmp(test_dir, pathCell));
    end
    if ~onPath
        addpath([main_dir filesep '..'])
        addpath(genpath(main_dir)); % dont rerun if already on path       
    end
    if nargout == 2
        data_dir_file = fullfile(main_dir,'slurm','data_dir.txt');
        if exist(data_dir_file,'file')
            data_dir = strtrim(fileread(data_dir_file));
        else
            data_dir = main_dir;
        end
        varargout = {data_dir};
    end
end