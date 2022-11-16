function python_exec_str = python_exec()
%PYTHON_EXEC Outputs string with to python executable 
%  
%   Inputs 
%   ------ 
%   Optional Inputs 
%   --------------- 
%   Outputs 
%   ------- 
%   Examples 
%   --------------- 

% AUTHOR    : Aman Aberra 
python_exec_str = ''; % replace with path to python executable
if isempty(python_exec_str)
    error('Specify path to python executable in ''python_exec_str'' within dnn_neuron_stim/mat_util/python_exec.m')
end
if ~exist(python_exec_str,'file')
    error('Incorrect path to python executable, no file found at %s',python_exec_str);
end