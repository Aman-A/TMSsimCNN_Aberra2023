function nrn_pop_file = getNrnPopFileName(nrn_pop_name,nrn_model_ver,reposition_mode)
%GETNRNPOPNAME Output name of neuron population file given input settings 
%  
%   Inputs 
%   ------ 
%   nrn_pop_name : string
%   nrn_model_ver : string
%   reposition_mode : string
%                     'off' (default) or repositionNeuronLayer mode: 'both','cell_length','fixed', or 'd2layer' 
%   Optional Inputs 
%   --------------- 
%   Outputs 
%   ------- 
%   Examples 
%   --------------- 

% AUTHOR    : Aman Aberra 
if any(strcmp(reposition_mode,{'off','none'})) || isempty(reposition_mode) 
    % default file output by generateNeuronPop
    nrn_pop_file = sprintf('%s_%s',nrn_pop_name,nrn_model_ver);
else % nrn pop file after running replace (if over_write was 0)
    nrn_pop_file = sprintf('%s_%s_%s_r',nrn_pop_name,nrn_model_ver,reposition_mode);
    fprintf('Loading repositioned population: %s\n',nrn_pop_file); 
end