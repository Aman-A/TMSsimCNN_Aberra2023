function data_out = apply_func(data_in,func)
% APPLY_FUNC apply function to array of data across replicates to collapse
% into vector
%  
%   Inputs 
%   ------
%   data_in : num_samples x num_replicates array
%             rows are different samples (e.g. cell positions for neural
%             simulation data), columns are replicates to apply function
%             over (e.g. clones/rotations for neural simulation data)
%   func : string or numeric
%          Function to apply to compartment data from each simulation,
%          either: 'max', 'maxabs', 'mean', 'median','min','maxabs_ind'.
%          'maxabs_ind' extracts value from each simulation at the single 
%          compartment with the maxabs value across all simulations. If
%          numeric and < 1, calculates quantile of func. Input 'none' if
%          you just want to extract data in comp_type/sectype compartments
%          without applying a function to collapse data across compartments
if ischar(func)
    if strcmp(func,'max')
        data_out = max(data_in,[],2,'omitnan');
    elseif strcmp(func,'maxabs')
        [~,max_inds] = max(abs(data_in),[],2,'linear','omitnan');
        data_out = data_in(max_inds);
    elseif strcmp(func,'mean')
        data_out = mean(data_in,2,'omitnan');
    elseif strcmp(func,'median')
        data_out = median(data_in,2,'omitnan');
    elseif strcmp(func,'min')
        data_out = min(data_in,[],2,'omitnan');
    elseif strcmp(func,'maxabs_ind')
        [~,max_ind] = max(max(abs(data_in),[],1,'linear','omitnan'));
        data_out = data_in(:,max_ind);
        inds_i = find(inds); sec_ind = inds_i(max_ind);
        fprintf('Extracting data from %s\n',cell_datai.secnames{sec_ind});
    elseif strcmp(func,'max_min_logratio') % Ratio of max value to min abs value
        data_max = max(data_in,[],2,'omitnan');
        data_min = min(data_in,[],2,'omitnan');        
        data_out = log10(abs(data_max)./abs(data_min));
    elseif strcmp(func,'std')
        data_out = std(data_in,0,2,'omitnan');
    elseif strcmp(func,'range')
        data_out = range(data_in,2);
    elseif strcmp(func,'none')
        data_out = data_in; 
    else
        error('func: %s does not exist\n',func); 
    end
elseif isnumeric(func) && func < 1 % quantile
    data_out = quantile(data_in,func,2);
elseif isinteger(func) % func is index to extract
    data_out = data_in(:,func);
end