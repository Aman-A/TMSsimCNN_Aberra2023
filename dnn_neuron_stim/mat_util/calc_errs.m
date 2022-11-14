function results = calc_errs(y_actual,y_pred,axis)
if nargin < 3
    axis = 1; % default compute for each column
end
errs = y_pred - y_actual;
results.errs = errs; 
results.mae = mean(abs(errs),axis,'omitnan');  % mean absolute error
results.mse = mean(errs.^2,axis,'omitnan'); % mean square error
results.rmse = sqrt(results.mse); % root mean square err
results.nrmse = results.rmse./range(y_actual,axis); % norm root mean square error (normalized by range)
results.mape = 100*mean(abs(errs)./abs(y_actual),axis,'omitnan');  % mean absolute percent error
results.medape = 100*median(abs(errs)./abs(y_actual),axis,'omitnan'); % median abs percent error
results.mane = mean(abs(errs)./max(abs(y_actual),[],axis,'omitnan'),axis,'omitnan'); % mean absolute normalized error (normalized by max abs of validation set)
results.mane2 = mean(abs(errs)./quantile(abs(y_actual),0.975,axis),axis,'omitnan'); % mean absolute normalized error (normalized by 0.975 quantile of validation set)
results.medane = median(abs(errs)./max(abs(y_actual),[],axis,'omitnan'),axis,'omitnan'); % median absolute normalized error (normalized by max abs of validation set)
results.medane2 = median(abs(errs)./quantile(abs(y_actual),0.975,axis),axis,'omitnan'); % median absolute normalized error (normalized by 0.975 quantile of validation set)
rsquared = zeros(1,size(y_actual,2)); 
p_value = zeros(1,size(y_actual,2)); 
for i = 1:size(y_actual,2)
    [~,~,~,~,stats] = regress(y_actual(:,i),[y_pred(:,i), ones(size(y_pred,1),1)]);
    rsquared(i) = stats(1); 
    p_value(i) = stats(3); 
end
results.rsquared = rsquared; 
results.p_value = p_value; 
end