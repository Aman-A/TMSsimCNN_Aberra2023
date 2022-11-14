%% Plots all maps parameters should already be set
save_figs = 0;
cbar_mode = 1; % 1 horz 2 vert
save_data = 1; % saves generated data for faster replotting
% CNN
plotCNNThreshLayer_ernie_test_L2_Fig(save_figs,cbar_mode,save_data)
plotCNNThreshLayer_ernie_test_L3_Fig(save_figs,cbar_mode,save_data)
plotCNNThreshLayer_ernie_test_L4_Fig(save_figs,cbar_mode,save_data)
% Unif E
plotUnifEMapThreshLayer_ernie_test_L2_Fig(save_figs,cbar_mode,save_data)
plotUnifEMapThreshLayer_ernie_test_L3_Fig(save_figs,cbar_mode,save_data)
plotUnifEMapThreshLayer_ernie_test_L4_Fig(save_figs,cbar_mode,save_data)