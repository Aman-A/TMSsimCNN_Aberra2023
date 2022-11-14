function Eopts = Eopts_handler(Eopts)
%EOPTS_HANDLER ... 
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
Eopts_default.drop_phi = 0; % drop azimuth for 3D_sph mode in processEfieldSamples, default keep
Eopts_default.theta_mode = 'norm180'; % normalize by 180 [0,1]
Eopts_default.phi_mode = 'norm180'; % normalize by 180 [-1,1]
Eopts_default.print_level = 0; 
Eopts_default.reorder_E = 0; % For older Ecells in which 1st row is center E
Eopts = sl.in.processVarargin(Eopts_default,Eopts); 
valid_theta_modes = {'norm180','norm','norm90','cos','cosine'}; 
valid_phi_modes = {'cos','cosine','sin','sine','norm','norm360','norm360_180'};
assert(any(strcmp(Eopts.theta_mode,valid_theta_modes)),...
    sprintf('theta_mode %s not defined',Eopts.theta_mode))
assert(any(strcmp(Eopts.phi_mode,valid_phi_modes)),...
    sprintf('phi_mode %s not defined',Eopts.phi_mode))