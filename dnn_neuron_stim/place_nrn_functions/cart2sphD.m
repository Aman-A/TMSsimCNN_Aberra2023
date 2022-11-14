function [thetas,phis,mags] = cart2sphD(x,y,z,varargin)
%CART2SPHD Transform Cartesian to spherical coordinates in degrees
%  
%   Inputs 
%   ------ 
%   Optional Inputs 
%   --------------- 
%   Outputs 
%   ------- 
%   thetas : vector 
%            polar angles in degrees (0 to +180 deg) from positive z axis
%   phis : vector
%          azimuthal angle in degrees (0 to +360 deg) from positive x axis
%   mags : vector
%          radius from origin (vector magnitude)
%   Examples 
%   --------------- 

% AUTHOR    : Aman Aberra 
in.phi_mode = 1; % 0 - phi is same as output of cart2sph (in deg)
                 % 1 - phi is 0 to 360 deg (all positive)
                 % 2 - phi is -180 to +180 deg
in = sl.in.processVarargin(in,varargin);
[phis_rad,lambdas_rad,mags] = cart2sph(x,y,z); % outputs azimuth and elevation of vectors in rad, and magnitude
thetas = 90 - lambdas_rad*180/pi; % convert to deg and shift range from [-90,+90] to [0,180]
if in.phi_mode == 0
    phis = phis_rad*180/pi; % convert to deg
elseif in.phi_mode == 1
    phis = phis_rad*180/pi; % convert to deg
    phis(phis < 0) = 360 + phis(phis < 0); 
elseif in.phi_mode == 2
    phis = phis_rad*180/pi; % convert to deg
    phis(phis >= 180) = phis(phis >= 180) - 360; 
    phis(phis <= -180) = phis(phis <= -180) + 360; 
else
    error('phi_mode must be 0, 1, or 2');
end
