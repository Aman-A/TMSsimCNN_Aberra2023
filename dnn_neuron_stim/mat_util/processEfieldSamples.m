
function [Efield_samples_out, varargout] = processEfieldSamples(Efield_samples_in,...
                                                                process_mode,N,...
                                                                varargin)
% Takes in cell array of E-field samples at each position, outputs feature
% vector at each position consisting of E-field components in either
% cartesian, spherical coordinates and either direct or transformed
% (described below)
%**Inputs:
% Efield_samples_in : cell array of size (Number of positions,1)
%                  Each entry must contain a (Number of Samples,3) array
%                  containing E-field samples centered at soma and rotated
%                  into local coordinate space
% process_mode : string, one of 'cart','sph','cart_tr','sph_tr'
% process_mode == 'cart' --> Cartesian
% process_mode == 'sph' --> Spherical.
% process_mode == 'cart_tr' -->  Cartesian Transformed
% process_mode == 'sph_tr'--> Spherical Transformed.
% process_mode == '3D_cart'
% process_mode == '3D_sph'
% process_mode == '3D_sphPD' - 3D spherical with phi (azimuthal angle) differences
% process_mode == '3D_sphNorm' - 3D Spherical w normalized polar and
% azimuthal angles
%**Outputs:
% Efield_samples_out : cell array of size (Number of positions,1)
%                      Contains transformed Efield vectors
%% Processing steps
%note each Efield sample should be of same size
%Keep Cartesian
in.drop_phi = 0; % remove azimuthal angle for 3D_sph mode, and stat_method ~= 'cell'
in.theta_mode = 'norm90'; % 'norm' or 'norm180' to normalize to [0,1], 'norm90'
                        % to normalized to [-1 1], 'cosine' to take cos(theta)
in.phi_mode = 'norm180';  % 'cosine'/'sine' (same thing), or 'norm180'
in.reorder_E = 0; % previously required reordering Efield samples to grid order, set to 1 for old datasets
in.print_level = 1;
in = sl.in.processVarargin(in,varargin);
if in.print_level > 0
    fprintf('Processing E-field samples with mode: %s\n',process_mode);
    fprintf('WARNING: If using old Ecell data, check row order. New data (after 3/12/21) uses default column order\n'); 
end
Efield_samples_out = cell(size(Efield_samples_in));
switch process_mode
    case '3D_cart'
        for i = 1:length(Efield_samples_in)
            Ei = Efield_samples_in{i};
            if in.reorder_E
                % since we moved [0 0 0] to the first row, re-order to move it
                % back to middle row
                nrows = uint64(size(Ei,1)); % make int to use in index below (rounds up)
                % Reorder to same output as meshgrid
                Ei = [Ei(2:nrows/2,:);Ei(1,:);Ei((nrows/2 + 1):end,:)];
            end
            % Reshape
            EX = reshape(Ei(:,1),N,N,N);
            EY = reshape(Ei(:,2),N,N,N);
            EZ = reshape(Ei(:,3),N,N,N);
            % Concatenate
            E4d = cat(4,EX,EY,EZ); % make the final NxNxNx3 tensor
            Efield_samples_out{i} = E4d;
        end
    case '3D_sph'
        for i = 1:length(Efield_samples_in)
            Ei = Efield_samples_in{i};
            if in.reorder_E
                % since we moved [0 0 0] to the first row, re-order to move it
                % back to middle row
                nrows = uint64(size(Ei,1)); % make int to use in index below (rounds up)
                % Reorder to same output as meshgrid
                Ei = [Ei(2:nrows/2,:);Ei(1,:);Ei((nrows/2 + 1):end,:)];
            end
            [thetas,phis,mags] = cart2sphD(Ei(:,1),Ei(:,2),Ei(:,3),'phi_mode',1); % outputs angles 0-360
            thetaTransform = norm_angles(thetas,in.theta_mode);

            % Reshape
            MAGS = reshape(mags,N,N,N); % Magnitudes
            THETA = reshape(thetaTransform,N,N,N); % cos of polar angles (theta)
            if in.drop_phi == 0
                if strcmp(in.phi_mode,'cos') || strcmp(in.phi_mode,'sin')
                    % Include azimuthal angles
                    phiTransform1 = cosd(phis); % include azimuthal angle of each E-field vector
                    phiTransform2 = sind(phis);
                    COSPHI = reshape(phiTransform1,N,N,N); % cos of azimuthal angles (phi)
                    SINPHI = reshape(phiTransform2,N,N,N); % sin of azimuthal angles (phi)
                    % Concatenate
                    E4d = cat(4,MAGS,THETA,COSPHI,SINPHI); % make the final NxNxNx4 tensor
                else
                    phiNORM = norm_angles(phis,in.phi_mode);
                    PHI = reshape(phiNORM,N,N,N);
                    E4d = cat(4,MAGS,THETA,PHI);
                end
            elseif in.drop_phi == 1
                % Concatenate
                E4d = cat(4,MAGS,THETA); % make the final NxNxNx2 tensor
            end
            Efield_samples_out{i} = E4d;
        end
    case '3D_sphPD'
        phi_cinds = zeros(length(Efield_samples_in),1);
        for i = 1:length(Efield_samples_in)
            Ei = Efield_samples_in{i};
            if in.reorder_E
                % since we moved [0 0 0] to the first row, re-order to move it
                % back to middle row
                nrows = uint64(size(Ei,1)); % make int to use in index below (rounds up)
                % Reorder to same output as meshgrid
                Ei = [Ei(2:nrows/2,:);Ei(1,:);Ei((nrows/2 + 1):end,:)];
            end
            [thetas,phis,mags] = cart2sphD(Ei(:,1),Ei(:,2),Ei(:,3),'phi_mode',1); % outputs angles 0-360
            thetaTransform = norm_angles(thetas,in.theta_mode);
            phi_cinds(i) = phis(cind); % Save phi at center for testing unif E approximation on this dataset
            phi_diffs = phis - phis(cind); % phi is 0-360
            phi_diffs(phi_diffs >= 180) = phi_diffs(phi_diffs >= 180) - 360; % [-180,180]
            phi_diffs(phi_diffs <= -180) = phi_diffs(phi_diffs <= -180) + 360;
            if any(abs(phi_diffs) >= 180)
               fprintf('WARNING: %g phis > 180 deg from center\n',sum(abs(phi_diffs) >= 180));
            end
            % Reshape to 4D tensor
            MAGS = reshape(mags,N,N,N); % Magnitudes
            THETA = reshape(thetaTransform,N,N,N); % cos of polar angles (theta)
            if strcmp(in.phi_mode,'norm180') || strcmp(in.phi_mode,'norm')
                phi_diffs = phi_diffs / 180; % normalize azimuthal differences to [-1 1]
                NORMPHIDIFF = reshape(phi_diffs,N,N,N); % cos of azimuthal angles (phi)
                % Concatenate
                E4d = cat(4,MAGS,THETA,NORMPHIDIFF); % make the final NxNxNx3 tensor
            elseif strncmp('cos',in.phi_mode,3) || strcmp('sin',in.phi_mode,3)
                phi_diffsTransform1 = cosd(phi_diffs); % cos of azimuthal angle differences (vs. center)
                phi_diffsTransform2 = sind(phi_diffs); % sin of azimuthal angle differences (vs. center)
                % Reshape to 4D tensor
                COSPHIDIFFS = reshape(phi_diffsTransform1,N,N,N);
                SINPHIDIFFS = reshape(phi_diffsTransform2,N,N,N);
                % Concatenate
                E4d = cat(4,MAGS,THETA,COSPHIDIFFS,SINPHIDIFFS); % make the final NxNxNx3 tensor
            end
            Efield_samples_out{i} = E4d;
        end
    case 'cart'
        % Keep cartesian, convert to linear vectors
        for i = 1:length(Efield_samples_in)
            % transpose so output is linear vector with components of same
            % vector are kept together [Ex_1;Ey_1;Ez_1;Ex_2;Ey_2;Ez_2;...]
            Ei = Efield_samples_in{i}';
            Efield_samples_out{i} = Ei(:);
        end
    case 'sph'
        % Convert to spherical coords and apply sin/cosine transforms
        for i = 1:length(Efield_samples_in)
            [thetas,phis,mags] = cart2sphD(Efield_samples_in{i}(:,1),Efield_samples_in{i}(:,2),Efield_samples_in{i}(:,3));
            thetaTransform = cosd(thetas);
            phiTransform1 = cosd(phis);
            phiTransform2 = sind(phis);
            processedEi = [mags,thetaTransform,phiTransform1,phiTransform2]'; % take transpose
            % transpose so output is linear vector:
            % [mag_1;thetaTransform_1;phiTransform1_1;phiTransform2_1;mag_2;thetaTransform_2;phiTransform1_2;phiTransform2_2;...]
            Efield_samples_out{i} = processedEi(:);
        end
    case 'cart_tr'
        % Keep cartesian, transform by taking difference vectors with center
        %extract Ec: each Efield sample should be same size, thus taking median
        %of 1st sample should give index of center Efield vector of all samples
        % *****
        % AA: check interpEfieldSample.m and samplePts.m to see how E-field samples
        % are ordered for each sample_method_struct.method, for 'box', currently,
        % center E-field is first row
        % *****
    %     Ec_ind = ceil(length(Efield_samples_in{1})./2);
        Nvecs = size(Efield_samples_in{1},1); % assume same for all positions
        Ec_ind = 1;
        off_c_inds = setdiff(1:Nvecs,Ec_ind);
        for i = 1:length(Efield_samples_in)
            Ec = Efield_samples_in{i}(Ec_ind,:); % center E-field vector
            Ec_mag = norm(Ec);
            unit_Ec = Ec/Ec_mag;
            Ei = Efield_samples_in{i}(off_c_inds,:); % off center E-field vectors
            diff_Ei = (Ei - Ec)/Ec_mag; % get normalized difference vectors
            % transpose so output is linear vector:
            % [mag_c;Ex_c;Ey_c;Ez_c;diffEx_1;diffEy_1;diffEz_1;diffEx_2;diffEy_2;diffEz_2;...]
            diff_Ei = diff_Ei'; % transpose so outputs are
            processedEi = [Ec_mag;
                                unit_Ec(:);
                                diff_Ei(:)];
            Efield_samples_out{i} = processedEi;
        end
    case 'sph_tr'
        % Convert to spherical coords and apply sin/cosine transforms,
        % transform by taking difference vectors with center
        Nvecs = size(Efield_samples_in{1},1); % assume same for all positions
        Ec_ind = 1;
        off_c_inds = setdiff(1:Nvecs,Ec_ind);
        for i = 1:length(Efield_samples_in)
            Ec = Efield_samples_in{i}(Ec_ind,:); % center E-field vector
            [Ec_th,Ec_ph,Ec_mag] = cart2sphD(Ec(1),Ec(2),Ec(3));
            % this is almost same as converting back to cartesian. Formulas for
            % x and y in cartesian are: x = sin(theta)*cos(phi);
            % y = sin(theta)*sin(phi)
            % so might not do much, but we can test it out...
            Ec_dir = [cosd(Ec_th);cosd(Ec_ph);sind(Ec_ph)];
            Ei = Efield_samples_in{i}(off_c_inds,:); % off center E-field vectors
            diff_Ei = (Ei - Ec)/Ec_mag; % get normalized difference vectors
            % convert norm. diff vectors to spherical coords
            [diff_Ei_th,diff_Ei_ph,diff_Ei_mag] = cart2sphD(diff_Ei(:,1),...
                                                            diff_Ei(:,2),...
                                                            diff_Ei(:,3));
            % apply sin/cos transformations and transpose so output is
            % linear vector with components from same vector adjacent elements:
            diff_Ei_sph = [diff_Ei_mag,... % magnitude of norm. diff vectors
                           cosd(diff_Ei_th),...
                           cosd(diff_Ei_ph),...
                           sind(diff_Ei_ph)];
            % [mag_c;cos_th_c;cos_ph_c;sin_ph_c;diffE_mag_1;diffE_cos_th_1;diffE_cos_ph_1;diffE_sin_ph_1;...]
            diff_Ei_sph = diff_Ei_sph'; % transpose so outputs are
            processedEi = [Ec_mag;
                                Ec_dir;
                                diff_Ei_sph(:)];
            Efield_samples_out{i} = processedEi;
        end
end
if nargout == 2
    if strcmp(process_mode,'3D_sphPD')
        varargout = {phi_cinds}; % also output azimuthal angle of center point
                                 % if using differences for testing uniform
                                 % E-field approximation on same dataset
    else
        varargout = {[]}; % empty output
    end
end
end
function angleTransform = norm_angles(angles,norm_mode)
if any(strcmp(norm_mode,{'norm','norm180'}))
    if (max(abs(angles)) <= 180)
        angleTransform = angles/180; % normalize polar angle
    else
        angleTransform = (180 - angles)/180;
    end
elseif strcmp(norm_mode,'norm90') % normalize to [0 1] (0 up, 1 down)
    angleTransform = (90 - angles)/90; % normalized to [-1 1] (1 up, -1 down)
elseif any(strcmp(norm_mode,{'cosine','cos'}))
    angleTransform = cosd(angles);
elseif any(strcmp(norm_mode,{'norm360'}))
    angleTransform = angles/360;
elseif any(strcmp(norm_mode,{'norm360_180'}))
            angleTransform = (180 - angles)/180;
end
end