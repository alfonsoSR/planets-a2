clear
clc
HOME = pwd;
addpath([HOME '/input']);
addpath([HOME '/Tools']);
addpath([HOME '/input/custom']);

% Model setup
Model = struct();
Model.name = "BouguerInversion";
Model.GM = 398600441500000.0;
Model.Re = 6378137.0;
Model.Re_analyse = Model.Re;
Model.geoid = 'none';
Model.nmax = 90;
Model.number_of_layers = 2;

% Define boundaries
D = -35e3;
[upper_bound, lon, lat] = gmt2matrix(load("topography.gmt"));
middle_bound = ones(size(upper_bound)) * D;

% Define latitude and longitude ranges
lon = lon(1, :);
lat = lat(:, 1);
dlon = lon(2) - lon(1);
dlat = lat(1) - lat(2);
lon_range = [min(lon), max(lon), dlon];
lat_range = [min(lat), max(lat) + dlat, dlat];

% Generate reference gravity field
[obs_gravity, ~, ~] = gmt2matrix(load("bouguer_anomalies.gmt"));
% obs_coeffs = load("gravity_coeffs.gmt");
% obs_coeffs(1, 3) = 0.;
% obs_coeffs(3, 3) = 0.;
% [obs_data] = model_SH_synthesis(...
%     lon_range, lat_range, 0., [0, Model.nmax], obs_coeffs, Model);
% obs_gravity = flip(sqrt(...
%     obs_data.vec.X .* obs_data.vec.X + ...
%     obs_data.vec.Y .* obs_data.vec.Y + ...
%     obs_data.vec.Z .* obs_data.vec.Z ...
% ));

% Generate boundary
alpha = 0.25;
while true

    % Compute coefficients with current boundary
    [gcoeffs] = segment_2layer_model(...
        upper_bound, middle_bound, -200e3, 2670., 3300., 25e3, Model);
    gcoeffs(1, 3) = 0.;
    gcoeffs(3, 3) = 0.;

    % Generate gravity field from computed coefficients
    [gfield_data] = model_SH_synthesis(...
        lon_range, lat_range, 0., [0, Model.nmax], gcoeffs, Model);
    computed_gravity = flip(sqrt(...
        gfield_data.vec.X .* gfield_data.vec.X + ...
        gfield_data.vec.Y .* gfield_data.vec.Y + ...
        gfield_data.vec.Z .* gfield_data.vec.Z ... 
    ));

    % Compute residual
    res = (obs_gravity - computed_gravity) * 1e5;
    val = max(max(abs(res)));
    disp(val);
    if val < 100
        break
    else
        disp(val);
    end

    % Update interface
    middle_bound = middle_bound + alpha * res;

end
