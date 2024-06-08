clear;
clc;

HOME = pwd;
addpath([HOME '/input']);
addpath([HOME '/Tools']);

visual_gmtfile("crust1.bd2.gmt");
visual_gmtfile("btfield.gmt");

% Load Bouguer field
bcfield_data = load("bcfield.gmt");
[bcfield, Lon, Lat] = gmt2matrix(bcfield_data);
bcfield = Europe_centered(bcfield);
lon = Lon(1, :);
lat = Lat(:, 1);
lon = lon - 180;

figure;
hold on
imagesc(lon, lat, bcfield);
hold off