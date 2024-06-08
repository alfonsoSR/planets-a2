% plot results
clear;
close all;
clc;

%% tutorial data
%%% insert output data file from Results here!!!%%%
source_dir = "./output/";
load(source_dir + "BouguerInversion_field.mat");
root = "./output/struct/";
% 
% load("../../output/data_Crust10_crust_0_179_03-Jun-2024 10:08:13.mat")
% % load('Results/data_Crust10_crust_3_179_26-Mar-2019 13:14:45.mat')
% root = "../../output/struct/";

%% plot different maps of the data
lon = data.grd.lon(1,:);
lats = data.grd.lat(:,1);

% pot = data.pot;
% lat = data.grd.lat;
% lon = data.grd.lon;
% radius = data.grd.r;
% R = data.vec.R;
% T = data.vec.T;
% L = data.vec.L;
% X = data.vec.X;
% Y = data.vec.Y;
% Z = data.vec.Z;
% Trr = data.ten.Trr;
% Ttt = data.ten.Ttt;
% Tll = data.ten.Tll;
% Trt = data.ten.Trt;
% Trl = data.ten.Trl;
% Ttl = data.ten.Ttl;
% Txx = data.ten.Txx;
% Tyy = data.ten.Tyy;
% Tzz = data.ten.Tzz;
% Txy = data.ten.Txy;
% Txz = data.ten.Txz;
% Tyz = data.ten.Tyz;
% latLim = data.latLim;
% lonLim = data.lonLim;
% height = data.height;
% save(strcat(root, "pot.txt"), "pot", "-ascii");
% save(strcat(root, "lat.txt"), "lat", "-ascii");
% save(strcat(root, "lon.txt"), "lon", "-ascii");
% save(strcat(root, "radius.txt"), "radius", "-ascii");
% save(strcat(root, "R.txt"), "R", "-ascii");
% save(strcat(root, "T.txt"), "T", "-ascii");
% save(strcat(root, "L.txt"), "L", "-ascii");
% save(strcat(root, "X.txt"), "X", "-ascii");
% save(strcat(root, "Y.txt"), "Y", "-ascii");
% save(strcat(root, "Z.txt"), "Z", "-ascii");
% save(strcat(root, "Trr.txt"), "Trr", "-ascii");
% save(strcat(root, "Ttt.txt"), "Ttt", "-ascii");
% save(strcat(root, "Tll.txt"), "Tll", "-ascii");
% save(strcat(root, "Trt.txt"), "Trt", "-ascii");
% save(strcat(root, "Trl.txt"), "Trl", "-ascii");
% save(strcat(root, "Ttl.txt"), "Ttl", "-ascii");
% save(strcat(root, "Txx.txt"), "Txx", "-ascii");
% save(strcat(root, "Tyy.txt"), "Tyy", "-ascii");
% save(strcat(root, "Tzz.txt"), "Tzz", "-ascii");
% save(strcat(root, "Txy.txt"), "Txy", "-ascii");
% save(strcat(root, "Txz.txt"), "Txz", "-ascii");
% save(strcat(root, "Tyz.txt"), "Tyz", "-ascii");
% save(strcat(root, "latLim.txt"), "latLim", "-ascii");
% save(strcat(root, "lonLim.txt"), "lonLim", "-ascii");
% save(strcat(root, "height.txt"), "height", "-ascii");

figure;
subplot(2,2,1)
imagesc(lon,lats,((data.pot)));c=colorbar;
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Potential gravity field'])
ylabel(c,'m*m/s/s')
set(gca,'YDir','normal')

subplot(2,2,2)
imagesc(lon,lats,((data.vec.Z)).*1e5);c=colorbar;
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Z-component of gravity vector'])
ylabel(c,'mGal')
set(gca,'YDir','normal')

subplot(2,2,3)
imagesc(lon,lats,((data.vec.X)).*1e5);c=colorbar;
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['X-component of gravity vector (North-South)'])
ylabel(c,'mGal')
set(gca,'YDir','normal')

subplot(2,2,4)
imagesc(lon,lats,((data.vec.Y)).*1e5);c=colorbar;
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Y-component of gravity vector (East-West)'])
ylabel(c,'mGal')
set(gca,'YDir','normal')

%% Tensor

figure;
subplot(3,3,1)
imagesc(lon,lats,((data.ten.Tzz).*1e9));c=colorbar;
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Tzz-component of gravity gradient tensor'])
ylabel(c,'Eotvos')
set(gca,'YDir','normal')

subplot(3,3,2)
imagesc(lon,lats,((data.ten.Txz).*1e9));c=colorbar;
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Txz-component of gravity gradient tensor'])
ylabel(c,'Eotvos')
set(gca,'YDir','normal')

subplot(3,3,3)
imagesc(lon,lats,((data.ten.Tyz).*1e9));c=colorbar;
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Tyz-component of gravity gradient tensor'])
ylabel(c,'Eotvos')
set(gca,'YDir','normal')

subplot(3,3,5)
imagesc(lon,lats,((data.ten.Txx).*1e9));c=colorbar;
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Txx-component of gravity gradient tensor'])
ylabel(c,'Eotvos')
set(gca,'YDir','normal')

subplot(3,3,6)
imagesc(lon,lats,((data.ten.Txy).*1e9));c=colorbar;
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Txy-component of gravity gradient tensor'])
ylabel(c,'Eotvos')
set(gca,'YDir','normal')

subplot(3,3,9)
imagesc(lon,lats,((data.ten.Tyy).*1e9));c=colorbar;
hold on
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Tyy-component of gravity gradient tensor'])
ylabel(c,'Eotvos')
set(gca,'YDir','normal')
