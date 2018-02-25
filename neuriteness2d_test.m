%% clear
clc; clear all; close all;

%% path
addpath('./lib')
addpath('../vesselness2d/lib')

%% load image
im = imread ('../vesselness2d/im/jellyfish.png');

%% vesselness
sigma = 3; 
gamma = 2; 

[imn,vx,vy] = neuriteness2d(im,sigma,gamma);

%% plot
figure; imagesc(im); colormap gray; 
set(gca,'ytick',[]); set(gca,'xtick',[]); axis image; axis tight;

figure; imagesc(imn); colormap gray; 
set(gca,'ytick',[]); set(gca,'xtick',[]); axis image; axis tight;