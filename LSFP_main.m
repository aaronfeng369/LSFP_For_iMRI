
% proposed algorithm LSFP for simulation
% March-29-2021==Zhao He==Original code
% April-06-2021==Zhao He==Revise and add comments 
clc; clear; close all;
addpath(genpath('./functions'));

%% load data
load('IMRI_Simu2_400.mat'); 
Nsample = size(kdata,1); % number of sampling in readout direction
overSamp = 2; % oversampling in readout direction
mC = 0; % parameter for cropping data

% parameters
Ntviews = 2000; % total number of spokes
Nspokes = 10; % spokes per frame
Ng = 5; % frames per group

tp = 30; % specific group for quantitative evaluation
RA = 'LSFP'; % reconstruction algorithm

%% Reconstruction

Nt = Ntviews/Nspokes;% total frames
Mg = Nt/Ng; % total groups
indG = reshape(1:Nt, [Ng, Mg])'; % spokes grouping

% kdata grouping
for ii = 1:Nt
    kdatau(:,:,:,ii) = kdata(:,(ii-1)*Nspokes+1:ii*Nspokes,:);
    ku(:,:,ii) = k(:,(ii-1)*Nspokes+1:ii*Nspokes);
end
wu = abs(ku);

% reconstruct specific group
position = tp; 
kdata1 = kdatau(:,:,:,indG(position,:));
ku1 = ku(:,:,indG(position,:));
wu1 = wu(:,:,indG(position,:));
param.E = MCNUFFT(ku1,wu1,b1); 
param.d = kdata1;
recon_nufft = param.E'*param.d; 

param.lambda_L = 0.025;
param.lambda_S = 0.5*max(abs(recon_nufft(:)));
param.lambda_spatial_L = 5e-5;
param.lambda_spatial_S = 5e-5;
param.nite = 40;
param.tol = 0.0025;

% LSFP reconstruction
tic;
[L1,S1] = lsfp(param);
toc;
LS  = flipdim(L1+S1,1);

% crop data 
LplusS_crop = flipud(crop(LS,Nsample/overSamp-mC,Nsample/overSamp-mC,Ng));

%% save data 
savename = [RA,'nsam',num2str(Nsample),'_nspokes',num2str(Nspokes),'_ng',num2str(Ng),...
            '_nt',num2str(Ng),'_f',num2str(position*Ng),'_',RA,'_it',num2str(param.nite)];
save([savename,'.mat'],'LS');
savegif([savename,'_LpS.gif'],LplusS_crop,10);






