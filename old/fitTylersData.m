
% Fit Tyler's speed prior
tic
clear all
close all

%% Load in within screen datasets

% 0.5deg/s (ref: left column, test: right column)
load('/home/tyler/Documents/MATLAB/cooperLab/3-Data/priorCoordsRDS/TSM/TSM_RDS_05dps/TSM-20220623T160924-1.mat','pa','ds');

screenDists = ds.screenDistance;
within = pa.refScreen == pa.testScreen;

% Ind 1: Ref, ind 2: test
stimV1 = [pa.response(:,5) pa.response(:,3)];
stimC1 = [pa.response(:,4) pa.response(:,2)];
if within
    stimD1 = [screenDists(pa.response(:,6)) screenDists(pa.response(:,6))];
else
    stimD1 = [screenDists(~(pa.response(:,6)-1)+1) screenDists(pa.response(:,6))];
end

% convert to ref chosen y/n 
r1     = pa.response(:,9) == pa.response(:,8);

% Completion time
t1 = pa.currentTime;

clear pa screenDists within

% 1deg/s
load('/home/tyler/Documents/MATLAB/cooperLab/3-Data/priorCoordsRDS/TSM/TSM_RDS_1dps/TSM-20220623T164058-1.mat','pa','ds');

screenDists = ds.screenDistance;
within = pa.refScreen == pa.testScreen;

stimV2 = [pa.response(:,5) pa.response(:,3)];
stimC2 = [pa.response(:,4) pa.response(:,2)];
if within
    stimD2 = [screenDists(pa.response(:,6)) screenDists(pa.response(:,6))];
else
    stimD2 = [screenDists(~(pa.response(:,6)-1)+1) screenDists(pa.response(:,6))];
end
r2     = pa.response(:,9) == pa.response(:,8);

% Completion time
t2 = pa.currentTime;

clear pa

% 2deg/s
load('/home/tyler/Documents/MATLAB/cooperLab/3-Data/priorCoordsRDS/TSM/TSM_RDS_2dps/TSM-20220623T170410-1.mat','pa','ds');

screenDists = ds.screenDistance;
within = pa.refScreen == pa.testScreen;

stimV3 = [pa.response(:,5) pa.response(:,3)];
stimC3 = [pa.response(:,4) pa.response(:,2)];
if within
    stimD3 = [screenDists(pa.response(:,6)) screenDists(pa.response(:,6))];
else
    stimD3 = [screenDists(~(pa.response(:,6)-1)+1) screenDists(pa.response(:,6))];
end
r3     = pa.response(:,9) == pa.response(:,8);

% Completion time
t3 = pa.currentTime;

clear pa screenDists within

% 4deg/s
load('/home/tyler/Documents/MATLAB/cooperLab/3-Data/priorCoordsRDS/TSM/TSM_RDS_4dps/TSM-20220623T172007-1.mat','pa','ds');

screenDists = ds.screenDistance;
within = pa.refScreen == pa.testScreen;

stimV4 = [pa.response(:,5) pa.response(:,3)];
stimC4 = [pa.response(:,4) pa.response(:,2)];
if within
    stimD4 = [screenDists(pa.response(:,6)) screenDists(pa.response(:,6))];
else
    stimD4 = [screenDists(~(pa.response(:,6)-1)+1) screenDists(pa.response(:,6))];
end
r4     = pa.response(:,9) == pa.response(:,8);

% Completion time
t4 = pa.currentTime;

clear pa screenDists within

% 8deg/s
load('/home/tyler/Documents/MATLAB/cooperLab/3-Data/priorCoordsRDS/TSM/TSM_RDS_8dps/TSM-20220623T173416-1.mat','pa','ds');

screenDists = ds.screenDistance;
within = pa.refScreen == pa.testScreen;

stimV5 = [pa.response(:,5) pa.response(:,3)];
stimC5 = [pa.response(:,4) pa.response(:,2)];
if within
    stimD5 = [screenDists(pa.response(:,6)) screenDists(pa.response(:,6))];
else
    stimD5 = [screenDists(~(pa.response(:,6)-1)+1) screenDists(pa.response(:,6))];
end
r5     = pa.response(:,9) == pa.response(:,8);

% Completion time
t5 = pa.currentTime;

clear pa screenDists within

%% Load in between screen datasets

% 0.25m
load('/home/tyler/Documents/MATLAB/cooperLab/3-Data/priorCoordsRDS/TSM/TSM_gabor_0_25m/TSM-20220614T174612-2.mat','pa','ds');

screenDists = ds.screenDistance;
within = pa.refScreen == pa.testScreen;

stimV6 = [pa.response(:,5) pa.response(:,3)];
stimC6 = [pa.response(:,4) pa.response(:,2)];
if within
    stimD6 = [screenDists(pa.response(:,6)) screenDists(pa.response(:,6))];
else
    stimD6 = [screenDists(~(pa.response(:,6)-1)+1) screenDists(pa.response(:,6))];
end

% Distance manipulation via size?
if isfield(pa,'perspective')
    if pa.perspective == 0
        farDist = min(screenDists);
    elseif (pa.perspective ~= 1) && (pa.perspective ~= 0)
        farDist = 1/pa.perspective;
    else
        farDist = 1;
    end
end
stimD6(stimD6==max(screenDists)) = farDist;

r6     = pa.response(:,9) == (~(pa.response(:,6)-1)+1);

% Completion time
t6 = pa.currentTime;

clear pa screenDists within farDist

% 0.5m
load('/home/tyler/Documents/MATLAB/cooperLab/3-Data/priorCoordsRDS/TSM/TSM_gabor_retinalSize/TSM-20220613T173303-2.mat','pa','ds');

screenDists = ds.screenDistance;
within = pa.refScreen == pa.testScreen;

stimV7 = [pa.response(:,5) pa.response(:,3)];
stimC7 = [pa.response(:,4) pa.response(:,2)];
if within
    stimD7 = [screenDists(pa.response(:,6)) screenDists(pa.response(:,6))];
else
    stimD7 = [screenDists(~(pa.response(:,6)-1)+1) screenDists(pa.response(:,6))];
end

% Distance manipulation via size?
if isfield(pa,'perspective')
    if pa.perspective == 0
        farDist = min(screenDists);
    elseif (pa.perspective ~= 1) && (pa.perspective ~= 0)
        farDist = 1/pa.perspective;
    else
        farDist = 1;
    end
else
    farDist = min(screenDists);
end
stimD7(stimD7==max(screenDists)) = farDist;

r7     = pa.response(:,9) == (~(pa.response(:,6)-1)+1);

% Completion time
t7 = pa.currentTime;

clear pa screenDists within farDist

% % 0.75m
% load('/home/tyler/Documents/MATLAB/cooperLab/3-Data/priorCoordsRDS/TSM/TSM_gabor_0_75m/TSM-20220614T140927-2.mat','pa','ds');
% 
% screenDists = ds.screenDistance;
% within = pa.refScreen == pa.testScreen;
% 
% stimV8 = [pa.response(:,5) pa.response(:,3)];
% stimC8 = [pa.response(:,4) pa.response(:,2)];
% if within
%     stimD8 = [screenDists(pa.response(:,6)) screenDists(pa.response(:,6))];
% else
%     stimD8 = [screenDists(~(pa.response(:,6)-1)+1) screenDists(pa.response(:,6))];
% end
% 
% % Distance manipulation via size?
% if isfield(pa,'perspective')
%     if pa.perspective == 0
%         farDist = min(screenDists);
%     elseif (pa.perspective ~= 1) && (pa.perspective ~= 0)
%         farDist = 1/pa.perspective;
%     else
%         farDist = 1;
%     end
% end
% stimD8(stimD8==max(screenDists)) = farDist;
% 
% r8     = pa.response(:,9) == (~(pa.response(:,6)-1)+1);
% 
% % Completion time
% t8 = pa.currentTime;
% 
% clear pa screenDists within farDist

% 1.25m
load('/home/tyler/Documents/MATLAB/cooperLab/3-Data/priorCoordsRDS/TSM/TSM_gabor_1_25m/TSM-20220614T162824-2.mat','pa','ds');

screenDists = ds.screenDistance;
within = pa.refScreen == pa.testScreen;

stimV9 = [pa.response(:,5) pa.response(:,3)];
stimC9 = [pa.response(:,4) pa.response(:,2)];
if within
    stimD9 = [screenDists(pa.response(:,6)) screenDists(pa.response(:,6))];
else
    stimD9 = [screenDists(~(pa.response(:,6)-1)+1) screenDists(pa.response(:,6))];
end

% Distance manipulation via size?
if isfield(pa,'perspective')
    if pa.perspective == 0
        farDist = min(screenDists);
    elseif (pa.perspective ~= 1) && (pa.perspective ~= 0)
        farDist = 1/pa.perspective;
    else
        farDist = 1;
    end
end
stimD9(stimD9==max(screenDists)) = farDist;

r9     = pa.response(:,9) == (~(pa.response(:,6)-1)+1);

% Completion time
t9 = pa.currentTime;

clear pa screenDists within farDist

% 1.0m
load('/home/tyler/Documents/MATLAB/cooperLab/3-Data/priorCoordsRDS/TSM/TSM_gabor_worldSize/TSM-20220613T164036-2.mat','pa','ds');

screenDists = ds.screenDistance;
within = pa.refScreen == pa.testScreen;

stimV10 = [pa.response(:,5) pa.response(:,3)];
stimC10 = [pa.response(:,4) pa.response(:,2)];
if within
    stimD10 = [screenDists(pa.response(:,6)) screenDists(pa.response(:,6))];
else
    stimD10 = [screenDists(~(pa.response(:,6)-1)+1) screenDists(pa.response(:,6))];
end

% Distance manipulation via size?
if isfield(pa,'perspective')
    if pa.perspective == 0
        farDist = min(screenDists);
    elseif (pa.perspective ~= 1) && (pa.perspective ~= 0)
        farDist = 1/pa.perspective;
    else
        farDist = 1;
    end
else
    farDist = max(screenDists);
end
stimD10(stimD10==max(screenDists)) = farDist;

r10     = pa.response(:,9) == (~(pa.response(:,6)-1)+1);

% Completion time
t10 = pa.currentTime;

clear pa screenDists within farDist


%% Combine data
% stimV = [stimV1;stimV2;stimV3;stimV4;stimV5];
% stimC = [stimC1;stimC2;stimC3;stimC4;stimC5];
% stimD = [stimD1;stimD2;stimD3;stimD4;stimD5];
% r = [r1;r2;r3;r4;r5];

t = [t1 t2 t3 t4 t5 t6 t7 t9 t10];

stimV = [stimV6;stimV7;stimV9;stimV10];
stimC = [stimC6;stimC7;stimC9;stimC10];
stimD = [stimD6;stimD7;stimD9;stimD10];
r = [r6;r7;r9;r10];


%% Run prior/likelihood estimation

xgrid = getLogXform(linspace(0,16,100),0.3);

geoOn = 0;
% geoOn = 1;

numStartPts = 5;

% solverOpts.fxnForm = 'numerical';
solverOpts.fxnForm = 'analytical';

solverOpts.solver = 'fmincon';
% solverOpts.solver = 'bads';

% Run for single Gaussian
nB = 1;
[outMatSG] = fit2AFCData_TSM(stimV,stimC,stimD,r,nB,xgrid,geoOn,numStartPts,solverOpts);

% Run for MoG
nB = 2;
[outMatMoG] = fit2AFCData_TSM(stimV,stimC,stimD,r,nB,xgrid,geoOn,numStartPts,solverOpts);

% Save fits
save(['/home/tyler/Documents/MATLAB/cooperLab/3-Data/priorCoordsRDS/TSM/outMatMog_geo',num2str(geoOn)],'outMatMoG');
save(['/home/tyler/Documents/MATLAB/cooperLab/3-Data/priorCoordsRDS/TSM/outMatSG_geo',num2str(geoOn)],'outMatSG');


%% Plot psychometric functions with Bayesian and Gaussian fits
% expType = 'within';
expType = 'between';
% expType = 'sizeControl';
[stimStruct,figs] = fitPsychFxns(stimV,stimC,stimD,r,outMatSG,outMatMoG,expType);


%% Plot Bayesian parameters 
plotBayesFit(outMatSG,outMatMoG,stimStruct,geoOn);


toc