function [f1] = plotThisData(subjID,exp,stim)

% Plot staircases and psychometric curve from one stimlus condition


%% Load in and sort data
% Set up directories
splPath = regexp(which('plotThisData'),filesep,'split');

dataDir = [filesep,fullfile(splPath{1:numel(splPath)-4}),filesep,'3-Data/priorCoordsRDS/'];

load([dataDir,'xSub',filesep,subjID,'.mat'],'stimStruct');

% cull empty fields
incInds    = arrayfun(@(x) ~isempty(x.testCont),stimStruct);
stimStruct = stimStruct(incInds);


%% Plot condition of interest

% Find indices that match input arguments
for ii = 1:numel(stimStruct)
%     expMatch(ii,1) = strcmp(exp,stimStruct(ii).blockType);
    expMatch(1,ii) = strcmp(exp,stimStruct(ii).blockType);
end

testConts = arrayfun(@(x) x.testCont,stimStruct);
% testSz = arrayfun(@(x) x.testSz,stimStruct);
testDists = arrayfun(@(x) x.testDist,stimStruct);
refVels   = arrayfun(@(x) x.refVel,stimStruct);

%%%% need to fix expmatch dimensions when > single dataset, needs
%%%% transposing
stimInd   = (testConts == stim(1)) & (testDists == stim(2)) & (refVels == stim(3)) & expMatch;
% stimInd   = (testConts == stim(1)) & (testDists == stim(2)) & (testDists == stim(3)) & (refVels == stim(4)) & expMatch;

thisData = stimStruct(stimInd);

% Grab info for plotting
numStairTypes = numel(thisData.stairVels);
numTrials     = cellfun(@(x) numel(x),thisData.stairVels);

% Plotting pars
sTickLin = [0.1 0.5 1 2 4 8 12 16 20];
speedTick = getLogXform(sTickLin,0.3);
for ii = 1:numel(sTickLin)
    sTickLab{ii} = sTickLin(ii);
end
colorMat  = colororder;


f1 = figure;
f1.Position = [100 100 1100 550];
hold on;

% Plot staircase values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1);
hold on;

for sid = 1:numStairTypes
    
    theseStairVels = getLogXform(thisData.stairVels{sid},0.3);
    lims(sid,:)    = [min(theseStairVels) max(theseStairVels)];
    
    % Plot all trials
    p(sid)    = scatter(1:numel(theseStairVels),theseStairVels,30,colorMat(sid,:));
    
    % Plot reversal trials in filled dots
    trialInds = 1:numel(theseStairVels);
    scatter(trialInds(thisData.reversalInds{sid}),theseStairVels(thisData.reversalInds{sid}),30,colorMat(sid,:),'filled');
    plot([0 max(numTrials)],getLogXform(thisData.thresh(sid),0.3)*[1 1],'color',colorMat(sid,:));
    
end

% Define speed limits for plotting
if stim(3) < 4
    speedLims    = getLogXform([0.1 8],0.3);
else 
    speedLims    = getLogXform([2 20],0.3);
end
% speedLims = [min(lims(:,1)) max(lims(:,2))];

xlabel('Trial');
ylabel('Test speed (deg/s)');
title(['V_{ref}=',num2str(stim(3)),'\circ/s; C_{test}=',num2str(stim(1)),...
       '; D_{test}=',num2str(stim(2)),'m']);
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],...
        'ytick',speedTick,'yticklabel',sTickLab,'ylim',speedLims,'xlim',[0 max(numTrials)]);

    
% Plot response probabilities and psychometric function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    theseStim = thisData.respData(:,1);
    respP  = thisData.respData(:,2)./thisData.respData(:,3);
    numInp = thisData.respData(:,3);
else
    % cull repeated entries (where are they coming from? eliminate.)
    [~,keepThese,~] = unique(round(thisData.respData(:,1),4));
    thisData.respData = thisData.respData(keepThese,:);
    thisData.pFxn = thisData.pFxn(keepThese,:);
    
    % bin/combine responses a bit
    stimEdges = linspace(thisData.respData(1,1),thisData.respData(end,1)+min(diff(thisData.respData(:,1))),13);
    theseStim = stimEdges(1:(numel(stimEdges)-1))+diff(stimEdges(1:2))/2;
    for ii = 1:numel(stimEdges)-1
        theseInds  = (thisData.respData(:,1) >= stimEdges(ii)) & (thisData.respData(:,1) < stimEdges(ii+1));
        
        respP(ii)  = sum(thisData.respData(theseInds,2))./sum(thisData.respData(theseInds,3));
        numInp(ii) = sum(thisData.respData(theseInds,3));
        
        if isnan(respP(ii))
            respP(ii) = 0;
        end
    end
    
    % Cull bins with no trials
    respP     = respP(numInp ~= 0);
    theseStim = theseStim(numInp ~= 0);
    numInp    = numInp(numInp ~= 0);
end

subplot(1,2,2);
hold on;

plot(thisData.respData(:,1),thisData.pFxn,'k','linewidth',2);

% % Calculate where each staircase would fall on psychometric function
% s1 = 0.5;
% s2 = sqrt(0.5);
% s3 = 1 - sqrt(0.5);

% Plot all response proportions,thresholds fit to staircases, and psychometric functions
scatter(theseStim,respP,30*numInp,'k');
% scatter(getLogXform(thisData.thresh(1),0.3),[s1],100,colorMat(1,:),'filled');
% scatter(getLogXform(thisData.thresh(2),0.3),[s2],100,colorMat(2,:),'filled');
% scatter(getLogXform(thisData.thresh(3),0.3),[s3],100,colorMat(3,:),'filled');

% Plot reference
plot(getLogXform(stim(3),0.3)*[1 1],[0 1],'--k');

% Plot PSE and confidence intervals
plot(thisData.PSE*[1 1],[0 0.5],'-.k','linewidth',2);
errorbar(thisData.PSE,0.5,[],[],thisData.PSE-thisData.PSECI(1),thisData.PSE-thisData.PSECI(2),'-.k','linewidth',2);

% Plot slope CIs
yLCI = thisData.slopeCI(2)*([-1 1]/4) + 0.5;
yUCI = thisData.slopeCI(1)*([-1 1]/4) + 0.5;
ySlope = thisData.slope*([-1 1]/4) + 0.5;
fill([-1 1 1 -1 -1]/4 + thisData.PSE,[yLCI fliplr(yUCI) yLCI(1)],'k','facealpha',0.25,'edgecolor','none');
plot([-1 1]/4 + thisData.PSE,ySlope,'k');

set(gca,'fontsize',20,'xtick',speedTick,'xticklabel',sTickLab,'xlim',speedLims,'ylim',[0 1]);
ylabel('p("test faster")');
xlabel('Test speed (deg/s)');
title(exp);
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);

end