% close all
clear all

%% Load in data

datDir = '/media/tyler/Data/MATLAB/cooperLab/3-Data/priorCoordsRDS/';
% datStr = 'tyler/tyler-20220311T170416-1.mat';
% datStr = 'tyler/tyler-20220314T165526-1.mat';
% datStr = 'CTW/CTW-20220316T144044-1.mat';
% datStr = 'CTW/CTW-20220316T151417-2.mat';
% datStr = 'tyler/tyler-20220318T182022-1.mat';
datStr = 'tyler/tyler-20220321T130518-1.mat';

load([datDir,datStr]);

% Response Mat
respMat       = pa.response;

% Sort data into stimulus conditions
numStair      = pa.numStaircases;

[contCombs,contInd,~] = unique(respMat(:,2));       % contrast combinations
[testScr,~,~]         = unique(respMat(:,6));       % presentation screen
[refVels,refVInd,~]   = unique(respMat(:,5));       % reference velocity
%%%%%% can take this part out in future
stairTypes = [1;2;3];
scInds     = makeCombos([numel(contInd) numel(testScr) numel(refVels) size(stairTypes,1)]);

for ii = 1:numStair
   
    respMat(respMat(:,10)==ii,13) = scInds(ii,4);
    
end
%%%%%%%%%%
stairTypes            = unique(respMat(:,13));      % staircase target %

% Group staircase triplets together
stimInds     = makeCombos([numel(contInd) numel(testScr) numel(refVels)]);
numUniqStim  = size(stimInds,1);

for ii = 1:numUniqStim
    
    theseInds = (respMat(:,2) == contCombs(stimInds(ii,1))) & ...
                (respMat(:,6) == testScr(stimInds(ii,2))) & ...
                (respMat(:,5) == refVels(stimInds(ii,3)));
    
    resps{ii} = respMat(theseInds,:);
    
end

% Which stimulus to plot?
stim         = 1; %[1 2 3 4]
theseResps   = resps{stim};


%% Plot test speeds for each staircase trial

stairVels    = cell(numel(stairTypes),1);
reversalInds = cell(numel(stairTypes),1);
thresh       = nan(numStair,1);
numTrials    = nan(numStair,1);
p            = nan(numStair,1);

for sid = 1:numel(stairTypes)
    
    % Get indices of all trials run for this staircase
%     stairInds{sid}  = find(theseResps(:,13) == sid);
    stairVels{sid}  = theseResps(theseResps(:,13) == sid,3);
%     stairResps{sid} = theseResps(theseResps(:,13) == sid,8) ~= theseResps(theseResps(:,13) == sid,9);
    
    temp1           = find(theseResps(:,13) == sid);
    thisStairInd    = theseResps(temp1(1),10);

%     trialInds      = [1:numel(stairVels{sid})]';
%     
%     % Find trial indices of all reversals
%     temp1          = [nan;diff(stairVels{sid})];                        % Find differences between trial speeds
%     temp2          = (temp1~=0) & ~isnan(temp1);                        % Find where speed actually changes
%     temp3          = [[nan;diff(sign(temp1(temp2)))] trialInds(temp2)]; % Find when the sign of these changes change
%     temp4          = temp3(abs(temp3(:,1)) > 0,2);                      
%     
%     reversalInds{sid} = trialInds(temp4);
%     reversalInds{sid} = find(abs([nan;diff(stairResps{sid})])>0);
    reversalInds{sid} = logical([get(gSCell{thisStairInd},'reversals')]);
    
    % Calculate threshold as mean of speeds after 1st reversal
    temp1       = find(reversalInds{sid}==1);
    thresh(sid) = mean(stairVels{sid}(temp1(end-5:end)));
    
    % Get number of trials run for this staircase
    numTrials(sid) = numel(stairVels{sid});
    
end

% Plotting pars
speedTick    = [0.1 0.5 1 2 4 8 12];
% speedLims    = [min([theseResps(:,3);theseResps(1,5)]) round(max([theseResps(:,3);theseResps(1,5)]))];
speedLims    = [min([theseResps(:,3);theseResps(1,5)]) 15];

colorMat     = colororder;

% Plot test speeds for all trials run in this staircase
f1 = figure;
f1.Position = [100 100 1900 500];
subplot(1,4,1);
hold on;

cind = 1;

for sid = 1:numel(stairTypes)
    % Plot all trials
    p(cind)     = scatter(1:numel(stairVels{sid}),stairVels{sid},30,colorMat(sid,:));
    
    % Plot reversal trials in filled dots
    trialInds = 1:numel(stairVels{sid});
    scatter(trialInds(reversalInds{sid}),stairVels{sid}(reversalInds{sid}),30,colorMat(sid,:),'filled');
    plot([0 max(numTrials)],thresh(sid)*[1 1],'color',colorMat(cind,:));
    cind = cind + 1;
end

p(4) = plot([0 max(numTrials)],theseResps(1,5)*[1 1],'--k');

switch theseResps(1,5)
    case 8
        locale = 'southeast';
    case 1
        locale = 'northeast';
end

legend([p(1) p(2) p(3) p(4)],{'1u/1d','2u/1d','1u/2d','Ref. vel.'},'location',locale);
xlabel('Trial');
ylabel('Test speed (deg/s)');
title(['V_{ref}=',num2str(theseResps(1,5)),'\circ/s; C_{ref}=',num2str(theseResps(1,2)),'; Dist.=',num2str(ds.screenDistance(theseResps(1,6))),'m']);
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'yscale','log',...
        'ytick',speedTick,'ylim',speedLims,'xlim',[0 max(numTrials)]);

    
%% Plot responses and psychometric functions

% Recode responses in terms of whether test stimulus was seen faster
[uniqueVel,~,uniInds] = unique(theseResps(:,3));
respDir   = theseResps(:,9);
refDir    = theseResps(:,8);

respStim  = respDir ~= refDir;

% Calculate proportion of "test faster" responses for each unique velocity
respP  = nan(numel(uniqueVel),1);
numInp = nan(numel(uniqueVel),1);

for ii = 1:numel(uniqueVel)
    
    inputs     = respStim(uniInds == ii);
    respP(ii)  = sum(inputs)/numel(inputs);
    numInp(ii) = numel(inputs);
    
end

% Fit a cumulative Gaussian to the data
[pFxn,PSE,sens] = fitCumGauss(uniqueVel,respP,numInp/sum(numInp));

% Fit to bootstrapped data
numRS = 40;

trueData = [theseResps(:,3) respStim];

numSampSize = 10;
sampSize = round(size(trueData,1)*linspace(0.1,1,numSampSize));

pFxnrs = cell(numRS,1);
PSErs  = nan(numRS,1);
sensrs = nan(numRS,1);

for ii = 1:numRS
    
    respsRS   = datasample(trueData,size(trueData,1),'Replace',true);
    [pFxnrs{ii},PSErs(ii),sensrs(ii)] = fitCumGauss(respsRS(:,1),respsRS(:,2),ones(size(respsRS,1),1),uniqueVel);
    
    for jj = 1:numSampSize
        respsRS2   = datasample(trueData,sampSize(jj),'Replace',true);
        [~,PSErs2(ii,jj),sensrs2(ii,jj)] = fitCumGauss(respsRS2(:,1),respsRS2(:,2),ones(size(respsRS2,1),1),uniqueVel);
    end
    
end

% Plot fit and resampled fits
subplot(1,4,2);
hold on;

for ii = 1:numRS
    
    plot(uniqueVel,pFxnrs{ii},'color',[0.65 0.65 0.65]);
    
end

plot(uniqueVel,pFxn,'k','linewidth',2);

% Calculate where each staircase would fall on psychometric function 
s1 = 0.5;
s2 = 1 - sqrt(0.5);
s3 = sqrt(0.5);

% Plot all response proportions,thresholds fit to staircases, and psychometric functions
scatter(uniqueVel,respP,30,'k');
scatter(thresh(1),[s1],100,colorMat(1,:),'filled');
scatter(thresh(2),[s2],100,colorMat(2,:),'filled');
scatter(thresh(3),[s3],100,colorMat(3,:),'filled');
plot(theseResps(1,5)*[1 1],[0 1],'--k');
set(gca,'xscale','log','fontsize',20,'xtick',speedTick,'xlim',speedLims);
ylabel('p("test faster")');
xlabel('Test speed (deg/s)');
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);


%% Plot estimates of PSE and slope with CIs

pseCI = 1.96*std(PSErs2,0,1)./sqrt(numRS);
sensCI = 1.96*std(sensrs2,0,1)./sqrt(numRS);

PSE2 = PSE*ones(numSampSize,1);
sens2 = sens*ones(numSampSize,1);

subplot(1,4,3);
errorbar(sampSize,PSE2,pseCI,'s','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
title('PSE');
set(gca,'ytick',speedTick,'ylim',speedLims,'yscale','log','fontsize',20,'plotboxaspectratio',[1 1 1]);
xlabel('Trial count');

subplot(1,4,4);
errorbar(sampSize,sens2,sensCI,'s','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
title('Slope');
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);
xlabel('Trial count');


