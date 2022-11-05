function [stats] = analyzeXsubRDS2(saveOn,varargin)

% Compare prior coords experiment across subjects (5x ref speeds, 4x simulated far dists)

%% Load in and sort data
% Set up directories
splPath = regexp(which('analyzeXsubRDS'),filesep,'split');

dataDir = [filesep,fullfile(splPath{1:numel(splPath)-4}),filesep,'3-Data/priorCoordsRDS/xSub/'];

if nargin == 0
    % Find datsets in subject directory
    files     = dir(dataDir);
    fileNames = arrayfun(@(x)x.name,files,'uniformoutput',false);
    datFiles  = fileNames( cellfun(@(x) contains(x,'mat'),fileNames) );
else
    % Specify which subjects to include in analysis
    datFiles  = varargin{1};
end

% Load in data and assign block names
for ii = 1:numel(datFiles)
    
    % make array of subjIDs
    if numel(datFiles{ii})>3
        subjID{ii} = datFiles{ii}(1:3);
    else
        subjID{ii} = datFiles{ii};
    end
    
    data(ii)            = load([dataDir,datFiles{ii}],'stimStruct');
    data(ii).stimStruct = data(ii).stimStruct(:);
    
    % Cull empty fields
    incInds = arrayfun(@(x) ~isempty(x.respMat),data(ii).stimStruct);
    data(ii).stimStruct = data(ii).stimStruct(incInds);
    
    % Reorganize fields into standard order
    blockType = arrayfun(@(x) x.blockType,data(ii).stimStruct(:),'uniformoutput',false);
    testCont  = arrayfun(@(x) x.testCont,data(ii).stimStruct(:));
    testDist  = arrayfun(@(x) x.testDist,data(ii).stimStruct(:));
    refDist   = arrayfun(@(x) x.refDist,data(ii).stimStruct(:));
    trueTestDist  = arrayfun(@(x) x.trueTestDist,data(ii).stimStruct(:));
    trueRefDist   = arrayfun(@(x) x.trueRefDist,data(ii).stimStruct(:));
    refVel    = arrayfun(@(x) x.refVel,data(ii).stimStruct(:));
    
    blockID = nan(numel(blockType),1);
    
    for jj = 1:numel(blockType)
        
        if strcmp(blockType{jj},'within')
            blockID(jj) = 1;
        elseif strcmp(blockType{jj},'between')
            blockID(jj) = 2;
        elseif strcmp(blockType{jj},'sizeCntl')
            blockID(jj) = 3;
        end
        
    end
    
    theseStimInds = [blockID testCont testDist refDist refVel trueTestDist trueRefDist];
    [~,sortInds]  = sortrows(theseStimInds);

    data(ii).stimStruct = data(ii).stimStruct(sortInds);
    stimInds{ii} = theseStimInds;

end

% Sort array of stim parameters
%%%%% This fails if the last subject doesn't have the proper
%%%%% number of conditions (29)
testCont = testCont(sortInds);
testDist = testDist(sortInds);
refDist  = refDist(sortInds);
trueTestDist = trueTestDist(sortInds);
trueRefDist  = trueRefDist(sortInds);
refVel   = refVel(sortInds);
blockID  = blockID(sortInds);

stimArr  = [blockID testCont testDist refDist refVel trueTestDist trueRefDist];

numSubs = numel(data);

% Make pairings 
% (near vs. far for each ref vel, contrast polarity, trial type)

% Within
withinTrials = blockID == 1;
withinDists  = unique(testDist(withinTrials));
withinPairs  = [find(testDist(withinTrials) == withinDists(1)) find(testDist(withinTrials) == withinDists(2))];

for ii = 1:size(withinPairs,1)
    pairPars_within{ii,1} = [blockID(withinPairs(ii,:)) testCont(withinPairs(ii,:)) testDist(withinPairs(ii,:)) refVel(withinPairs(ii,:))];
end

% Between
betweenPairs  = [find((testDist == 0.5) & (refDist ~= 0.5) & (blockID == 2)), ...
                 find((testDist ~= 0.5) & (refDist == 0.5) & (blockID == 2))];
betweenPairs  = [betweenPairs;...
                 find((testDist == 0.5) & (refDist == 0.5) & (trueRefDist ~= 0.5) & (blockID == 2)), ...
                 find((testDist == 0.5) & (refDist == 0.5) & (trueRefDist == 0.5) & (blockID == 2))];

for ii = 1:size(betweenPairs,1)
    pairPars_between{ii,1} = [blockID(betweenPairs(ii,:)) testDist(betweenPairs(ii,:)) refDist(betweenPairs(ii,:))];
end

% Size cntl
cntlTrial     = blockID == 3;

% Grab best fit PSEs and slopes
PSEs    = nan(numel(blockID),numSubs);
PSECILB = nan(numel(blockID),numSubs);
PSECIUB = nan(numel(blockID),numSubs);
slopes    = nan(numel(blockID),numSubs);
slopeCILB = nan(numel(blockID),numSubs);
slopeCIUB = nan(numel(blockID),numSubs);

for ii = 1:numSubs

    % Select dataset
    thisDataset = data(ii).stimStruct(:);
    
    if size(thisDataset,1) == numel(blockID)
        % Grab PSEs
        PSEs(:,ii)    = arrayfun(@(x) x.PSE,thisDataset);
        PSECILB(:,ii) = arrayfun(@(x) x.PSECI(1),thisDataset);
        PSECIUB(:,ii) = arrayfun(@(x) x.PSECI(2),thisDataset);

        % Grab slopes
        slopes(:,ii)    = arrayfun(@(x) x.slope,thisDataset);
        slopeCILB(:,ii) = arrayfun(@(x) x.slopeCI(1),thisDataset);
        slopeCIUB(:,ii) = arrayfun(@(x) x.slopeCI(2),thisDataset);

        % ID bad fits
        badFits(:,ii) = arrayfun(@(x) x.badFit,thisDataset);
    else
        % For subjects that don't have fits for all conditions

        % Parameters for this sub's unique trial types
        thesePars = stimInds{ii};
        theseInds = nan(size(thesePars,1),1);

        % Find how this sub's trial type inds match up with the group's inds
        for jj = 1:size(thesePars,1)
            theseInds(jj) = find(sum(abs(thesePars(jj,:) - stimArr),2)==0);
        end

        % Grab PSEs
        PSEs(theseInds,ii)    = arrayfun(@(x) x.PSE,thisDataset);
        PSECILB(theseInds,ii) = arrayfun(@(x) x.PSECI(1),thisDataset);
        PSECIUB(theseInds,ii) = arrayfun(@(x) x.PSECI(2),thisDataset);

        % Grab slopes
        slopes(theseInds,ii)    = arrayfun(@(x) x.slope,thisDataset);
        slopeCILB(theseInds,ii) = arrayfun(@(x) x.slopeCI(1),thisDataset);
        slopeCIUB(theseInds,ii) = arrayfun(@(x) x.slopeCI(2),thisDataset);

        % ID bad fits
        badFits(theseInds,ii) = arrayfun(@(x) x.badFit,thisDataset);
    end
end 

%% Detect and remove bad fits

numBadFits = sum(badFits);
numFits    = sum(~isnan(PSEs));
percBF     = numBadFits./numFits;

PSEs(badFits)      = nan;
PSECILB(badFits)   = nan;
PSECIUB(badFits)   = nan;
slopes(badFits)    = nan;
slopeCILB(badFits) = nan;
slopeCIUB(badFits) = nan;


%% Run group stats

%%%% note: everything we're sorting here is still log-xformed, need to
%%%% xform back before saving stats?

% Within
%-------------------------%
% PSE
meanPSEs_within = mean(PSEs(withinPairs(:),:),2,'omitnan');
stdPSEs_within  = std(PSEs(withinPairs(:),:),[],2,'omitnan');
meanPSEs_within = reshape(meanPSEs_within,size(withinPairs));
stdPSEs_within  = reshape(stdPSEs_within,size(withinPairs));

[~,~,groupPSECI_within{1}]  = ttest(PSEs(withinPairs(:,1),:),[],'dim',2);
[~,~,groupPSECI_within{2}]  = ttest(PSEs(withinPairs(:,2),:),[],'dim',2);
[~,pNormPSE_within]         = ttest(PSEs(withinPairs(:,1),:),PSEs(withinPairs(:,2),:),'dim',2);

for ii = 1:size(withinPairs,1)
    pSRPSE_within(ii)       = signrank(PSEs(withinPairs(ii,1),:),PSEs(withinPairs(ii,2),:));
end

% Slope
meanSlopes_within = mean(slopes(withinPairs(:),:),2,'omitnan');
stdSlopes_within  = std(slopes(withinPairs(:),:),[],2,'omitnan');
meanSlopes_within = reshape(meanSlopes_within,size(withinPairs));
stdSlopes_within  = reshape(stdSlopes_within,size(withinPairs));

[~,~,groupSlopesCI_within{1}]  = ttest(slopes(withinPairs(:,1),:),[],'dim',2);
[~,~,groupSlopesCI_within{2}]  = ttest(slopes(withinPairs(:,2),:),[],'dim',2);
[~,pNormSlopes_within]         = ttest(slopes(withinPairs(:,1),:),slopes(withinPairs(:,2),:),'dim',2);

for ii = 1:size(withinPairs,1)
    try
    pSRSlopes_within(ii)       = signrank(slopes(withinPairs(ii,1)),slopes(withinPairs(ii,2)));
    catch
    pSRSlopes_within(ii)       = nan;
    end
end


% Between
%-------------------------%
betweenDists = refDist(betweenPairs);
betweenTrueDists = trueRefDist(betweenPairs);
betweenRefV  = refVel(betweenPairs(1));

[~,betweenSort] = sort(betweenDists(:,1));
betweenSort = flipud(betweenSort);

betweenDists     = betweenDists(betweenSort,:);
betweenTrueDists = betweenTrueDists(betweenSort,:);
betweenPairs     = betweenPairs(betweenSort,:);

numBetPairs = size(betweenPairs,1);

% PSE
meanPSEs_between = mean(PSEs(betweenPairs(:),:),2,'omitnan');
stdPSEs_between  = std(PSEs(betweenPairs(:),:),[],2,'omitnan');
meanPSEs_between = reshape(meanPSEs_between,size(betweenPairs));
stdPSEs_between  = reshape(stdPSEs_between,size(betweenPairs));

[~,~,groupPSECI_between{1}]  = ttest(PSEs(betweenPairs(:,1),:),[],'dim',2);
[~,~,groupPSECI_between{2}]  = ttest(PSEs(betweenPairs(:,2),:),[],'dim',2);
[~,pNormPSE_between]         = ttest(PSEs(betweenPairs(:,1),:),PSEs(betweenPairs(:,2),:),'dim',2);

for ii = 1:size(betweenPairs,1)-1
    pSRPSE_between(ii)         = signrank(PSEs(betweenPairs(ii,1),:),PSEs(betweenPairs(ii,2),:));
end


% Size Control
%-------------------------%
meanPSEs_sc = mean(PSEs(cntlTrial,:),2,'omitnan');
stdPSEs_sc  = std(PSEs(cntlTrial,:),[],2,'omitnan');
[~,pNormPSE_sc,groupPSECI_sc] = ttest(PSEs(cntlTrial,:));


% Output stats
stats.meanPSEs_within      = meanPSEs_within;
stats.stdPSEs_within       = stdPSEs_within;
stats.groupPSECI_within    = groupPSECI_within;
stats.pNormPSE_within      = pNormPSE_within;
stats.pSRPSE_within        = pSRPSE_within;

stats.meanSlopes_within    = meanSlopes_within;
stats.stdSlopes_within     = stdSlopes_within;
stats.groupSlopesCI_within = groupSlopesCI_within;
stats.pNormSlopes_within   = pNormSlopes_within;
stats.pSRSlopes_within     = pSRSlopes_within;

stats.meanPSEs_between     = meanPSEs_between;
stats.stdPSEs_between      = stdPSEs_between;
stats.groupPSECI_between   = groupPSECI_between;
stats.pNormPSE_between     = pNormPSE_between;
stats.pSRPSE_between       = pSRPSE_between;

stats.meanPSEs_sc          = meanPSEs_sc;
stats.stdPSEs_sc           = stdPSEs_sc;
stats.groupPSECI_sc        = groupPSECI_sc;
stats.pNormPSE_sc          = pNormPSE_sc;


%% Plot 

% Setup some axes vars
spread   = linspace(-0.1,0.1,numSubs);
% colorMat = colororder;
ytixLin  = [0.001 0.1 0.5 1 2 4 8 16 32];
ytix     = getLogXform(ytixLin,0.3);
for ii = 1:numel(ytixLin)
    ytixlabel{ii} = num2str(ytixLin(ii));
end

yLims = getLogXform([0.1 32],0.3);

withinRefV  = unique(refVel);
withinDists = unique(testDist(withinPairs));
withinConts = unique(testCont(withinPairs));

subColor    = 0.75*[1 1 1];

slopeMax    = ceil(max(slopes(:)));

% Within
%-------------------------%
% PSEs
f1 = figure;
f1.Position = [2000 100 2400 400];

% Loop over ref vels
for ii = 1:numel(withinRefV)

    subplot(1,5,ii);
    hold on

    % Plot line for reference velocity
    thisRefVel = withinRefV(ii);
    plot([0 3],getLogXform(thisRefVel,0.3)*[1 1],'--k')

    % Plot individual subjects
    for jj = 1:numSubs    % Loop over subjects
        for hh = 1:2    % Loop over contrast polarity
            theseRows   = (thisRefVel == refVel(withinPairs(:,1))) & (withinConts(hh) == testCont(withinPairs(:,1)));
            thesePSEs   = PSEs(withinPairs(theseRows,:),jj);
            thesePSELBs = PSECILB(withinPairs(theseRows,:),jj);
            thesePSEUBs = PSECIUB(withinPairs(theseRows,:),jj);
            thisDist    = withinDists;

            errorbar(thisDist + spread(jj),thesePSEs,thesePSEs-thesePSELBs,thesePSEUBs-thesePSEs,'color',subColor,'linewidth',2);

            if hh == 1
                scatter(thisDist + spread(jj),thesePSEs,150,subColor,'filled');
            else
                scatter(thisDist + spread(jj),thesePSEs,150,subColor,'linewidth',2);
            end
        end
    end

    % Plot group mean and STD
    for hh = 1:2    % Loop over contrast polarity 
        theseRows       = (thisRefVel == refVel(withinPairs(:,1))) & (withinConts(hh) == testCont(withinPairs(:,1)));
        meanPSEs        = meanPSEs_within(theseRows,:);
        groupPSECI(1,:) = [groupPSECI_within{1}(theseRows,1) groupPSECI_within{2}(theseRows,1)];
        groupPSECI(2,:) = [groupPSECI_within{1}(theseRows,2) groupPSECI_within{2}(theseRows,2)];

        errorbar(withinDists,meanPSEs,meanPSEs-groupPSECI(1,:),groupPSECI(2,:)-meanPSEs,'color',[0 0 0],'linewidth',2);

        if hh == 1
            scatter(withinDists,meanPSEs,200,[0 0 0],'filled');
        else
            scatter(withinDists,meanPSEs,200,[0 0 0],'linewidth',2);
        end
    end
    
    set(gca,'fontsize',20,'xlim',[0 1.5],'xtick',[0.5 1],'ytick',ytix,'yticklabel',ytixlabel,'ylim',yLims);
    xlabel('Viewing Distance (m)');
    ylabel('PSE (\circ/s)');
    title(['V_{ref} = ',num2str(thisRefVel)]);
%     text(0.25,0.2*abs(diff(theseyLims))+theseyLims(1),['p = ',num2str(round(pNormPSE(hh),5))],'fontsize',20);

end


% Slopes
f2 = figure;
f2.Position = [2000 600 2400 400];

% Loop over ref vels
for ii = 1:numel(withinRefV)

    subplot(1,5,ii);
    hold on

    thisRefVel = withinRefV(ii);

    % Plot individual subjects
    for jj = 1:numSubs    % Loop over subjects
        for hh = 1:2    % Loop over contrast polarity
            theseRows     = (thisRefVel == refVel(withinPairs(:,1))) & (withinConts(hh) == testCont(withinPairs(:,1)));
            theseSlopes   = slopes(withinPairs(theseRows,:),jj);
            theseSlopeLBs = slopeCILB(withinPairs(theseRows,:),jj);
            theseSlopeUBs = slopeCIUB(withinPairs(theseRows,:),jj);
            thisDist      = withinDists;

            %%%%% gotta do errorbar per person if you want to connect lines
            errorbar(thisDist + spread(jj),theseSlopes,theseSlopes-theseSlopeLBs,theseSlopeUBs-theseSlopes,'color',subColor,'linewidth',2);

            if hh == 1
                scatter(thisDist + spread(jj),theseSlopes,150,subColor,'filled');
            else
                scatter(thisDist + spread(jj),theseSlopes,150,subColor,'linewidth',2);
            end
        end
    end

    % Plot group mean and STD
    for hh = 1:2    % Loop over contrast polarity 
        theseRows          = (thisRefVel == refVel(withinPairs(:,1))) & (withinConts(hh) == testCont(withinPairs(:,1)));
        meanSlopes         = meanSlopes_within(theseRows,:);
        groupSlopesCI(1,:) = [groupSlopesCI_within{1}(theseRows,1) groupSlopesCI_within{2}(theseRows,1)];
        groupSlopesCI(2,:) = [groupSlopesCI_within{1}(theseRows,2) groupSlopesCI_within{2}(theseRows,2)];

        errorbar(withinDists,meanSlopes,meanSlopes-groupSlopesCI(1,:),groupSlopesCI(2,:)-meanSlopes,'color',[0 0 0],'linewidth',2);

        if hh == 1
            scatter(withinDists,meanSlopes,200,[0 0 0],'filled');
        else
            scatter(withinDists,meanSlopes,200,[0 0 0],'linewidth',2);
        end
    end
    
    set(gca,'fontsize',20,'xlim',[0 1.5],'xtick',[0.5 1],'ylim',[0 slopeMax]);
    xlabel('Viewing Distance (m)');
    ylabel('Slope (\Delta p/\Delta log(V))');
    title(['V_{ref} = ',num2str(thisRefVel)]);
%     text(0.25,0.2*abs(diff(theseyLims))+theseyLims(1),['p = ',num2str(round(pNormPSE(hh),5))],'fontsize',20);

end


% Between
%-------------------------%
f3 = figure;
f3.Position = [2000 100 2000 400];

% Loop over simulated distances
for ii = 1:4

    theseDists = betweenTrueDists(ii,:);
    simDist    = betweenDists(ii,1);

    subplot(1,4,ii);
    hold on

    % Plot lines for reference velocity and double/half velocity
    thisRefVel = betweenRefV;
    plot([0 3],getLogXform(thisRefVel,0.3)*[1 1],'--k')
    plot([0 3],getLogXform(thisRefVel*0.5,0.3)*[1 1],'--k')
    plot([0 3],getLogXform(thisRefVel*2,0.3)*[1 1],'--k')

    % Plot individual subjects
    for jj = 1:numSubs    % Loop over subjects
            
        theseRows   = betweenPairs(ii,:);
        thesePSEs   = PSEs(theseRows,jj);
        thesePSELBs = PSECILB(theseRows,jj);
        thesePSEUBs = PSECIUB(theseRows,jj);
        thisDist    = theseDists;

        errorbar(thisDist + spread(jj),thesePSEs,thesePSEs-thesePSELBs,thesePSEUBs-thesePSEs,'color',subColor,'linewidth',2);

        scatter(thisDist + spread(jj),thesePSEs,150,subColor,'filled');

    end
    
    % Plot group mean and STD
    theseRows       = ii;
    meanPSEs        = meanPSEs_between(theseRows,:);
    groupPSECI(1,:) = [groupPSECI_between{1}(theseRows,1) groupPSECI_between{2}(theseRows,1)];
    groupPSECI(2,:) = [groupPSECI_between{1}(theseRows,2) groupPSECI_between{2}(theseRows,2)];

    errorbar(theseDists,meanPSEs,meanPSEs-groupPSECI(1,:),groupPSECI(2,:)-meanPSEs,'color',[0 0 0],'linewidth',2);

    scatter(theseDists,meanPSEs,200,[0 0 0],'filled');

    set(gca,'fontsize',20,'xlim',[0 1.5],'xtick',[0.5 1],'ytick',ytix,'yticklabel',ytixlabel,'ylim',yLims);
    xlabel('Reference Distance (m)');
    ylabel('PSE (\circ/s)');
    title(['D_{ref,sim} = ',num2str(simDist)]);

end


% Size Control
%-------------------------%
f4 = figure;
f4.Position = [100 100 450 400];
hold on;

% Plot line for reference velocity
thisRefVel = betweenRefV;
plot([0 3],getLogXform(thisRefVel,0.3)*[1 1],'--k')
thisDist    = 0.5;

for jj = 1:numSubs    % Loop over subjects

    thisRow     = cntlTrial;

    if sum(thisRow)>0
        thesePSEs   = PSEs(thisRow,jj);
        thesePSELBs = PSECILB(thisRow,jj);
        thesePSEUBs = PSECIUB(thisRow,jj);

        errorbar(thisDist + spread(jj),thesePSEs,thesePSEs-thesePSELBs,thesePSEUBs-thesePSEs,...
            'color',subColor,'linewidth',2);

        scatter(thisDist + spread(jj),thesePSEs,150,subColor,'filled');
    end
end
    
% Plot group mean and STD
meanPSEs   = meanPSEs_sc;
groupPSECI = groupPSECI_sc;

errorbar(thisDist,meanPSEs,meanPSEs-groupPSECI(1),groupPSECI(2)-meanPSEs,'color',[0 0 0],'linewidth',2);

scatter(thisDist,meanPSEs,200,[0 0 0],'filled');

set(gca,'fontsize',20,'xlim',[0.25 0.75],'xtick',[],'ytick',ytix,'yticklabel',ytixlabel,'ylim',yLims);
% xlabel('Reference Distance (m)');
ylabel('PSE (\circ/s)');
title('Size Control (Diam. Ref > test)');


% Plot how many psych functions culled per subject
cutoff = 1/3;
overcutoff  = find(percBF >= cutoff);
undercutoff = find(percBF < cutoff);

f5 = figure;
f5.Position = [100 100 1450 450];
hold on;

bar(overcutoff,percBF(overcutoff),'facecolor',[0.8 0 0]);
bar(undercutoff,percBF(undercutoff),'facecolor',[0 0 1]);
plot([0 numel(datFiles)+1],cutoff*[1 1],'--k','linewidth',2);
set(gca,'ylim',[0 1],'ytick',0:0.2:1,'xlim',[0 numel(datFiles)+1],'xtick',1:numel(datFiles),...
    'xticklabel',datFiles,'fontsize',20);
xtickangle(45);
ylabel('Percentage bad p. fxns');
xlabel('Subject ID');


%% Save plots

if saveOn
   
    if nargin>2
        subdir = [varargin{2},'/'];
    else
        subdir = '';
    end
    
    if ~exist([dataDir,'figures/',subdir])
        mkdir([dataDir,'figures/',subdir]);
    end
     
    saveas(f1,[dataDir,'figures/',subdir,'within_PSE.svg']);
    saveas(f2,[dataDir,'figures/',subdir,'within_slope.svg']);
    saveas(f3,[dataDir,'figures/',subdir,'between_PSE.svg']);
    saveas(f4,[dataDir,'figures/',subdir,'sc_PSE.svg']);

end


end